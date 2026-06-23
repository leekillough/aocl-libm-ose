/*
 * Copyright (C) 2008-2024 Advanced Micro Devices, Inc. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*
 * ISO-IEC-10967-2: Elementary Numerical Functions
 * Signature:
 *   double cbrt(double x)
 *
 * Spec:
 * NOTE: The algorithm is optimized from the assembly version of CBRT function
 * To calculate (x)^1/3
 * step 1) Extract exponent and mantissa from input.
 * step 2) If it is subnormal input
 *         Exponent would be 0 and mantissa is non zero
 *         Normalise subnormal number:
 *         Shifting mantissa bits to left until MSB is 1 and 
 *         Number of times bits are shifted will contribute to exponent
 * step 3) Reduce the input [1, 2)
            3.1) Replace exponent with 3ff i.e 1
            3.2) Or with the mantissa
 * step 4) Scaling factor <= exponent/3 and Cuberoot2Index <= remainder (exponent % 3)
 * step 5) Polynomial approximation on reduced input
 * step 6) Multiply result of Polynomial approximation to cube-root remainder and scale factor
 * step 7) Return : Check for proper sign  and return the result
 *
 * Mathmatical Explanation
 * (x)^(1/3) = (x_d * 2^n)^(1/3)
 *           =  x_d^(1/3) * 2^(n/3)
 *           =  x_d^(1/3) * 2^(Quotient) * 2^(Remainder/3), where x_d is reduced input between [1,2)
 *
 *
*/

#include <stdint.h>
#include <libm_util_amd.h>
#include <libm/alm_special.h>
#include <immintrin.h>

#include <libm_macros.h>
#include <libm/types.h>

#include <libm/typehelper.h>
#include <libm/amd_funcs_internal.h>
#include <libm/compiler.h>
#include <libm/alm_special.h>
#include <cbrt_data.h>

#define ONE_BY_512          0.001953125                         // 0x3f60000000000000

#define CBRT_EXP_COEFF_1    3.33333333333333314829616256247E-1  // 0x3fd5555555555555
#define CBRT_EXP_COEFF_2    -1.11111111111111104943205418749E-1 // 0xbfbc71c71c71c71c
#define CBRT_EXP_COEFF_3    6.17283950617283916351141215273E-2  // 0x3faf9add3c0ca458
#define CBRT_EXP_COEFF_4    -4.11522633744855967363740489873E-2 // 0xbfa511e8d2b3183b
#define CBRT_EXP_COEFF_5    3.01783264746227734842687340233E-2  // 0x3f9ee7113506ac13
#define CBRT_EXP_COEFF_6    -2.34720317024843770636888251602E-2 // 0xbf98090d6221a247

/*
 * cbrt(2^k) high and low parts for k in {-2, -1, 0, 1, 2}, indexed by k+2.
 * Stored as two parallel arrays so both loads hit the same cache line and
 * the compiler can emit a single indexed load for each, with no branch.
 * rem from biased_exp % 3 is in {-2,-1,0,1,2}; index = rem + 2.
 */
static const double cbrt_rem_h[5] = {
    6.299605071544647216796875E-1,   /* cbrt(2^-2) high  0x3FE428A2F0000000  k=-2 */
    7.93700516223907470703125E-1,    /* cbrt(2^-1) high  0x3FE965FEA0000000  k=-1 */
    1.0E0,                           /* cbrt(2^0)  high  0x3FF0000000000000  k= 0 */
    1.259921014308929443359375E0,    /* cbrt(2^1)  high  0x3FF428A2F0000000  k= 1 */
    1.58740103244781494140625E0,     /* cbrt(2^2)  high  0x3FF965FEA0000000  k= 2 */
};

static const double cbrt_rem_t[5] = {
    1.77929718607039166806688400583E-8,  /* cbrt(2^-2) low  0x3e531ae515c447bb */
    9.76019226667272715610794680662E-9,  /* cbrt(2^-1) low  0x3e44f5b8f20ac166 */
    0.0E0,                               /* cbrt(2^0)  low  0x0000000000000000 */
    3.55859437214078333613376801167E-8,  /* cbrt(2^1)  low  0x3e631ae515c447bb */
    1.95203845333454543122158936132E-8,  /* cbrt(2^2)  low  0x3e54f5b8f20ac166 */
};

static inline void cbrt_special(double x, U32 code) {
    flt64_t ix = {.d = x};

    switch (code){
    case AMD_F_INVALID:
        __alm_handle_error(ix.u | QNAN_MASK_64, AMD_F_INVALID);
        break;
    case ALM_E_OVERFLOW:
        __alm_handle_error(ix.u, AMD_F_OVERFLOW);
        break;
    default:
        break;
    }
}

double
ALM_PROTO_OPT(cbrt)(double x) {
    uint64_t ix = asuint64(x);

    uint64_t ixe = EXPBITS_DP64 & ix;
    uint64_t ixm = MANTBITS_DP64 & ix;

    if (unlikely(ixe == PINFBITPATT_DP64)) {
        if (ixm == 0)
            cbrt_special(x, AMD_F_OVERFLOW);
        else
            cbrt_special(x, AMD_F_INVALID);
        return x + x;
    }

    ixe >>= EXPSHIFTBITS_DP64;

    if (unlikely(ixe == 0)) {
        if (ixm == 0)
            return 0.0;
        /* Subnormal: normalise by reinterpreting as 1.mantissa - 1.0 */
        uint64_t abs_ix = ix & POS_BITSET_DP64;
        double tmp = asdouble(abs_ix | ONEEXPBITS_DP64) - 1.0;
        ix = asuint64(tmp);
        ixe = (ix & EXPBITS_DP64) >> EXPSHIFTBITS_DP64;
        ixm = ix & MANTBITS_DP64;
        ixe += (uint64_t)EMIN_DP64;
    }

    const uint64_t sign = asuint64(x) >> 63;

    int64_t biased_exp = (int64_t)ixe - 1023;

    /* rem is in {-2,-1,0,1,2}; index into cbrt_rem_h/t is rem+2. */
    int64_t quotient = biased_exp / 3;
    int64_t rem      = biased_exp % 3;

    uint64_t exponDouble = (uint64_t)(quotient + 1023) << 52;

    /* Reduced mantissa in [0.5, 1). */
    double r = asdouble((ix & MANTBITS_DP64) | HALFEXPBITS_DP64);

    /*
     * 9-bit table index: upper 9 mantissa bits, rounded to nearest.
     * Bit 43 is the rounding bit; bits 44..52 are the index.
     */
    uint64_t mant_idx = ((ixm >> 43) & 1) + ((ixm >> 44) | 0x100);

    /*
     * Convert mant_idx to double without vcvtsi2sd.
     * mant_idx is in [256, 512].  OR it into the mantissa of 2^52 then
     * subtract the magic constant — IEEE 754 exact integer representability
     * guarantees the result equals mant_idx exactly.
     */
    static const double midx_magic = 4503599627370496.0; /* 2^52 */
    uint64_t midx_bits = mant_idx | (uint64_t)0x4330000000000000ULL;
    double midx_f = asdouble(midx_bits) - midx_magic;

    r = r - midx_f * ONE_BY_512;
    r = asdouble(InverseTable[mant_idx - 256]) * r;

    /*
     * Degree-6 polynomial: c1*r + c2*r^2 + c3*r^3 + c4*r^4 + c5*r^5 + c6*r^6.
     * Evaluated in two independent chains (pairs of consecutive terms) so both
     * FMA execution units on Zen 5 stay busy simultaneously.
     *
     * Pair grouping: p12 = r*(c1 + c2*r), p34 = r^3*(c3 + c4*r), p56 = r^5*(c5 + c6*r)
     * Chains A and B each compute two of these pairs independently, then sum.
     */
    double r2 = r * r;
    double r3 = r2 * r;
    double r4 = r2 * r2;
    double r5 = r4 * r;
    double r6 = r3 * r3;

    /* Two independent accumulators: A takes terms 1,3,5; B takes terms 2,4,6 */
    double polyA = CBRT_EXP_COEFF_1 * r;
    double polyB = CBRT_EXP_COEFF_2 * r2;
    polyA += CBRT_EXP_COEFF_3 * r3;
    polyB += CBRT_EXP_COEFF_4 * r4;
    polyA += CBRT_EXP_COEFF_5 * r5;
    polyB += CBRT_EXP_COEFF_6 * r6;
    double poly = polyA + polyB;

    /* cbrt_rem_h/t are indexed by rem+2, covering rem in {-2,-1,0,1,2}. */
    double cbrtRem_h = cbrt_rem_h[rem + 2];
    double cbrtRem_t = cbrt_rem_t[rem + 2];

    uint64_t fidx = (mant_idx - 256) << 1;
    double cbrtF_t = asdouble(F_H_L[fidx]);
    double cbrtF_h = asdouble(F_H_L[fidx + 1]);

    double bH = cbrtF_h * cbrtRem_h;
    double bT = cbrtF_t * cbrtRem_t + cbrtF_t * cbrtRem_h + cbrtRem_t * cbrtF_h;

    double ans = (poly * bT + bT) + (poly * bH + bH);
    ans = ans * asdouble(exponDouble);

    /* Apply sign via bit manipulation — no branch, no multiply by -1. */
    return asdouble(asuint64(ans) | (sign << 63));
}
