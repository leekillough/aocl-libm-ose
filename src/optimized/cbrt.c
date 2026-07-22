/*
 * Copyright (C) 2008-2026 Advanced Micro Devices, Inc. All rights reserved.
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
 *         Reinterpret the absolute value as 1.mantissa in [1,2) then
 *         subtract 1.0 to isolate the mantissa; the resulting exponent
 *         field gives the shift count, from which the true exponent is
 *         derived.
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
#include <cbrt_data.h>

/* ONE_BY_512 = 2^-9 exactly; all others are nearest doubles to Taylor coefficients. */
#define ONE_BY_512          0x1.0000000000000p-9

#define CBRT_EXP_COEFF_1    0x1.5555555555555p-2   /*  1/3     */
#define CBRT_EXP_COEFF_2   -0x1.c71c71c71c71cp-4   /* -1/9     */
#define CBRT_EXP_COEFF_3    0x1.f9add3c0ca458p-5   /*  5/81    */
#define CBRT_EXP_COEFF_4   -0x1.511e8d2b3183bp-5   /* -10/243  */
#define CBRT_EXP_COEFF_5    0x1.ee7113506ac13p-6   /*  22/729  */
#define CBRT_EXP_COEFF_6   -0x1.8090d6221a247p-6   /* -154/6561 */

/*
 * cbrt(2^k) high and low parts for k in {-2,-1,0,1,2}, indexed by k+2.
 * Table lookup converts control dependence (branch on rem) into data
 * dependence (indexed load), which the out-of-order core hides better
 * than a branch when rem is uniformly random over {-2,-1,0,1,2}.
 */
static const double CbrtRemH[5] = {
    0x1.428a2f0000000p-1,  /* cbrt(2^-2) high  k=-2 */
    0x1.965fea0000000p-1,  /* cbrt(2^-1) high  k=-1 */
    0x1.0000000000000p+0,  /* cbrt(2^0)  high  k= 0 */
    0x1.428a2f0000000p+0,  /* cbrt(2^1)  high  k= 1 */
    0x1.965fea0000000p+0,  /* cbrt(2^2)  high  k= 2 */
};

static const double CbrtRemT[5] = {
    0x1.31ae515c447bbp-26, /* cbrt(2^-2) low   k=-2 */
    0x1.4f5b8f20ac166p-27, /* cbrt(2^-1) low   k=-1 */
    0x0.0000000000000p+0,  /* cbrt(2^0)  low   k= 0 */
    0x1.31ae515c447bbp-25, /* cbrt(2^1)  low   k= 1 */
    0x1.4f5b8f20ac166p-26, /* cbrt(2^2)  low   k= 2 */
};

double
ALM_PROTO_OPT(cbrt)(double x) {
    flt64_t  xdu = { .d = x };
    uint64_t ix  = xdu.u;
    uint64_t ixe = EXPBITS_DP64 & ix;
    uint64_t ixm = MANTBITS_DP64 & ix;

    /* +/-Inf and NaN: early exit, flat path for normals below. */
    if (unlikely(ixe == EXPBITS_DP64)) {
        if ((ixm != 0) && ((ixm & QNAN_MASK_64) == 0))
            return __alm_handle_error(ix | QNAN_MASK_64, AMD_F_INVALID);
        return x;  /* +/-Inf -> +/-Inf, qNaN -> qNaN */
    }

    ixe >>= EXPSHIFTBITS_DP64;

    /* +/-0 */
    if (unlikely((ixe | ixm) == 0))
        return x;

    int64_t biased_exp = (int64_t)ixe - 1023;

    if (unlikely(ixe == 0)) {
        /* Subnormal: reinterpret as 1.mantissa - 1.0 to extract true exponent. */
        flt64_t tmp = { .u = (ix & POS_BITSET_DP64) | ONEEXPBITS_DP64 };
        tmp.d -= 1.0;
        biased_exp = (int64_t)((tmp.u & EXPBITS_DP64) >> EXPSHIFTBITS_DP64)
                     + (EMIN_DP64 - 1023);
        ixm = tmp.u & MANTBITS_DP64;
    }

    /* The compiler strength-reduces the division to modular multiplication. */
    int64_t quotient = biased_exp / 3;
    int64_t rem      = biased_exp - quotient * 3;

    /* Reduced mantissa in [0.5, 1). */
    flt64_t rdu = { .u = ixm | HALFEXPBITS_DP64 };

    /*
     * 9-bit table index: upper 9 mantissa bits, rounded to nearest.
     * Bit 43 is the rounding bit; bits 44..52 are the unrounded index.
     */
    uint64_t mant_idx = ((ixm >> 43) & 1) + ((ixm >> 44) | 0x100);

    /*
     * Convert mant_idx to double without vcvtsi2sd: OR into the mantissa
     * of 2^52 then subtract the magic constant.  IEEE 754 exact integer
     * representability guarantees the result equals mant_idx exactly.
     */
    flt64_t midx = { .u = mant_idx | EXP_VAL_52_DP64 };
    flt64_t mant = { .u = InverseTable[mant_idx - 256] };

    double idx_frac = (midx.d - 0x1p52) * ONE_BY_512;
    double r = mant.d * (rdu.d - idx_frac);

    double r2 = r * r;
    double r3 = r2 * r;
    double r4 = r2 * r2;
    double r5 = r4 * r;
    double r6 = r3 * r3;

    /*
     * Two independent accumulator chains (odd A, even B) keep both FP
     * execution units on Zen 5 busy simultaneously.
     */
    double polyA = CBRT_EXP_COEFF_1 * r;
    double polyB = CBRT_EXP_COEFF_2 * r2;
    polyA += CBRT_EXP_COEFF_3 * r3;
    polyB += CBRT_EXP_COEFF_4 * r4;
    polyA += CBRT_EXP_COEFF_5 * r5;
    polyB += CBRT_EXP_COEFF_6 * r6;
    double poly = polyA + polyB;

    double cbrtRem_h = CbrtRemH[rem + 2];
    double cbrtRem_t = CbrtRemT[rem + 2];

    uint64_t fidx = (mant_idx - 256) << 1;
    flt64_t cbrtF_t = { .u = F_H_L[fidx] };
    flt64_t cbrtF_h = { .u = F_H_L[fidx + 1] };

    double bH = cbrtF_h.d * cbrtRem_h;
    double bT = (cbrtF_t.d * cbrtRem_t) + (cbrtF_t.d * cbrtRem_h) + (cbrtRem_t * cbrtF_h.d);

    /* Two independent (poly*b + b) pairs keep both FP units busy at the final stage. */
    double ans = (poly * bT + bT) + (poly * bH + bH);

    /* Scale by 2^quotient via integer-shift union load; no domain crossing. */
    flt64_t scale = { .u = (uint64_t)(quotient + 1023) << 52 };
    return copysign(ans * scale.d, x);
}
