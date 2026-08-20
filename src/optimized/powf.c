/*
 * Copyright (C) 2008-2020 Advanced Micro Devices, Inc. All rights reserved.
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
 * Contains implementation of powf()
 * Prototype :
 * float powf(float x, float y)
 *
 * Algorithm
 * x^y = e^(y*ln(x))
 *
 * Look in exp, log for the respective algorithms
 *

*/
#include <stdint.h>
#include <math.h>
#include <float.h>

#include "libm_macros.h"
#include "libm_util_amd.h"
#include <libm/typehelper.h>
#include <libm/amd_funcs_internal.h>
#include <libm/compiler.h>


#define N 10
#define TABLE_SIZE (1ULL << N)
#define MAX_POLYDEGREE  8

#if N == 8
#define POLY_DEGREE 4
extern const uint64_t log_256[];
extern const double log_f_inv_256[];
#define TAB_LOG(j)   asdouble(log_256[j])
#define TAB_F_INV log_f_inv_256
#define MANT_MASK_N  (0x000FF00000000000ULL)
#define MANT_MASK_N1 (0x0000080000000000ULL)
#elif N == 9
#define POLY_DEGREE 2
extern const double log_512[];
extern const double log_f_inv_512[];
#define TAB_LOG(j)   log_512[j]
#define TAB_F_INV log_f_inv_512
#define MANT_MASK_N  (0x000FFC0000000000ULL)
#define MANT_MASK_N1 (0x0000040000000000ULL)
#elif N == 10
#define POLY_DEGREE 2
extern const double log_1024[];
extern const double log_f_inv_1024[];
#define TAB_LOG(j)   log_1024[j]
#define TAB_F_INV log_f_inv_1024
#define MANT_MASK_N  (0x000FFC0000000000ULL)
#define MANT_MASK_N1 (0x0000020000000000ULL)
#endif

#define MANT_BITS_MASK (TABLE_SIZE - 1)
#define MANT1_BITS_MASK (1ULL << (N + 1))

#define EXPF_N 10
#define EXPF_POLY_DEGREE 3
#if EXPF_N == 6
#undef EXPF_POLY_DEGREE
#define EXPF_POLY_DEGREE 4
#elif EXPF_N == 5
#undef EXPF_POLY_DEGREE
#define EXPF_POLY_DEGREE 3
#elif EXPF_N == 4
#undef EXPF_POLY_DEGREE
#define EXPF_POLY_DEGREE 3
#endif

#define EXPF_TABLE_SIZE (1 << EXPF_N)
#define EXPF_MAX_POLY_DEGREE 4

static struct {
    double ALIGN(16) poly[MAX_POLYDEGREE];
    double_t ln2_lead, ln2_tail;
} log_data = {
#if N == 8
#endif
              //.ln2  = 0x1.62e42fefa39efp-1, /* ln(2) */
              .ln2_lead = 0x1.62e42e0000000p-1,
              .ln2_tail = 0x1.efa39ef35793cp-25,
              /*
               * Polynomial constants, 1/x! (reciprocal x)
               * To make better use of cache line,
               * we dont store 0! and 1!
               */
              .poly = { /* skip for 0/1 and 1/1 */
                       0x1.0000000000000p-1,    /* 1/2 */
                       0x1.5555555555555p-2,    /* 1/3 */
                       0x1.0000000000000p-2,    /* 1/4 */
                       0x1.999999999999ap-3,    /* 1/5 */
                       0x1.5555555555555p-3,    /* 1/6 */
                       0x1.2492492492492p-3,    /* 1/7 */
                       0x1.0000000000000p-3,    /* 1/8 */
                       0x1.c71c71c71c71cp-4,    /* 1/9 */
              },
};

#define C2 log_data.poly[0]
#define C3 log_data.poly[1]
#define C4 log_data.poly[2]
#define C5 log_data.poly[3]
#define C6 log_data.poly[4]
#define C7 log_data.poly[5]
#define C8 log_data.poly[6]
#define LN2_LEAD log_data.ln2_lead
#define LN2_TAIL log_data.ln2_tail

#include "expf_data.h"

extern double __two_to_jby128[];
extern double __two_to_jby256[];
extern double __two_to_jby512[];
extern double __two_to_jby1024[];

static struct expf_data expf_v2_data = {
#if EXPF_N == 10
    .ln2by_tblsz = 0x1.62e42fefa39efp-11,  /* ln(2)/1024 */
    .tblsz_byln2 = 0x1.71547652b82fep+10,  /* 1024/ln(2) */
#elif EXPF_N == 9
    .ln2by_tblsz = 0x1.62e42fefa39efp-10,  /* ln(2)/512 */
    .tblsz_byln2 = 0x1.71547652b82fep+9,   /* 512/ln(2) */
#elif EXPF_N == 8
    .ln2by_tblsz = 0x1.62e42fefa39efp-9,   /* ln(2)/256 */
    .tblsz_byln2 = 0x1.71547652b82fep+8,   /* 256/ln(2) */
#elif EXPF_N == 7
    .ln2by_tblsz = 0x1.62e42fefa39efp-8,   /* ln(2)/128 */
    .tblsz_byln2 = 0x1.71547652b82fep+7,   /* 128/ln(2) */
#else
    .ln2by_tblsz = 0x1.62e42fefa39efp-7,
    .tblsz_byln2 = 0x1.71547652b82fep+6,
#endif
    .Huge = 0x1.8000000000000p+52,
#if EXPF_N == 10
    .table_v3 = &__two_to_jby1024[0],
#elif EXPF_N == 9
    .table_v3 = &__two_to_jby512[0],
#elif EXPF_N == 8
    .table_v3 = &__two_to_jby256[0],
#elif EXPF_N == 7
    .table_v3 = &__two_to_jby128[0],
#elif EXPF_N == 6
    .table_v3 = &__two_to_jby64[0],
#elif EXPF_N == 5
    .table_v3 = (double*)L__two_to_jby32_table,
#endif

    .poly = {
        1.0,    /* 1/1! = 1 */
        0x1.0000000000000p-1,   /* 1/2! = 1/2    */
        0x1.5555555555555p-3,   /* 1/3! = 1/6    */
        0x1.cacccaa4ba57cp-5,   /* 1/4! = 1/24   */
    },
};

#define D1 expf_v2_data.poly[0]
#define D2 expf_v2_data.poly[1]
#define D3 expf_v2_data.poly[2]
#define D4 expf_v2_data.poly[3]

#define EXPF_LN2_BY_TBLSZ expf_v2_data.ln2by_tblsz
#define EXPF_TBLSZ_BY_LN2 expf_v2_data.tblsz_byln2
#define EXPF_HUGE expf_v2_data.Huge
#define EXPF_TABLE expf_v2_data.table

#define EXPF_FARG_MIN -0x1.9fe368p6     /* log(2^-150) ~= -103.97 */
#define EXPF_FARG_MAX 0x1.62e430p6     /* log(FLT_MAX) rounded up, ~88.7229 */
#define Ln2 0x1.62e42fefa39efp-1

struct log_table {
    double lead, tail;
};


/* Returns 0 if not int, 1 if odd int, 2 if even int.  The argument is
   the bit representation of a non-zero finite floating-point value.  */
static inline int
checkint (uint32_t iy)
{
    int32_t e = iy >> 23 & 0xff;
    if (e < 0x7f)
        return 0;
    if (e > 0x7f + 23)
        return 2;
    if (iy & (uint32_t)((1 << (0x7f + 23 - e)) - 1))
        return 0;
    if (iy & (uint32_t)(1 << (0x7f + 23 - e)))
        return 1;
    return 2;
}



static inline int
zeroinfnan (uint32_t ix)
{
    return 2 * ix - 1 >= 2u * POS_INF_F32 - 1;
}


/* Compute log(x) where x is in the float range (including subnormals).
 * Passing as double_t avoids float exponent field manipulation for subnormals:
 * (double)subnormal_float is exact and gives the correct biased exponent. */
static inline double_t
calculate_log(double_t x)
{
    double_t q, r;

    double_t dexpo, temp;

    /* Extract exponent and mantissa from the double representation.
     * For subnormal floats, the double holds the normalized form, so
     * the exponent and mantissa bits are correct without any adjustment.
     * The top 23 bits of the double mantissa match the float mantissa
     * for normal values, and give the correct normalized mantissa for
     * subnormals. */
    uint64_t udx = asuint64(x);

    int32_t expo = (int32_t)((udx >> EXPSHIFTBITS_DP64) &
                             (EXPBITS_DP64 >> EXPSHIFTBITS_DP64)) - EXPBIAS_DP64;

    uint32_t mant = (uint32_t)(udx >> (EXPSHIFTBITS_DP64 - EXPSHIFTBITS_SP32)) & MANTBITS_SP32;

    dexpo = (double)expo;

    uint32_t mant_n = mant & (MANTBITS_SP32 & ~((1u << (EXPSHIFTBITS_SP32 - N)) - 1u));

    /*
     * Step needed for better accuracy
    uint32_t mant_n1 = ux & 0x00004000;
    uint32_t j = (mant_n) + (mant_n1 << 1);
    */

    uint32_t j = (mant_n);

    mant |= HALFEXPBITS_SP32;               /* F */

    float j_times_half = asfloat(HALFEXPBITS_SP32 | j); /* Y */

    j >>= (EXPSHIFTBITS_SP32 - N);

    /* f = Y - F in double to avoid catastrophic cancellation */
    double_t f = (double_t)j_times_half - (double_t)asfloat(mant);

    r = f * TAB_F_INV[j];

    /* q = r + r^2*C2 */

    q = r + r * r * C2;

    /* m*log(2) + log(G) - poly */

    temp  = (dexpo * Ln2) + TAB_LOG(j);

    temp -= q;

    return temp;

}

static inline float calculate_exp(double_t x, uint64_t sign_bias)
{
    double_t poly, dn, r, z;
    uint64_t n, j;

    if (unlikely(x > EXPF_FARG_MAX)) {
        ALM_RAISE_FE_OVERFLOW();
        return asfloat((uint32_t)(sign_bias >> 32) | PINFBITPATT_SP32);
    }

    if (unlikely(x < EXPF_FARG_MIN)) {
        ALM_RAISE_FE_UNDERFLOW();
        return asfloat((uint32_t)(sign_bias >> 32));
    }

    z = x *  EXPF_TBLSZ_BY_LN2;

    /*
     * n  = (int) scale(x)
     * dn = (double) n
     */
#undef FAST_INTEGER_CONVERSION
#define FAST_INTEGER_CONVERSION 1
#if FAST_INTEGER_CONVERSION
    dn = z + EXPF_HUGE;

    n = asuint64(dn);

    dn -=  EXPF_HUGE;
#else
    n = z;

    dn = cast_i32_to_float(n);

#endif

    r = x - dn * EXPF_LN2_BY_TBLSZ;

    j = n % EXPF_TABLE_SIZE;

    /* polynomial = r + r^2*D2 + r^3*D3 */

    double_t tbl = asdouble(sign_bias | (asuint64(expf_v2_data.table_v3[j]) + (n << (52 - EXPF_N))));

#if EXPF_N >= 7
    poly = r + r * r * D2;              /* degree-2 sufficient for |r| <= ln2/128 or smaller */
#else
    poly = r + r * r * (D2 + r * D3);
#endif

    double_t result = tbl + tbl * poly;

    return (float_t)(result);

}

float ALM_PROTO_OPT(powf)(float x, float y)
{
    double_t logx, ylogx, result, q, r;
    uint32_t ux, uy;

    ux = asuint32(x);

    uy = asuint32(y);

    uint64_t sign_bias = 0;

    if (unlikely(zeroinfnan(uy))) {
        if (2 * uy == 0) {
            ALM_KEEP_ALIVE_SP32(x + x); // raise FE_INVALID if x is sNaN
            return 1.0f;
        }

        if (ux == ONEEXPBITS_SP32) {
            ALM_KEEP_ALIVE_SP32(y + y);  // raise FE_INVALID if y is sNaN
            return 1.0f;
        }

        if (2 * ux > 2u * POS_INF_F32 || 2 * uy > 2u * POS_INF_F32)
            return x + y;

        // y is +/-Inf
        if (ux == (ONEEXPBITS_SP32 | SIGNBIT_SP32))
            return 1.0f;

        if ((2 * ux < 2u * ONEEXPBITS_SP32) == !(uy & SIGNBIT_SP32))
            return 0.0f; /* |x|<1 && y==+inf, or |x|>1 && y==-inf */

        return asfloat(POS_INF_F32);
    } else if (unlikely(ux - 1 >= POS_INF_F32 - 1)) {
        // x is zero, negative, Inf or NaN

        if (unlikely(zeroinfnan(ux))) {
            // zero, NaN or Inf
            if (2 * ux > 2 * POS_INF_F32)
                return x + x;  // NaN: propagate qNaN, raise FE_INVALID for sNaN

            // x is +/-0 or +/-Inf
            // powf(+/-0, y): result is 0 for y>=0, Inf for y<0
            // powf(+/-Inf, y): result is Inf for y>=0, 0 for y<0
            // result negative for powf(-0, odd integer) and powf(-Inf, odd integer)
            return asfloat(((ux & SIGNBIT_SP32) && checkint(uy) == 1 ? SIGNBIT_SP32 : 0) |
                           // x is negative     && y is odd integer: negative else positive
                           (!(uy & SIGNBIT_SP32) == !(2*ux) ? 0 : POS_INF_F32)
                           //  y is nonnegative  == x is zero: zero else infinity
                );
        }

        // negative x
        int yint = checkint(uy);

        if (yint == 0)
            return (float)sqrt(x); /* x < 0, y non-integer: NaN + FE_INVALID */

        if (yint == 1)
            sign_bias = SIGNBIT_DP64;

        ux &= POS_BITSET_F32; /* x is negative, y is integer */
        x = asfloat(ux);
    }

    /* Near-1 path: x in (0.9375, 1.0625) = (0x3F700000, 0x3F880000).
     * Unsigned wrap trick: (0x3F880000 - ux) < (0x3F880000 - 0x3F700000)
     * is equivalent to ux in (0x3F700000, 0x3F880000). */
#define NEAR1_LO 0x3F700000u   /* 0.9375f */
#define NEAR1_HI 0x3F880000u   /* 1.0625f */

    if ((NEAR1_HI - ux) < (NEAR1_HI - NEAR1_LO)) {

        double dx = (double_t)x;

        double_t  u, u2, u3, u7;
        double_t  A1, A2, B1, B2, R1, R2;

        static const double ca[5] = {

                0x1.55555555554e6p-4, /* 1/2^2 * 3 */
                0x1.9999999bac6d4p-7, /* 1/2^4 * 5 */
                0x1.2492307f1519fp-9, /* 1/2^6 * 7 */
                0x1.c8034c85dfff0p-12 /* 1/2^8 * 9 */
        };

        /*
         * Less than threshold, no table lookup
         *
         */

        r = dx - 1.0;

        double_t u_by_2 = r / (2.0 + r);

        q = u_by_2 * r;  /* correction */

        u = u_by_2 + u_by_2;

#define CA1 ca[0]
#define CA2 ca[1]
#define CA3 ca[2]
#define CA4 ca[3]

        u2 = u * u;

        A1 = CA2 * u2 + CA1;

        A2 = CA4 * u2 + CA3;

        u3 = u2 * u;

        B1 = u3 * A1;

        u7 = u * (u3 * u3);

        B2 = u7 * A2;

        R1 = B1 + B2;

        R2 = R1 - q;

        logx = r + R2;

        ylogx = (double)y * logx;

        result = calculate_exp(ylogx, sign_bias);

        return (float)result;

    }


    logx = calculate_log((double_t)x);

    ylogx = (double)y * logx;

    result = calculate_exp(ylogx, sign_bias);

    return (float)result;
}
