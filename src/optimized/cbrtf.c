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
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*
 * ISO-IEC-10967-2: Elementary Numerical Functions
 * Signature:
 *   float cbrtf(float x)
 *
 * cbrt(x) = cbrt(m * 2^n)
 *         = cbrt(m) * 2^quotient * CbrtfRem[rem+2]
 * where m in [1,2), quotient = trunc(n/3), rem = n - 3*quotient in {-2..2}.
 *
 * cbrt(m) ~ (1 + t) * CubeRootTable[k]
 * where k = top 8 mantissa bits (0..255), t = 3-term poly on
 * r = m * DoubleReciprocalTable[k] - 1: t = r/3 - r^2/9 + 5*r^3/81.
 * All table lookups and arithmetic are done in double precision so that the
 * only rounding step is the final (float) cast, keeping the error at <=1 ULP
 * for all 2^32 inputs (18194 inputs are exactly 1 ULP from correctly rounded).
 */

#include <stdint.h>
#include <libm_util_amd.h>
#include <libm/alm_special.h>

#include <libm_macros.h>
#include <libm/types.h>

#include <libm/typehelper.h>
#include <libm/amd_funcs_internal.h>
#include <libm/compiler.h>
#include <cbrtf_data.h>

/*
 * cbrt(2^k) for k in {-2,-1,0,1,2}, indexed by k+2.
 * Stored in double so that their full precision reaches the final (float) cast.
 */
static const double CbrtfRem[5] = {
    0x1.428a2f98d728bp-1,   /* cbrt(2^-2)  k=-2 */
    0x1.965fea53d6e3dp-1,   /* cbrt(2^-1)  k=-1 */
    0x1.0000000000000p+0,   /* cbrt(2^0)   k= 0 */
    0x1.428a2f98d728bp+0,   /* cbrt(2^1)   k= 1 */
    0x1.965fea53d6e3dp+0,   /* cbrt(2^2)   k= 2 */
};

float
ALM_PROTO_OPT(cbrtf)(float x) {
    flt32_t  xdu    = { .f = x };
    uint32_t ix     = xdu.u;
    uint32_t ixe    = EXPBITS_SP32 & ix;
    uint32_t ixm    = MANTBITS_SP32 & ix;
    float    result = x;   /* cbrtf(x) -> x if x is +/-0, +/-Inf, qNaN */

    if (likely(ixe != EXPBITS_SP32)) {
        /* Not +/-Inf, NaN */
        ixe >>= EXPSHIFTBITS_SP32;
        if (likely((ixe | ixm) != 0)) {
            /* ixe is in [0, 254], so the conversion to int32_t is safe */
            int32_t biased_exp = (int32_t)ixe - 127;

            if (unlikely(ixe == 0)) {
                /* Subnormal: normalise via 1.mantissa - 1.0f self-subtraction trick. */
                flt32_t tmp = { .u = (ix & POS_BITSET_F32) | ONEEXPBITS_SP32 };
                tmp.f -= 1.0f;
                /* Extracted biased exponent is in [104, 126], within int32_t range. */
                biased_exp = (int32_t)((tmp.u & EXPBITS_SP32) >> EXPSHIFTBITS_SP32)
                             + (EMIN_SP32 - 127);
                ixm = tmp.u & MANTBITS_SP32;
            }

            /* The compiler strength-reduces the division to modular multiplication */
            int32_t quotient = biased_exp / 3;
            int32_t rem      = biased_exp - quotient * 3;

            /* Mantissa in [1, 2): set exponent field to 127. */
            flt32_t mfdu = { .u = ixm | ONEEXPBITS_SP32 };

            /* 8-bit table index: top 8 bits of the 23-bit mantissa. */
            uint32_t tidx = ixm >> 15;

            /*
             * All arithmetic in double so that the only rounding step is the final
             * (float) cast.  DoubleReciprocalTable and CubeRootTable hold 53-bit
             * accurate values; CbrtfRem is also double.  The conversion of mf
             * float->double is exact.
             */
            /* rd via FMA: product is exact internally, one rounding at the end. */
            double rd = fma((double)mfdu.f, DoubleReciprocalTable[tidx], -1.0);

            /*
             * 3-term poly: cbrt(1+r)-1 ~= r/3 - r^2/9 + 5*r^3/81.
             * Inner Horner step uses FMA to eliminate the intermediate rounding of
             * (1/3 + r*(-1/9)).  r^3 for the 3rd term is computed in parallel on
             * the r^2 chain; the correction adds one FMA to the critical path.
             */
            double r2 = rd * rd;
            double td = rd * fma(rd, -0x1.c71c71c71c71cp-4, 0x1.5555555555555p-2);
            td = fma(r2 * rd, 0x1.f9add3c0ca458p-5, td);

            flt32_t scaledu = { .u = (uint32_t)(quotient + 127) << 23 };
            double scale = CbrtfRem[rem + 2] * (double)scaledu.f;

            /* ans = (1+td)*cs = cs + td*cs; FMA avoids rounding the (1+td) sum. */
            double cs  = CubeRootTable[tidx] * scale;
            double ans = fma(td, cs, cs);

            result = copysignf((float)ans, x);
        }
    } else {
        /* +/-Inf, NaN */
        if ((ixm != 0) && ((ixm & QNAN_MASK_32) == 0)) {
            /* sNaN: quiet the NaN and raise FE_INVALID */
            result = __alm_handle_errorf(ix | QNAN_MASK_32, AMD_F_INVALID);
        }
    }

    return result;
}
