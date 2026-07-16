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

/******************************************
 * Implementation Notes:
 *
 * Prototype:
 * double fmod(double x, double y)
 *
 * Algorithm:
 * fmod(x, y) = x - n*y, where n = trunc(x/y).
 *
 * Fast path (normal x and y with exponent difference <= 52):
 *   n = trunc(|x|/|y|) fits in a uint64_t < 2^53, so n is exactly
 *   representable; the FMA then computes |x| - n*|y| with a single
 *   rounding rather than two.
 *
 * General path (subnormals, or exponent difference > 52):
 *   Reduce in 52-bit chunks: find the largest w = |y| * 2^(52*k) such
 *   that w <= |x|, do one Dekker-exact reduction step, divide w by 2^52,
 *   repeat until w == |y|. Isolated in a cold helper to keep the fast
 *   path's XMM register pressure low.
 *
 */

#if defined(__clang__) || defined(_MSC_VER)
#pragma STDC FENV_ACCESS ON
#endif

#include <stdint.h>
#include <math.h>
#include <fenv.h>
#include <float.h>

#include "libm_macros.h"
#include "libm_util_amd.h"
#include <libm/alm_special.h>
#include <libm/typehelper.h>
#include <libm/amd_funcs_internal.h>
#include <libm/compiler.h>

#define SCALE_2_POW_52    0x1p52   /* 2^52: scale up by one 52-bit chunk */
#define SCALE_2_POW_N52   0x1p-52  /* 2^-52: scale down by one 52-bit chunk */

/* General path for subnormals or exponent difference > 52.
 * Kept out-of-line and cold so the fast path saves no XMM registers. */
NOINLINE_COLD
static double FmodGeneral(double adx, const double ady)
{
    double w = ady;
    double t = adx * SCALE_2_POW_N52;
    while (w <= t)
    {
        w *= SCALE_2_POW_52;
    }
    for (;;)
    {
        double tw = w <= ady ? ady : w;
        double r = (double)(uint64_t)(adx / tw);
        /* Veltkamp split tw into hy + ty (each at most 27 significant bits).
         * splitter*tw can overflow for very large tw (around 2^997), so we
         * conservatively scale tw by 2^-28 before splitting and scale hy back
         * up. Scaling down by 2^-28 cannot underflow here because tw is
         * extremely large. */
        const double splitter = 0x1.0000002p+27; /* 2^27 + 1 */
        double hy;
        if (unlikely(tw > 0x1p996)) {
            double tw_sc = tw * 0x1p-28;
            double ctw = splitter * tw_sc;
            hy = (ctw - (ctw - tw_sc)) * 0x1p28;
        } else {
            double ctw = splitter * tw;
            hy = ctw - (ctw - tw);
        }
        double cr = splitter * r;
        double hr = cr - (cr - r);
        double ty = tw - hy;
        double tr = r - hr;
        double c = r*tw;
        /* Use fma throughout to avoid intermediate overflow: hy and hr are
         * Veltkamp high parts and can individually round up, so hy*hr may
         * exceed max_double (e.g. fmod(max_double, 1.0)) even though r*tw
         * <= adx. fma computes each product+accumulate in infinite precision,
         * preventing overflow and reducing rounding error. */
        double cc = fma(tr, ty, fma(ty, hr, fma(hy, tr, fma(hy, hr, -c))));
        double v = adx - c;
        double res = (((adx - v) - c) - cc) + v;
        adx = res < 0 ? res + tw : res;
        if (w <= ady)
        {
            break;
        }
        w *= SCALE_2_POW_N52;
    }
    return adx;
}

double ALM_PROTO_OPT(fmod)(double x, double y)
{
    uint64_t ay = asuint64(y) & ~SIGNBIT_DP64;
    double result = x;

    /* Check if y is NaN. If yes return NaN */
    if (unlikely(ay > POS_INF_F64))
    {
        result = x * y;
    }
    else
    {
        uint64_t ax = asuint64(x) & ~SIGNBIT_DP64;

        /* Check if x is NaN or INF */
        if (unlikely((~ax & EXPBITS_DP64) == 0))
        {
            /* x is NaN: return NaN. qNaN must not raise FE_INVALID. */
            if (ax > POS_INF_F64)
            {
                result = x + x;
            }
            else
            {
                /* x is INF. Return NaN and raise exception */
                result = __alm_handle_error(ay | QNANBITPATT_DP64, AMD_F_INVALID);
            }
        }
        /* Check if y is Zero. If yes, return NaN and raise exception */
        else if (unlikely(ay == 0))
        {
            result = __alm_handle_error(QNANBITPATT_DP64, AMD_F_INVALID);
        }
        else if (ax == ay)
        {
            result = copysign(0.0, x);
        }
        else
        {
            double adx = asdouble(ax);
            double ady = asdouble(ay);

            if (adx >= ady)
            {
                int  except = fetestexcept(FE_ALL_EXCEPT);
                uint64_t xe = ax >> EXPSHIFTBITS_DP64;
                uint64_t ye = ay >> EXPSHIFTBITS_DP64;

                if (unlikely((xe == 0) || (ye == 0) || (xe > (ye + EXPSHIFTBITS_DP64))))
                {
                    adx = FmodGeneral(adx, ady);
                }
                else
                {
                    /* Fast path: normal x and y with diff_exp <= 52.
                     * n = floor(|x|/|y|) < 2^53, so n is exactly representable; the FMA
                     * then computes |x| - n*|y| with a single rounding rather than two. */
                    adx = fma(-(double)(uint64_t)(adx / ady), ady, adx);
                    if (adx < 0)
                    {
                        // Division rounds up in FE_TONEAREST/FE_UPWARD; correct by one ady
                        adx += ady;
                    }
                }
                result = copysign(adx, x);

                // fmod is exact; suppress any spurious FE_INEXACT or FE_UNDERFLOW
                // from intermediate division and scaling operations.
                int to_clear = ~except & (FE_INEXACT | FE_UNDERFLOW);
                if (to_clear != 0) {
                    feclearexcept(to_clear);
                }
            }

        }
    }

    return result;
}
