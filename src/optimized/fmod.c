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
 * Fast path (normal x and y with exponent difference d <= 52):
 *   Extract 53-bit integer significands Mx, My.  Compute rem = Mx * 2^d mod My
 *   using __uint128_t (at most 105 bits; GCC and Clang).  Pack rem back into
 *   a double.  No floating-point operations are performed, so no exceptions are
 *   raised.  x86_64 uses inline divq in Rem128 for a single-instruction modulo.
 *
 * General path (subnormals, or exponent difference > 52):
 *   FmodGeneral: iterative Veltkamp-Dekker reduction.  Wrapped with
 *   fetestexcept / feclearexcept to suppress spurious FE_INEXACT / FE_UNDERFLOW,
 *   since fmod is exact.
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

#define SCALE_2_POW_52    0x1p52
#define SCALE_2_POW_N52   0x1p-52

NOINLINE_COLD
static double FmodGeneral(double adx, double ady)
{
    /* General path: suppress spurious FE_INEXACT/FE_UNDERFLOW. */
    int to_clear = ~fetestexcept(FE_ALL_EXCEPT) & (FE_INEXACT | FE_UNDERFLOW);
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
        const double splitter = 0x1.0000002p+27;
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
    if (to_clear) {
        feclearexcept(to_clear);
    }
    return adx;
}

#ifdef __x86_64__

static inline uint64_t Rem128(uint64_t Mx, uint64_t My, int d)
{
    uint64_t quot;
    uint64_t rem;
    uint64_t hi = (d != 0) ? Mx >> (64 - d) : 0;
    uint64_t lo = Mx << d;
    __asm__("divq %[divisor]"
            : "=a"(quot), "=d"(rem)
            : "0"(lo), "1"(hi), [divisor] "r"(My));
    return rem;
}

#elif defined(__SIZEOF_INT128__)

static inline uint64_t Rem128(uint64_t Mx, uint64_t My, int d)
{
    return (uint64_t)(((__uint128_t)Mx << d) % My);
}

#endif

double ALM_PROTO_OPT(fmod)(double x, double y)
{
    uint64_t ay = asuint64(y) & ~SIGNBIT_DP64;
    double result = x;

    if (unlikely(ay > POS_INF_F64))
    {
        result = x * y;
    }
    else
    {
        uint64_t ax = asuint64(x) & ~SIGNBIT_DP64;

        if (unlikely((~ax & EXPBITS_DP64) == 0))
        {
            if (ax > POS_INF_F64)
            {
                result = x + x;
            }
            else
            {
                result = __alm_handle_error(INDEFBITPATT_DP64, AMD_F_INVALID);
            }
        }
        else if (unlikely(ay == 0))
        {
            result = __alm_handle_error(INDEFBITPATT_DP64, AMD_F_INVALID);
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
                uint64_t xe = ax >> EXPSHIFTBITS_DP64;
                uint64_t ye = ay >> EXPSHIFTBITS_DP64;
                int       d = (int)(xe - ye);

                if (unlikely((xe == 0) || (ye == 0) || (d > EXPSHIFTBITS_DP64)))
                {
                    // General path
                    adx = FmodGeneral(adx, ady);
                }
                else
                {
                    // Fast path

#if (defined(__GNUC__) || defined(__clang__)) && defined(__SIZEOF_INT128__)
                    // Integer fast path (GCC/Clang): no FP ops, no exceptions.
                    // rem = Mx * 2^d mod My; at most 53+52 = 105 bits.
                    uint64_t rem = Rem128((ax & MANTBITS_DP64) | IMPBIT_DP64,
                                          (ay & MANTBITS_DP64) | IMPBIT_DP64, d);
                    uint64_t rbits = 0;
                    if (rem != 0) {
                        int k = __builtin_clzll(rem) - (64 - MANTLENGTH_DP64);
                        rbits = ((int) ye > k) ?
                            ((ye - k) << EXPSHIFTBITS_DP64) | ((rem << k) & MANTBITS_DP64) :
                            rem << ((int)ye - 1);
                    }
                    adx = asdouble(rbits);
#else
                    // Fallback (non-GCC/Clang compiler, or no __uint128_t):
                    // single FP division with exception suppression.
                    // FmodGeneral is not needed since d <= 52 and both inputs are normal;
                    // one truncated division + FMA gives the exact remainder.
                    // FE_INEXACT comes from the division; FE_UNDERFLOW if the result
                    // is subnormal (possible when ye is near the minimum normal exponent).
                    int to_clear = ~fetestexcept(FE_ALL_EXCEPT) & (FE_INEXACT | FE_UNDERFLOW);
                    adx = fma(-(double)(uint64_t)(adx / ady), ady, adx);

                    // Division rounds up in FE_TONEAREST/FE_UPWARD; correct by one ady
                    if (adx < 0.0) {
                        adx += ady;
                    }
                    if (to_clear != 0) {
                        feclearexcept(to_clear);
                    }
#endif
                }

                result = copysign(adx, x);
            }
        }
    }

    return result;
}
