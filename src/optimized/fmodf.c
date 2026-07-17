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
 * float fmodf(float x, float y)
 *
 * Algorithm:
 * fmodf(x, y) = x - n*y, where n = trunc(x/y).
 *
 * Integer fast path (both x and y normal, exponent difference d <= 104):
 *   Extract 24-bit float significands Mx, My.  Compute rem = Mx * 2^d mod My
 *   using __uint128_t (at most 128 bits).  Pack rem back into a float.
 *   No floating-point operations, so no exceptions are raised.
 *   FE_UNDERFLOW is raised explicitly for subnormal results on Linux.
 *
 * Slow path (subnormal inputs, or exponent difference d > 104):
 *   Double-precision iterative reduction.  Wrapped with fetestexcept /
 *   feclearexcept to suppress spurious FE_INEXACT, since fmodf is exact.
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

#define FMODF_CHUNK_EXP 0x180000000000000u  /* 24 * 2^52 */

float ALM_PROTO_OPT(fmodf)(float x, float y)
{
    uint32_t fay = asuint32(y) & ~SIGNBIT_SP32;
    float result = x;

    if (unlikely(fay > POS_INF_F32))
    {
        result = x * y;
    }
    else
    {
        uint32_t fax = asuint32(x) & ~SIGNBIT_SP32;

        if (unlikely((~fax & EXPBITS_SP32) == 0))
        {
            if (fax > POS_INF_F32)
            {
                result = x + x;
            }
            else
            {
                result = __alm_handle_errorf(fay | QNANBITPATT_SP32, AMD_F_INVALID);
            }
        }
        else if (unlikely(fay == 0))
        {
            result = __alm_handle_errorf(QNANBITPATT_SP32, AMD_F_INVALID);
        }
        else if (fax == fay)
        {
            result = copysignf(0.0f, x);
        }
        else
        {
            uint64_t ax = asuint64((double) x) & POS_BITSET_DP64;
            uint64_t ay = asuint64((double) y) & POS_BITSET_DP64;
            if (ax >= ay)
            {
                int xe_f = (int)(fax >> 23);   /* float biased exponent of |x| */
                int ye_f = (int)(fay >> 23);   /* float biased exponent of |y| */
                int d    = xe_f - ye_f;        /* exponent diff >= 0 since ax>=ay */

#if (defined(__GNUC__) || defined(__clang__)) && defined(__SIZEOF_INT128__)
                if (xe_f > 0 && ye_f > 0 && d <= 104)
                {
                    uint32_t    Mx = (fax & MANTBITS_SP32) | IMPBIT_SP32;
                    uint32_t    My = (fay & MANTBITS_SP32) | IMPBIT_SP32;
                    uint32_t   rem = (uint32_t)(((__uint128_t)Mx << d) % My);
                    uint32_t rbits = 0;
                    if (rem != 0) {
                        int k = __builtin_clz(rem) - 8;
                        if (ye_f > k)
                        {
                            rbits = ((uint32_t)(ye_f - k) << EXPSHIFTBITS_SP32)
                                | ((rem << k) & MANTBITS_SP32);
                        } else {
                            rbits = rem << (ye_f - 1); // Subnormal float result
#ifndef WINDOWS
                            feraiseexcept(FE_UNDERFLOW);
#endif
                        }
                    }
                    result = asfloat(rbits);
                }
                else
#endif

                {
                    // Slow path: double-precision loop; suppress spurious FE_INEXACT.
                    // Used for subnormals, d > 104, or targets without __uint128_t.
                    int except = fetestexcept(FE_ALL_EXCEPT);
                    uint64_t quo = (ax - ay) / FMODF_CHUNK_EXP;
                    double   adx = asdouble(ax);
                    do
                    {
                        double ady = asdouble(quo * FMODF_CHUNK_EXP + ay);
                        adx = fma(-(double)(uint64_t)(adx / ady), ady, adx);
                        if (adx < 0)
                        {
                            adx += ady;
                        }
                    }
                    while (unlikely(quo-- != 0));
                    result = (float) adx;

#ifndef WINDOWS
                    uint32_t fbits = asuint32(result);
                    if (unlikely((fbits < 0x00800000u) && (fbits != 0) &&
                                 ((except & FE_UNDERFLOW) == 0))) {
                        feraiseexcept(FE_UNDERFLOW);
                    }
#endif
                    if ((except & FE_INEXACT) == 0) {
                        feclearexcept(FE_INEXACT);
                    }
                }

                result = copysignf(result, x);
            }
        }
    }

    return result;
}
