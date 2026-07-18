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
 * Integer fast path (exponent difference d <= 40):
 *   Extract 24-bit float significands Mx, My.  Compute rem = Mx * 2^d mod My
 *   using 64-bit arithmetic (Mx < 2^24 and d <= 40, so Mx * 2^d < 2^64).
 *   No floating-point operations, so no exceptions are raised.
 *   FE_UNDERFLOW is raised explicitly for subnormal results on Linux.
 *
 * Slow path (exponent difference d > 40):
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
#include <limits.h>

#include "libm_macros.h"
#include "libm_util_amd.h"
#include <libm/alm_special.h>
#include <libm/typehelper.h>
#include <libm/amd_funcs_internal.h>
#include <libm/compiler.h>

#define GCC_OR_CLANG    (defined(__GNUC__) || defined(__clang__))
#define FMODF_CHUNK_EXP 0x180000000000000u  /* 24 * 2^52 */

// Slow path: double-precision loop; suppress spurious FE_INEXACT.
// Used for d > 40, or compilers without __builtin_clz.
NOINLINE_COLD
static uint32_t FmodfGeneral(uint64_t ax, uint64_t ay, uint32_t xsign)
{
    int   except = fetestexcept(FE_ALL_EXCEPT);
    uint64_t quo = (ax - ay) / FMODF_CHUNK_EXP;
    double   adx = asdouble(ax);
    do
    {
        double ady = asdouble(quo * FMODF_CHUNK_EXP + ay);
        adx = fma(-(double)(uint64_t)(adx / ady), ady, adx);
        if (adx < 0.0)
        {
            adx += ady;
        }
    }
    while (unlikely(quo-- != 0));
    uint32_t result = asuint32((float) adx);

#ifdef __linux__
    if (unlikely((result-1 < POS_HDENORM_F32) &&
                 ((except & FE_UNDERFLOW) == 0))) {
        feraiseexcept(FE_UNDERFLOW);
    }
#endif

    if ((except & FE_INEXACT) == 0) {
        feclearexcept(FE_INEXACT);
    }

    return result | xsign;
}

#if GCC_OR_CLANG

typedef struct
{
    int e;
    uint32_t m;
} F32ExpMan;

static inline F32ExpMan F32Extract(uint32_t fax)
{
    int lz;
    return unlikely(fax < POS_LNORMAL_F32) ?
        lz = __builtin_clz(fax),
        (F32ExpMan) {
            .e = (int)(sizeof(uint32_t) * CHAR_BIT - MANTLENGTH_SP32 + 1) - lz,
            .m = fax << (lz - (int)(sizeof(uint32_t) * CHAR_BIT - MANTLENGTH_SP32))
        } :
        (F32ExpMan) {
            .e = (int)(fax >> EXPSHIFTBITS_SP32),
            .m = (fax & MANTBITS_SP32) | IMPBIT_SP32
        };
}

#endif // GCC_OR_CLANG

float ALM_PROTO_OPT(fmodf)(float x, float y)
{
    uint32_t   fax = asuint32(x);
    uint32_t   fay = asuint32(y) & POS_BITSET_F32;
    float   result = x;
    uint32_t xsign = fax & SIGNBIT_SP32;
    fax &= POS_BITSET_F32;

    if (unlikely(fay > POS_INF_F32))
    {   // |y| NaN
        result = x * y;
    }
    else if (unlikely(fax > POS_INF_F32))
    {   // |x| NaN
        result = x + x;
    }
    else if (unlikely(fax == POS_INF_F32) || (fay == 0))
    {   // |x| == Inf || y == 0
        result = __alm_handle_errorf(INDEFBITPATT_SP32, AMD_F_INVALID);
    }
    else if (fax >= fay)
    {   // |x| >= |y|
#if GCC_OR_CLANG
        F32ExpMan fpx = F32Extract(fax);
        F32ExpMan fpy = F32Extract(fay);
        int         d = fpx.e - fpy.e;

        if (unlikely(d > (int)(sizeof(uint64_t) * CHAR_BIT - MANTLENGTH_SP32)))
        {
#endif // GCC_OR_CLANG

            uint64_t ax = asuint64((double)x) & POS_BITSET_DP64;
            uint64_t ay = asuint64((double)y) & POS_BITSET_DP64;
            result = asfloat(FmodfGeneral(ax, ay, xsign));

#if GCC_OR_CLANG
        } else {
            uint32_t   rem = (uint32_t)(((uint64_t)fpx.m << d) % fpy.m);
            uint32_t rbits = 0;
            if (rem != 0) {
                int k = __builtin_clz(rem) + MANTLENGTH_SP32
                    - sizeof(rem) * CHAR_BIT;
                if (fpy.e > k)
                {
                    rbits = ((uint32_t)(fpy.e - k) << EXPSHIFTBITS_SP32)
                        | ((rem << k) & MANTBITS_SP32);
                } else {
                    // Subnormal float result
                    rbits = (fpy.e > 0) ? rem << (fpy.e - 1) : rem >> (1 - fpy.e);
#ifdef __linux__
                    feraiseexcept(FE_UNDERFLOW);
#endif
                }
            }
            result = asfloat(rbits | xsign);
        }
#endif // GCC_OR_CLANG

    }
    return result;
}
