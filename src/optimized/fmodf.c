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
 * Integer fast path (exponent difference shift <= 40):
 *   Extract 24-bit float significands Mx, My.  Compute rem = Mx * 2^d mod My
 *   using 64-bit arithmetic (Mx < 2^24 and shift <= 40, so Mx * 2^d < 2^64).
 *   No floating-point operations, so no exceptions are raised.
 *   FE_UNDERFLOW is raised explicitly for subnormal results on Linux.
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

#if defined(__GNUC__) || defined(__clang__)

#define CLZ32(x) __builtin_clz(x)

#elif defined(_MSC_VER)

#include <intrin.h>
static inline int alm_clz32(uint32_t x)
{
    unsigned long idx;
    _BitScanReverse(&idx, x);
    return 31 - (int)idx;
}
#define CLZ32(x) alm_clz32(x)

#else

#error "Intrinsic for counting leading zeroes not found"

#endif

typedef struct
{
    int e;
    uint32_t m;
} F32ExpMan;

// Extract a 32-bit floating point into exponent and mantissa, handling subnormals
static inline F32ExpMan F32Extract(uint32_t fax)
{
    int lz;
    return unlikely(fax < POS_LNORMAL_F32) ?
        lz = CLZ32(fax),
        (F32ExpMan) {
            .e = (int)(sizeof(uint32_t) * CHAR_BIT - MANTLENGTH_SP32 + 1) - lz,
            .m = fax << (lz - (int)(sizeof(uint32_t) * CHAR_BIT - MANTLENGTH_SP32))
        } :
        (F32ExpMan) {
            .e = (int)(fax >> EXPSHIFTBITS_SP32),
            .m = (fax & MANTBITS_SP32) | IMPBIT_SP32
        };
}

float ALM_PROTO_OPT(fmodf)(float x, float y)
{
    uint32_t   fax = asuint32(x);
    uint32_t   fay = asuint32(y) & POS_BITSET_F32;
    float   result = x;
    uint32_t xsign = fax & SIGNBIT_SP32;
    fax &= POS_BITSET_F32;

    if (unlikely(((fay - 1) | fax) >= POS_INF_F32))
    {
        if (fay > POS_INF_F32)
        {   // |y| NaN
            result = x * y;
        }
        else if (fax > POS_INF_F32)
        {   // |x| NaN
            result = x + x;
        }
        else if ((fax == POS_INF_F32) || (fay == 0))
        {   // |x| == Inf || y == 0
            result = __alm_handle_errorf(INDEFBITPATT_SP32, AMD_F_INVALID);
        }
        else if (fax >= fay)
        {   // |x| >= |y|
            goto normal;
        }
    }
    else if (fax >= fay)
    {   // |x| >= |y|
    normal:
        F32ExpMan fpx = F32Extract(fax);
        F32ExpMan fpy = F32Extract(fay);
        int     shift = fpx.e - fpy.e;
        uint32_t  rem = fpx.m;
        const int maxshift = sizeof(uint64_t) * CHAR_BIT - MANTLENGTH_SP32;

        while (unlikely(shift > maxshift)) {
            rem = (uint32_t)(((uint64_t)rem << maxshift) % fpy.m);
            shift -= maxshift;
        }
        rem = (uint32_t)(((uint64_t)rem << shift) % fpy.m);
        if (likely(rem != 0)) {
            int k = CLZ32(rem) + MANTLENGTH_SP32
                - (int)(sizeof(rem) * CHAR_BIT);
            if (fpy.e > k)
            {
                rem = ((uint32_t)(fpy.e - k) << EXPSHIFTBITS_SP32)
                    | ((rem << k) & MANTBITS_SP32);
            } else {
                // Subnormal float result
                rem = (fpy.e > 0) ? rem << (fpy.e - 1) : rem >> (1 - fpy.e);
#ifdef __linux__
                feraiseexcept(FE_UNDERFLOW);
#endif
            }
        }
        result = asfloat(rem | xsign);
    }
    return result;
}
