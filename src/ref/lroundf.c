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

#include "libm_util_amd.h"
#include <libm/alm_special.h>
#include <libm/amd_funcs_internal.h>
#include <limits.h>

/*
 * Overflow thresholds for lroundf(float), derived from LONG_MIN and LONG_MAX.
 *
 * For float inputs, ULP at the magnitude of LONG_MIN and LONG_MAX is far larger
 * than 0.5f in both the 32-bit and 64-bit long cases:
 *   32-bit long: ULP at 2^31 in float is 2^8 = 256;
 *     (float)LONG_MIN - 0.5f rounds to (float)LONG_MIN = -2^31.
 *     (float)LONG_MAX + 0.5f: (float)LONG_MAX rounds to 2^31; +0.5f stays 2^31.
 *   64-bit long: ULP at 2^63 in float is 2^40;
 *     (float)LONG_MIN - 0.5f rounds to (float)LONG_MIN = -2^63.
 *     (float)LONG_MAX + 0.5f rounds to 2^63.
 *
 * In both cases LROUNDF_MIN = (float)LONG_MIN (exactly representable as a power
 * of two and a valid lroundf input), so >= is required on the lower bound.
 * LROUNDF_MAX = 2^(N-1) is itself overflowing, so < (strict) is used on the upper.
 * No #if on LONG_MAX is needed: the same condition form applies to both sizes.
 *
 * NaN: ordered comparisons raise FE_INVALID for NaN per IEEE 754 / C Annex F;
 * __alm_handle_errorf raises it again for Inf and out-of-range finite values.
 * LONG_MIN is returned for all out-of-range inputs regardless of sign.
 *
 * Floats with |x| >= 2^23 are already exact integers; adding 0.5f would land on
 * a halfway case that rounds to even, yielding the wrong result for odd integers.
 * Those values are handled by direct cast in the else-if branch.
 */
#define LROUNDF_MIN  ((float)LONG_MIN - 0.5f)  /* rounds to (float)LONG_MIN */
#define LROUNDF_MAX  ((float)LONG_MAX + 0.5f)  /* rounds to 2^(N-1) */
#define LROUNDF_INRANGE(x) ((x) >= LROUNDF_MIN && (x) < LROUNDF_MAX)

/* 2^23 as a float bit-pattern: floats with |x| >= 2^23 are exact integers. */
#define LROUNDF_INT_BITS    0x4B000000U

/* long is in the definition of the lroundf API, and is not chosen for its size */
long ALM_PROTO_REF(lroundf)(float x)
{
    UT32 u = { .f32 = x };

    long result = 0;

    if (unlikely(!LROUNDF_INRANGE(x))) {
        /* NaN, Inf, or x outside [LONG_MIN, LONG_MAX]: out of range.
           LONG_MIN is returned for all such inputs regardless of sign. */
        __alm_handle_errorf(EXPBITS_SP32 | QNAN_MASK_32, AMD_F_INVALID);
        result = LONG_MIN;
    } else if (unlikely((u.u32 & POS_BITSET_F32) >= LROUNDF_INT_BITS)) {
        /* |x| >= 2^23: already an exact integer; adding 0.5f would create a
           halfway case that rounds to even, yielding a wrong result. */
        result = (long)x;
    } else {
        UT32 half = { .u32 = (u.u32 & SIGNBIT_SP32) | HALFEXPBITS_SP32 };
        result = (long)(x + half.f32);
    }

    return result;
}
