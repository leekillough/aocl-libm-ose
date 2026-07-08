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
 * Overflow thresholds for llround(double), derived from LLONG_MIN and LLONG_MAX.
 *
 * LLONG_MIN = -2^63 is exactly representable as double; subtracting 0.5 rounds
 * back to -2^63 (ULP at 2^63 is 2048, which swamps 0.5).  So LLROUND_MIN =
 * -2^63 = LLONG_MIN as double, which is itself a valid input; use >= on the lower.
 * LLONG_MAX = 2^63-1 rounds to 2^63 in double; adding 0.5 stays 2^63.  So
 * LLROUND_MAX = 2^63 is itself overflowing; use < (strict) on the upper bound.
 *
 * NaN: ordered comparisons raise FE_INVALID for NaN per IEEE 754 / C Annex F;
 * __alm_handle_error raises it again for Inf and out-of-range finite values.
 * LLONG_MIN is returned for all out-of-range inputs regardless of sign.
 */
#define LLROUND_MIN       ((double)LLONG_MIN - 0.5)   /* rounds to -2^63 = LLONG_MIN */
#define LLROUND_MAX       ((double)LLONG_MAX + 0.5)   /* rounds to 2^63 (overflows) */
#define LLROUND_INRANGE(x) ((x) >= LLROUND_MIN && (x) < LLROUND_MAX)

/* 2^52 as a double bit-pattern: doubles with |x| >= 2^52 are already exact integers. */
#define LLROUND_INT_BITS   0x4330000000000000ULL
/* |0.5| in double, used to construct copysign(0.5, x). */
#define LLROUND_HALF_BITS    0x3FE0000000000000ULL

long long ALM_PROTO_REF(llround)(double x)
{
    UT64 u = { .f64 = x };
    long long result = 0;

    if (unlikely(!LLROUND_INRANGE(x))) {
        /* NaN, Inf, or x outside [-2^63, 2^63): out of long long range. */
        __alm_handle_error(EXPBITS_DP64 | QNAN_MASK_64, AMD_F_INVALID);
        result = LLONG_MIN;
    } else if (unlikely((u.u64 & POS_BITSET_DP64) >= LLROUND_INT_BITS)) {
        /* |x| >= 2^52: already an exact integer; adding 0.5 would create a
         * halfway case that rounds to even, yielding a wrong result. */
        result = (long long)x;
    } else {
        UT64 half = { .u64 = (u.u64 & SIGNBIT_DP64) | LLROUND_HALF_BITS };
        result = (long long)(x + half.f64);
    }

    return result;
}
