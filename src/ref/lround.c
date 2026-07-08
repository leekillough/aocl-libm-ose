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
 * Overflow thresholds for lround(double), derived from LONG_MIN and LONG_MAX.
 *
 * 32-bit long (LONG_MAX == 0x7fffffff):
 *   LONG_MIN = -2^31 and LONG_MAX = 2^31-1 are both exactly representable as double.
 *   (double)LONG_MAX + 0.5 = 2^31 - 0.5  (exact; this is 0x1.fffffffep+30).
 *   (double)LONG_MIN - 0.5 = -(2^31+0.5) (exact; this is -0x1.00000001p+31).
 *   Both bounds are exactly representable, so strict > and < correctly exclude
 *   the overflowing boundary values.
 *
 * 64-bit long (LONG_MAX == 0x7fffffffffffffff):
 *   LONG_MIN = -2^63 is exactly representable; LONG_MAX = 2^63-1 rounds to 2^63.
 *   (double)LONG_MAX + 0.5 rounds to 2^63  (ULP at 2^63 is 2048, swamps 0.5).
 *   (double)LONG_MIN - 0.5 rounds to -2^63 (ULP at 2^63 swamps 0.5).
 *   LROUND_MIN = -2^63 = LONG_MIN as double, which is a valid input, so >= is
 *   required on the lower bound.  LROUND_MAX = 2^63 is itself overflowing, so
 *   < (strict) is required on the upper bound.
 *
 * NaN: ordered comparisons raise FE_INVALID for NaN per IEEE 754 / C Annex F;
 * __alm_handle_error raises it again for Inf and out-of-range finite values.
 * LONG_MIN is returned for all out-of-range inputs regardless of sign, matching
 * the x86 integer-indefinite value from cvtsd2si.
 *
 * Doubles with |x| >= 2^52 are already exact integers; adding 0.5 would land on
 * a halfway case that rounds to even, yielding the wrong result for odd integers.
 * Those values are handled by direct cast in the else-if branch.
 */
#if LONG_MAX == 0x7fffffff
#define LROUND_MIN    ((double)LONG_MIN - 0.5)    /* -(2^31+0.5), exact */
#define LROUND_MAX    ((double)LONG_MAX + 0.5)    /* 2^31-0.5, exact */
#define LROUND_INRANGE(x)  ((x) > LROUND_MIN && (x) < LROUND_MAX)
#else
#define LROUND_MIN    ((double)LONG_MIN - 0.5)    /* rounds to -2^63 = LONG_MIN */
#define LROUND_MAX    ((double)LONG_MAX + 0.5)    /* rounds to 2^63 (overflows) */
#define LROUND_INRANGE(x)  ((x) >= LROUND_MIN && (x) < LROUND_MAX)
#endif

/* 2^52 as a double bit-pattern: doubles with |x| >= 2^52 are exact integers. */
#define LROUND_INT_BITS   0x4330000000000000ULL

/* long is in the definition of the lround API, and is not chosen for its size */
long ALM_PROTO_REF(lround)(double x)
{
    UT64 u = { .f64 = x };

    long result = 0;

    if (unlikely(!LROUND_INRANGE(x))) {
        /* NaN, Inf, or x outside [LONG_MIN, LONG_MAX]: out of range.
         * LONG_MIN is returned for all such inputs regardless of sign. */
        __alm_handle_error(EXPBITS_DP64 | QNAN_MASK_64, AMD_F_INVALID);
        result = LONG_MIN;
    } else if (unlikely((u.u64 & POS_BITSET_DP64) >= LROUND_INT_BITS)) {
        /* |x| >= 2^52: already an exact integer; adding 0.5 would create a
         * halfway case that rounds to even, yielding a wrong result. */
        result = (long)x;
    } else {
        UT64 half = { .u64 = (u.u64 & SIGNBIT_DP64) | HALFEXPBITS_DP64 };
        result = (long)(x + half.f64);
    }

    return result;
}
