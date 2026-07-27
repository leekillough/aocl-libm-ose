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
#include <libm/amd_funcs_internal.h>
#include <libm/typehelper.h>
#include <fenv.h>
#include <limits.h>

/*
 * Overflow thresholds for llroundf(float), derived from LLONG_MIN and LLONG_MAX.
 *
 * LLONG_MIN = -2^63 is exactly representable as float; subtracting 0.5f rounds
 * back to -2^63 (ULP at 2^63 in float is 2^40, which swamps 0.5f).  So
 * LLROUNDF_MIN = -2^63 = LLONG_MIN as float, which is itself a valid input;
 * use >= on the lower bound.  LLONG_MAX = 2^63-1 rounds to 2^63 in float;
 * adding 0.5f stays 2^63.  So LLROUNDF_MAX = 2^63 overflows; use < (strict).
 *
 * NaN: ordered comparisons raise FE_INVALID for NaN per IEEE 754 / C Annex F;
 * __alm_handle_errorf raises it again for Inf and out-of-range finite values.
 * LLONG_MIN is returned for all out-of-range inputs regardless of sign.
 */
#define LLROUNDF_MIN        ((float)LLONG_MIN - 0.5f)   /* rounds to -2^63 = LLONG_MIN */
#define LLROUNDF_MAX        ((float)LLONG_MAX + 0.5f)   /* rounds to 2^63 (overflows) */
#define LLROUNDF_INRANGE(x) (((x) >= LLROUNDF_MIN) && ((x) < LLROUNDF_MAX))

/* long long is in the definition of the llroundf API, and is not chosen for its size */
long long ALM_PROTO_REF(llroundf)(float x)
{
    long long result = LLONG_MIN;

    if (unlikely(!LLROUNDF_INRANGE(x)))
    {
        /* NaN, Inf, or x outside [-2^63, 2^63): out of long long range. */
        feraiseexcept(FE_INVALID);
    }
    else
    {
        uint32_t ux = asuint32(x);
        if (likely((ux & POS_BITSET_F32) < EXP_VAL_23_F32)) {
            /* Only add if |x| < 2^23; if |x| >= 2^23: already an exact integer;
               adding 0.5f would create a halfway case that depends on rounding
               mode, yielding a wrong result. */
            x += asfloat((ux & SIGNBIT_SP32) | HALFEXPBITS_SP32);
        }
        result = (long long)x;
    }

    return result;
}
