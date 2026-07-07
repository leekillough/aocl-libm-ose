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
#include <string.h>

/* Bounds for the valid llround(double) input range. */
#define LLROUND_D_MAX        0x1.0p+63          /* 2^63: first double too large for long long */
#define LLROUND_D_MIN       -0x1.0p+63          /* -2^63: long long minimum as double */
/* Bit pattern of 2^52: doubles with |x| >= 2^52 are already exact integers. */
#define LLROUND_D_INT_BITS   0x4330000000000000ULL
/* Bit pattern returned for out-of-range input: LLONG_MIN = 0x8000000000000000. */
#define LLROUND_OOR_BITS     0x8000000000000000ULL
/* Sign-bit mask and exponent of ±0.5 in double (used to construct signed half). */
#define LLROUND_SIGN_MASK    0x8000000000000000ULL
#define LLROUND_HALF_BITS    0x3FE0000000000000ULL  /* |0.5| in double */
/* Absolute-value mask for double. */
#define LLROUND_ABS_MASK     0x7fffffffffffffffULL

long long ALM_PROTO_REF(llround)(double x)
{
    uint64_t ui;
    memcpy(&ui, &x, sizeof(ui));
    long long result;

    if (unlikely(!((x >= LLROUND_D_MIN) && (x < LLROUND_D_MAX)))) {
        __alm_handle_error(LLROUND_OOR_BITS, AMD_F_NONE);
        result = (long long)LLROUND_OOR_BITS;
    } else if ((ui & LLROUND_ABS_MASK) >= LLROUND_D_INT_BITS) {
        /* |x| >= 2^52: already an exact integer; adding 0.5 would create a
         * halfway case that rounds to even, corrupting exact odd integers. */
        result = (long long)x;
    } else {
        ui = (ui & LLROUND_SIGN_MASK) | LLROUND_HALF_BITS;
        double half;
        memcpy(&half, &ui, sizeof(half));
        result = (long long)(x + half);
    }

    return result;
}
