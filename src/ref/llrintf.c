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

#include "fn_macros.h"
#include "libm_util_amd.h"
#include <libm/alm_special.h>
#include <libm/amd_funcs_internal.h>

long long ALM_PROTO_REF(llrintf)(float x)
{
#if (defined(__GNUC__) || defined(__clang__)) && defined(__SSE2__) && \
    (defined(__x86_64__) || defined(_M_X64))
    long long r;
    __asm__("cvtss2si %1, %0" : "=r"(r) : "x"(x));
    return r;
#else
    UT32 checkbits;
    UT32 val_2p23;
    long long result;

    checkbits.f32 = x;

    if ((checkbits.u32 & POS_BITSET_F32) > EXP_VAL_23_F32) {
        __alm_handle_errorf(SIGNBIT_DP64, AMD_F_INVALID);
        result = (long long)SIGNBIT_DP64;
    } else {
        val_2p23.u32 = (checkbits.u32 & SIGNBIT_SP32) | EXP_VAL_23_F32;
        result = (long long)((x + val_2p23.f32) - val_2p23.f32);
    }

    return result;
#endif
}
