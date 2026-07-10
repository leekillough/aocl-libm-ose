/*
 * Copyright (C) 2008-2025 Advanced Micro Devices, Inc. All rights reserved.
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

#include <stdint.h>
#include <math.h>
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

    /* Check if y is NaN. If yes, return NaN */
    if (unlikely(fay > POS_INF_F32))
    {
        result = x * y;
    }
    /* Check if y is Zero. If yes, return NaN and raise exception*/
    else if (unlikely(fay == 0))
    {
        result = __alm_handle_errorf(QNANBITPATT_SP32, AMD_F_INVALID);
    }
    else
    {
        uint32_t fax = asuint32(x) & ~SIGNBIT_SP32;

        /* Check if x is NaN or INF */
        if (unlikely((~fax & EXPBITS_SP32) == 0))
        {
            /* x is NaN. Return NaN */
            if (fax > POS_INF_F32)
            {
                /*
                  The old Windows path that called __alm_handle_errorf for x
                  NaN was wrong; it unconditionally raised FE_INVALID for qNaN
                  inputs, which violates IEEE 754.
                */
                result = x + x;
            }
            else
            {
                /* x is INF. Return NaN and raise exception */
                result = __alm_handle_errorf(fay | QNANBITPATT_SP32, AMD_F_INVALID);
            }
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
                /* quo = floor(log2( |x|/|y| ) / 24 ) */
                uint64_t quo = (ax - ay) / FMODF_CHUNK_EXP;
                double   adx = asdouble(ax);
                do
                {
                    /* |y| * 2^(24*quo) */
                    double ady = asdouble(quo * FMODF_CHUNK_EXP + ay);

                    /* Subtract floor( |x|/|y| ) * |y| from |x| */
                    adx = fma(-(double)(uint64_t)(adx / ady), ady, adx);

                    /* Division rounds up in FE_TONEAREST/FE_UPWARD; correct by one ady */
                    if (adx < 0)
                        adx += ady;
                }
                while (unlikely(quo-- != 0));

                /* Convert reduced |x| to float and copy original x sign */
                result = copysignf((float) adx, x);
            }
        }
    }

    return result;
}
