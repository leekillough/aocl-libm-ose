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

/* Reduction chunk size: 2^24 steps per iteration, matching float mantissa width. */
#define FMODF_CHUNK_BITS    24
/* Scale step: multiply by 2^-24 to shrink w by one chunk per iteration. */
#define FMODF_CHUNK_DOWN    0x1p-24
/* IEEE 754 double bias, used when constructing the initial scale exponent. */
#define FMODF_DOUBLE_BIAS   1023
/* Double mantissa width in bits, used when packing the scale exponent. */
#define FMODF_MANTISSA_BITS 52

float ALM_PROTO_OPT(fmodf)(float x, float y)
{
    uint32_t fay = asuint32(y) & ~SIGNBIT_SP32;
    float result = 0.0f;

    /*Check if y in NaN. If yes, return NaN */
    if (unlikely(fay > POS_INF_F32))
    {
        result = x * y;
    }

    /* Check if y is Zero. If yes, return NaN and raise exception*/
    else if (unlikely(fay == 0))
    {
        result = __alm_handle_errorf(fay | QNANBITPATT_SP32, AMD_F_INVALID);
    }
    else
    {
        uint32_t fax = asuint32(x) & ~SIGNBIT_SP32;

        /* Check if x is NaN or INF */
        if (unlikely((fax & EXPBITS_SP32) >= EXPBITS_SP32))
        {
            /* x is NaN. Return NaN */
            if (fax > POS_INF_F32)
            {
                /*
                  The old Windows path that called __alm_handle_error for x
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

            double adx = asdouble(ax);
            double ady = asdouble(ay);

            if (adx < ady)
            {
                result = x;
            }
            else
            {
                uint64_t xe = (EXPBITS_DP64 & ax) >> FMODF_MANTISSA_BITS;
                uint64_t ye = (EXPBITS_DP64 & ay) >> FMODF_MANTISSA_BITS;

                int64_t scale = (int64_t)(FMODF_DOUBLE_BIAS) << FMODF_MANTISSA_BITS;
                int64_t quo = 0;

                if (ye < xe)
                {
                    quo = (int64_t)(xe - ye) / FMODF_CHUNK_BITS;
                    scale = (FMODF_CHUNK_BITS * quo + FMODF_DOUBLE_BIAS) << FMODF_MANTISSA_BITS;
                }

                double w = asdouble((uint64_t)scale) * ady;
                while (quo > 0)
                {
                    quo--;
                    adx -= (double)(uint64_t)(adx / w) * w;
                    w *= FMODF_CHUNK_DOWN;
                }

                adx -= (double)(uint64_t)(adx / w) * w;
                result = copysignf((float)adx, x);
            }
        }
    }

    return result;
}
