/*
 * Copyright (C) 2008-2024 Advanced Micro Devices, Inc. All rights reserved.
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
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*
 * ISO-IEC-10967-2: Elementary Numerical Functions
 * Signature:
 *   float cbrtf(float x)
 *
 * Rewritten to work entirely in single precision.  memcpy-based bit
 * reinterpretation helpers compile to vmovd on x86 with -O2/-O3 and are
 * portable to non-x86 targets.
 *
 * cbrt(x) = cbrt(m * 2^n)
 *         = cbrt(m) * 2^quotient * cbrtf_rem[rem+2]
 * where m in [1,2), quotient = trunc(n/3), rem = n - 3*quotient in {-2..2}.
 *
 * cbrt(m) ~ (1 + t) * FloatCubeRootTable[k]
 * where k = top 8 mantissa bits, t = Horner 2-term poly on r = m*Recip[k]-1.
 */

#include <stdint.h>
#include <string.h>
#include <libm_util_amd.h>
#include <libm/alm_special.h>

#include <libm_macros.h>
#include <libm/types.h>

#include <libm/typehelper.h>
#include <libm/amd_funcs_internal.h>
#include <libm/compiler.h>
#include <libm/alm_special.h>
#include <cbrtf_data.h>

static inline uint32_t F2U(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }
static inline float    U2F(uint32_t u) { float f; memcpy(&f, &u, 4); return f; }

/*
 * cbrt(2^k) for k in {-2,-1,0,1,2}, indexed by k+2.
 */
static const float cbrtf_rem[5] = {
    6.299605249474365823E-1f,   /* cbrt(2^-2)  k=-2 */
    7.937005259840997374E-1f,   /* cbrt(2^-1)  k=-1 */
    1.0f,                       /* cbrt(2^0)   k= 0 */
    1.2599210498948731648f,     /* cbrt(2^1)   k= 1 */
    1.5874010519681994748f,     /* cbrt(2^2)   k= 2 */
};

float
ALM_PROTO_OPT(cbrtf)(float x) {
    uint32_t ix  = F2U(x);
    uint32_t ixe = EXPBITS_SP32 & ix;
    uint32_t ixm = MANTBITS_SP32 & ix;

    if (unlikely(ixe == PINFBITPATT_SP32)) {
        if (ixm == 0)
            __alm_handle_errorf(ix, AMD_F_OVERFLOW);
        else
            __alm_handle_errorf(ix | QNAN_MASK_32, AMD_F_INVALID);
        return x + x;
    }

    ixe >>= EXPSHIFTBITS_SP32;

    if (unlikely(ixe == 0)) {
        if (ixm == 0)
            return x;
        /* Subnormal: normalise via 1.mantissa - 1.0f self-subtraction trick. */
        uint32_t tmp_u = (ix & POS_BITSET_F32) | ONEEXPBITS_SP32;
        tmp_u = F2U(U2F(tmp_u) - 1.0f);
        ixe = ((tmp_u & EXPBITS_SP32) >> EXPSHIFTBITS_SP32) + (uint32_t)EMIN_SP32;
        ixm = tmp_u & MANTBITS_SP32;
    }

    int32_t biased_exp = (int32_t)ixe - 127;

    /* Signed divide-by-3 via multiply-shift; single imulq + sar + sub. */
    int32_t quotient = (int32_t)(((int64_t)biased_exp * 0x55555556LL) >> 32) - (biased_exp >> 31);
    int32_t rem      = biased_exp - quotient * 3;

    /* Mantissa in [1, 2): set exponent field to 127. */
    float mf = U2F(ixm | ONEEXPBITS_SP32);

    /* 8-bit table index: top 8 bits of the 23-bit mantissa. */
    uint32_t tidx = ixm >> 15;

    /* Reciprocal-reduce: r = mf * Recip[tidx] - 1, |r| < 1/512. */
    float r = mf * FloatReciprocalTable[tidx] - 1.0f;

    /* Horner 2-term: cbrt(1+r) - 1 ~= r*(1/3 + r*(-1/9)). */
    float t = r * (1.0f/3.0f + r * (-1.0f/9.0f));

    /* Scale and remainder correction overlap the FP chain on the critical path. */
    float scale = cbrtf_rem[rem + 2] * U2F((uint32_t)(quotient + 127) << 23);

    /* Reconstruct mantissa cube-root. */
    float ans = (1.0f + t) * FloatCubeRootTable[tidx];

    return copysignf(ans * scale, x);
}
