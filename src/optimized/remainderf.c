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
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*
 * float remainderf(float x, float y)
 *
 * Returns x - n*y where n = rint(x/y) (nearest-even integer).
 *
 * Algorithm (matches Intel IMF z0 structure):
 *
 *   Single combined special-case check:
 *     (ax - 1u) >= 0x7f7fffff || (ay - 1u) >= 0x7f7fffff
 *   fires iff any of: x=NaN, x=Inf, y=NaN, y=Inf, y=0, x=0.
 *   Normal finite non-zero inputs have ax,ay in [1, 0x7f7fffff], so
 *   each of ax-1u, ay-1u is at most 0x7f7ffffe < 0x7f7fffff.
 *   ay=0: ay-1u wraps to 0xffffffff >= 0x7f7fffff.
 *   ax=0: ax-1u wraps to 0xffffffff >= 0x7f7fffff.
 *   x/y=Inf (0x7f800000): ax/ay-1u = 0x7f7fffff >= 0x7f7fffff.
 *   x/y=NaN (> 0x7f800000): ax/ay-1u >= 0x7f800000 > 0x7f7fffff.
 *   Note: OR of two values each <= 0x7f7ffffe can reach 0x7f7fffff
 *   (false positive for FLT_MAX pairs), so test each operand separately.
 *
 *   Work with ax = |x|, ay = |y|; apply sign of x at end.
 *
 *   Fast path (ax <= ay): n is 0 or 1, no division needed.
 *     2*ax <= ay  -> n=0, result = x
 *     2*ax >  ay  -> n=1, result = sign(x)*(ax - ay) [exact by Sterbenz]
 *
 *   Single-precision path (ax > ay, exponent diff <= ~24):
 *     n = rintf(ax/ay) via vdivss + vroundss nearest-even (XMM-only).
 *     r = ax - n*ay via vfnmadd231ss.
 *     If r in [0, ay): result is exact, return immediately.
 *     If r < 0 (rintf rounded up) or r >= ay (exponent diff > ~24,
 *     n_f lost low bits): fall back to double.
 *
 *   Large-exponent path (d > 52): 24-bit chunk reduction, then single step.
 *
 * Special cases (IEEE 754):
 *   remainder(x,   0) -> NaN, invalid
 *   remainder(Inf, y) -> NaN, invalid  (even if y is a quiet NaN)
 *   remainder(SNaN, y) -> NaN, invalid
 *   remainder(x, SNaN) -> NaN, invalid
 *   remainder(QNaN, y) -> QNaN         (quiet, no exception)
 *   remainder(x, QNaN) -> QNaN         (quiet, no exception)
 *   remainder(x, Inf) -> x             (finite x)
 *   remainder(0,   y) -> 0             (sign preserved)
 */

#include <stdint.h>
#include <math.h>
#include <string.h>

#include "libm_macros.h"
#include "libm_util_amd.h"
#include <libm/alm_special.h>
#include <libm/typehelper.h>
#include <libm/amd_funcs_internal.h>
#include <libm/compiler.h>

static inline uint32_t F2U(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }

float ALM_PROTO_OPT(remainderf)(float x, float y)
{
    uint32_t ix = F2U(x);
    uint32_t iy = F2U(y);
    uint32_t ax = ix & 0x7fffffffu;
    uint32_t ay = iy & 0x7fffffffu;

    /*
     * Single branch for all special cases.  For normal finite non-zero
     * inputs ax,ay in [1, 0x7f7fffff], so ax-1u and ay-1u are each at most
     * 0x7f7ffffe -- both strictly below 0x7f7fffffu.  The branch is taken
     * only when at least one of x, y is zero, Inf, or NaN.
     */
    if (unlikely((ax - 1u) >= 0x7f7fffffu || (ay - 1u) >= 0x7f7fffffu)) {
        /*
         * Signaling NaN (SNaN): exponent=0xff, quiet bit (bit 22) clear,
         * mantissa non-zero.  Uint32 range [0x7f800001, 0x7fbfffff].
         * Per IEEE 754, any SNaN operand raises FE_INVALID.
         * Test: subtract 0x7f800001; result < 0x003fffff iff SNaN.
         */
        if ((ax - 0x7f800001u) < 0x003fffffu ||
            (ay - 0x7f800001u) < 0x003fffffu)
            return __alm_handle_errorf(0xffc00000u, AMD_F_INVALID);
        /*
         * x=Inf is invalid regardless of y (even if y is a quiet NaN):
         * check before QNaN propagation so remainder(Inf, QNaN) -> FE_INVALID.
         */
        if (ax == 0x7f800000u)
            return __alm_handle_errorf(0xffc00000u, AMD_F_INVALID);
        /* Quiet NaN: propagate without raising an exception. */
        if (ax > 0x7f800000u) return x;              /* x=QNaN */
        if (ay > 0x7f800000u) return y;              /* y=QNaN */
        /* y=0: invalid (x is confirmed finite non-zero here). */
        if (ay == 0u)
            return __alm_handle_errorf(0xffc00000u, AMD_F_INVALID);
        return x;                                     /* y=Inf or x=0 */
    }

    /*
     * Reconstruct |x| and |y| as floats.  IEEE 754 positive floats are
     * monotone in their integer representation, so ax <= ay iff |x| <= |y|.
     */
    float fax, fay;
    memcpy(&fax, &ax, 4);
    memcpy(&fay, &ay, 4);

    /*
     * Fast path: |x| <= |y|.  n is 0 or 1; no division needed.
     * 2*fax is exact (no overflow: ax < 0x7f800000).
     * fax - fay is exact by Sterbenz (fay/2 < fax <= fay for n=1 case).
     */
    if (likely(ax <= ay)) {
        float ax2 = fax + fax;
        if (ax2 <= fay)
            return x;              /* n=0 (includes tie 2|x|==|y|: rounds to 0) */
        float r = fax - fay;       /* exact by Sterbenz, n=1 */
        if (r == 0.0f)
            return copysignf(0.0f, x);
        return (ix & 0x80000000u) ? -r : r;
    }

    /*
     * |x| > |y|.  Attempt single-precision reduction first.
     * rintf rounds the quotient to the nearest integer with ties-to-even,
     * which is exactly the n required by IEEE 754 remainder.  The compiler
     * emits vdivss + vroundss (nearest-even mode) + vfnmadd231ss, staying
     * entirely in XMM registers.
     *
     * r_f can be negative when rintf rounds up (q fractional part > 0.5),
     * and can be >= fay only when the exponent difference exceeds ~24 bits
     * (single precision runs out of mantissa bits for the quotient).
     * Both cases fall through to the double-precision path below.
     */
    float q_f = fax / fay;                   /* vdivss */
    float n_f = rintf(q_f);                  /* vroundss nearest-even */
    float r_f = fmaf(-n_f, fay, fax);        /* vfnmadd231ss */

    /*
     * Fast return: r_f is the exact remainder when it falls in (-fay, fay).
     * Using the unsigned integer trick: for positive IEEE 754 floats the
     * bit pattern is monotone, so F2U(r_f) < ay iff 0.0f <= r_f < fay.
     * Negative r_f has its sign bit set, so F2U gives a large value >= ay.
     */
    if (likely(F2U(r_f) < ay)) {
        if (r_f == 0.0f)
            return copysignf(0.0f, x);
        return (ix & 0x80000000u) ? -r_f : r_f;
    }

    /*
     * Float path inaccurate: r_f < 0 (rintf rounded up, or fay is a tiny
     * denormal causing q_f to overflow to inf) or r_f >= fay (exponent
     * difference > ~24 bits, n_f lost low bits).
     *
     * Recompute entirely in double, using rint for nearest-even rounding.
     */
    double adx = (double)fax;
    double ady = (double)fay;
    uint64_t adx_bits = asuint64(adx);
    uint64_t ady_bits = asuint64(ady);
    int32_t xe_d = (int32_t)(adx_bits >> 52);
    int32_t ye_d = (int32_t)(ady_bits >> 52);
    int32_t d    = xe_d - ye_d;

    if (likely(d <= 52)) {
        double n_d = rint(adx / ady);
        double r = fma(-n_d, ady, adx);
        /* r can be slightly outside (-ady, ady) due to double rounding. */
        if (unlikely(r >= ady))       r -= ady;
        else if (unlikely(r < -ady))  r += ady;
        if (r == 0.0)
            return copysignf(0.0f, x);
        float rf = (float)(r < 0.0 ? -r : r);
        uint32_t sign_r = (r < 0.0) ? 0x80000000u : 0u;
        uint32_t result_bits = F2U(rf) | ((ix & 0x80000000u) ^ sign_r);
        float ret;
        memcpy(&ret, &result_bits, 4);
        return ret;
    }

    /*
     * d > 52: multi-step 24-bit chunk reduction.
     */
    {
        int32_t nsteps = d / 24;
        double two_p_minus24 = asdouble(UINT64_C(0x3e70000000000000));  /* 2^-24 */
        uint64_t wu = ady_bits + ((uint64_t)(24 * nsteps) << 52);
        double w = asdouble(wu);
        for (int32_t i = 0; i < nsteps; i++) {
            uint64_t q = (uint64_t)(adx / w);
            adx -= (double)q * w;
            w *= two_p_minus24;
        }
        double n = rint(adx / w);
        adx = fma(-n, w, adx);
        if (adx == 0.0)
            return copysignf(0.0f, x);
        float rf = (float)(adx < 0.0 ? -adx : adx);
        uint32_t sign_r = (adx < 0.0) ? 0x80000000u : 0u;
        uint32_t result_bits = F2U(rf) | ((ix & 0x80000000u) ^ sign_r);
        float ret;
        memcpy(&ret, &result_bits, 4);
        return ret;
    }
}
