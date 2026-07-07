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

// float remainderf(float x, float y)
//
// Returns x - n*y where n = rint(x/y) (nearest-even integer).
//
// Algorithm:
//
//   Single combined special-case check:
//     (ax - 1u) >= 0x7f7fffff || (ay - 1u) >= 0x7f7fffff
//   fires iff any of: x=NaN, x=Inf, y=NaN, y=Inf, y=0, x=0.
//   Normal finite non-zero inputs have ax,ay in [1, 0x7f7fffff], so
//   each of ax-1u, ay-1u is at most 0x7f7ffffe < 0x7f7fffff.
//   ay=0: ay-1u wraps to 0xffffffff >= 0x7f7fffff.
//   ax=0: ax-1u wraps to 0xffffffff >= 0x7f7fffff.
//   x/y=Inf (0x7f800000): ax/ay-1u = 0x7f7fffff >= 0x7f7fffff.
//   x/y=NaN (> 0x7f800000): ax/ay-1u >= 0x7f800000 > 0x7f7fffff.
//   Note: OR of two values each <= 0x7f7ffffe can reach 0x7f7fffff
//   (false positive for FLT_MAX pairs), so test each operand separately.
//
//   Work with ax = |x|, ay = |y|; apply sign of x at end.
//
//   Fast path (ax <= ay): n is 0 or 1, no division needed.
//     2*ax <= ay  -> n=0, result = x
//     2*ax >  ay  -> n=1, result = sign(x)*(ax - ay) [exact by Sterbenz]
//
//   Single-precision path (ax > ay, exponent diff <= ~24):
//     n = RneF(ax/ay): vroundss imm=8 (SSE4.1, rounding-mode independent) or
//     trunc(q+0.5) fallback.  r = ax - n*ay via vfnmadd231ss.
//     If r in [0, ay): result is exact, return immediately.
//     If r < 0 (n rounded up) or r >= ay (exponent diff > ~24): fall back to double.
//
//   Large-exponent path (d > 52): 24-bit chunk reduction, then single step.
//
// Special cases (IEEE 754):
//   remainder(x,   0)  -> NaN, invalid
//   remainder(Inf, y)  -> NaN, invalid  (even if y is a quiet NaN)
//   remainder(SNaN, y) -> NaN, invalid
//   remainder(x, SNaN) -> NaN, invalid
//   remainder(QNaN, y) -> QNaN          (quiet, no exception)
//   remainder(x, QNaN) -> QNaN          (quiet, no exception)
//   remainder(x, Inf)  -> x             (finite x)
//   remainder(0,   y)  -> 0             (sign preserved)

#include <stdint.h>
#include <math.h>
#include <string.h>

#include "libm_macros.h"
#include "libm_util_amd.h"
#include <libm/alm_special.h>
#include <libm/typehelper.h>
#include <libm/amd_funcs_internal.h>
#include <libm/compiler.h>

// RneF / RneD: round q to nearest integer, q nonnegative at every call site.
// SSE4.1 path: vroundss/vroundsd imm=8 -- nearest-even, rounding-mode independent.
// Non-SSE4.1 path: trunc(q+0.5) = round-half-away-from-zero (not nearest-even).
//   Half-integer ties that land on the wrong side are corrected at each call site
//   with: if (r + r == -ady && (int64_t)n & 1) r += ady;

#ifdef __SSE4_1__

// With SSE4.1: vroundss/vroundsd with imm=8 (nearest-even, suppress inexact,
// ignores MXCSR rounding-mode bits) -- one instruction, exact nearest-even.
#include <immintrin.h>
static inline float RneF(float q) {
    __m128 v = _mm_set_ss(q);
    return _mm_cvtss_f32(_mm_round_ss(v, v, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
}
static inline double RneD(double q) {
    __m128d v = _mm_set_sd(q);
    return _mm_cvtsd_f64(_mm_round_sd(v, v, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
}

#else // __SSE4_1__

// Without SSE4.1: trunc(q + 0.5) is identical to round(q) for nonnegative q
// (both are round-half-away-from-zero) but compiles to two instructions
// (add + trunc) instead of four, because the sign-handling overhead that
// round()/roundf() carries drops out when q is known to be nonnegative.
static inline float  RneF(float  q) { return truncf(q + 0.5f); }
static inline double RneD(double q) { return trunc(q + 0.5); }

#endif // __SSE4_1__

// memcpy() type-punning on float creates single register-register instruction,
// while union type-punning on float spills registers to the stack
static inline uint32_t FloatToUint(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }

float ALM_PROTO_OPT(remainderf)(float x, float y)
{
    uint32_t ix = FloatToUint(x);
    uint32_t iy = FloatToUint(y);
    uint32_t ax = ix & 0x7fffffffu;
    uint32_t ay = iy & 0x7fffffffu;

    // Single branch for all special cases.  For normal finite non-zero
    // inputs ax,ay in [1, 0x7f7fffff], so ax-1u and ay-1u are each at most
    // 0x7f7ffffe -- both strictly below 0x7f7fffffu.  The branch is taken
    // only when at least one of x, y is zero, Inf, or NaN.
    if (unlikely((ax - 1u) >= 0x7f7fffffu || (ay - 1u) >= 0x7f7fffffu)) {

        // Signaling NaN (SNaN): exponent=0xff, quiet bit (bit 22) clear,
        // mantissa non-zero.  Uint32 range [0x7f800001, 0x7fbfffff].
        // Per IEEE 754, any SNaN operand raises FE_INVALID.
        // Test: subtract 0x7f800001; result < 0x003fffff iff SNaN.
        if ((ax - 0x7f800001u) < 0x003fffffu ||
            (ay - 0x7f800001u) < 0x003fffffu)
            return __alm_handle_errorf(INDEFBITPATT_SP32, AMD_F_INVALID);

        // x=Inf is invalid regardless of y (even if y is a quiet NaN):
        // check before QNaN propagation so remainder(Inf, QNaN) -> FE_INVALID.
        if (ax == 0x7f800000u)
            return __alm_handle_errorf(INDEFBITPATT_SP32, AMD_F_INVALID);

        // Quiet NaN: propagate without raising an exception.
        if (ax > 0x7f800000u) return x;  // x=QNaN
        if (ay > 0x7f800000u) return y;  // y=QNaN

        // y=0: invalid (x is confirmed finite non-zero here).
        if (ay == 0u)
            return __alm_handle_errorf(INDEFBITPATT_SP32, AMD_F_INVALID);

        // y=Inf or x=0
        return x;
    }

    // Reconstruct |x| and |y| as floats.  IEEE 754 nonnegative floats are
    // monotone in their integer representation, so ax <= ay iff |x| <= |y|.
    float fax, fay;
    memcpy(&fax, &ax, 4);
    memcpy(&fay, &ay, 4);

    // Fast path: |x| <= |y|.  n is 0 or 1; no division needed.
    // 2*fax is exact (no overflow: ax < 0x7f800000).
    // fax - fay is exact by Sterbenz (fay/2 < fax <= fay for n=1 case).
    if (likely(ax <= ay)) {
        float ax2 = fax + fax;
        if (ax2 <= fay)
            return x;  // n=0 (includes tie 2|x|==|y|: rounds to 0)
        float r = fax - fay;  // exact by Sterbenz, n=1
        if (r == 0.0f)
            return copysignf(0.0f, x);
        return ix & 0x80000000u ? -r : r;
    }

    // |x| > |y|.  Attempt single-precision reduction first.
    // RneF rounds q to the nearest-even integer, rounding-mode independently,
    // via vroundss imm=8 (SSE4.1) or trunc(q+0.5) (fallback, nonnegative q).
    // The compiler emits vdivss + round + vfnmadd231ss in XMM registers.
    //
    // r_f can be negative when RneF rounds up (fractional part > 0.5, or a
    // half-integer tie rounded toward the odd neighbour in the fallback path),
    // and can be >= fay when the exponent difference exceeds ~24 bits.
    // Both cases fall through to the double path below.
    float q_f = fax / fay;             // vdivss
    float n_f = RneF(q_f);             // nearest-even (SSE4.1) or round-half-up (fallback)
    float r_f = fmaf(-n_f, fay, fax);  // vfnmadd231ss

    // Fast return: r_f is the exact remainder when it falls in [0, fay/2).
    // The correct remainder range is (-fay/2, fay/2]; values in [fay/2, fay)
    // mean n_f should have been one higher, so those fall through to the double
    // path.  Using the unsigned integer trick: for nonnegative IEEE 754 floats
    // the bit pattern is monotone, so FloatToUint(r_f) < FloatToUint(fay*0.5f)
    // iff 0.0f <= r_f < fay/2.  Negative r_f has its sign bit set, giving a
    // large unsigned value that always fails the test.
    uint32_t half_ay = FloatToUint(fay * 0.5f);
    if (likely(FloatToUint(r_f) < half_ay)) {
        if (r_f == 0.0f)
            return copysignf(0.0f, x);
        return ix & 0x80000000u ? -r_f : r_f;
    }

    // Float path inaccurate: r_f < 0 (RneF rounded up, or fay is a tiny
    // denormal causing q_f to overflow to inf), r_f in [fay/2, fay) (n_f too
    // small -- true remainder is negative), or r_f >= fay (exponent difference
    // > ~24 bits, n_f lost low bits).
    //
    // Recompute entirely in double using RneD (rounding-mode independent).
    double adx = (double)fax;
    double ady = (double)fay;
    uint64_t adx_bits = asuint64(adx);
    uint64_t ady_bits = asuint64(ady);
    int32_t xe_d = (int32_t)(adx_bits >> 52);
    int32_t ye_d = (int32_t)(ady_bits >> 52);
    int32_t d    = xe_d - ye_d;

    if (likely(d <= 52)) {
        double n_d = RneD(adx / ady);
        double r = fma(-n_d, ady, adx);
        if (unlikely(r >= ady))      r -= ady;
        else if (unlikely(r < -ady)) r += ady;

        // Fallback tie correction: trunc(q+0.5) gives n_d = N+1 (odd) for an
        // even-floor half-integer tie; RneD gives n_d = N (even).  Correct
        // the fallback case; harmless in the SSE4.1 path (n_d is already even).
        else if (unlikely(r + r == -ady && (int64_t)n_d & 1)) r += ady;
        if (r == 0.0)
            return copysignf(0.0f, x);
        float rf = (float)(r < 0.0 ? -r : r);
        uint32_t sign_r = r < 0.0 ? 0x80000000u : 0u;
        uint32_t result_bits = FloatToUint(rf) | ((ix & 0x80000000u) ^ sign_r);
        float ret;
        memcpy(&ret, &result_bits, 4);
        return ret;
    }

    // d > 52: multi-step 24-bit chunk reduction.
    {
        int32_t nsteps = d / 24;
        uint64_t wu = ady_bits + ((uint64_t)(24 * nsteps) << 52);
        double w = asdouble(wu);
        for (int32_t i = 0; i < nsteps; i++) {
            uint64_t q = (uint64_t)(adx / w);
            adx -= (double)q * w;
            w *= 0x1p-24;  // 2^-24
        }
        double n = RneD(adx / w);
        adx = fma(-n, w, adx);
        if (unlikely(adx >= ady))      adx -= ady;
        else if (unlikely(adx < -ady)) adx += ady;
        else if (unlikely(adx + adx == -ady && (int64_t)n & 1)) adx += ady;
        if (adx == 0.0)
            return copysignf(0.0f, x);
        float rf = (float)(adx < 0.0 ? -adx : adx);
        uint32_t sign_r = adx < 0.0 ? 0x80000000u : 0u;
        uint32_t result_bits = FloatToUint(rf) | ((ix & 0x80000000u) ^ sign_r);
        float ret;
        memcpy(&ret, &result_bits, 4);
        return ret;
    }
}
