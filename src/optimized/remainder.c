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

// double remainder(double x, double y)
//
// Returns x - n*y where n = roundTiesToEven(x/y) (nearest-even integer).
//
// Algorithm:
//
//   Single combined special-case check:
//     (ax - 1) >= 0x7fefffffffffffff || (ay - 1) >= 0x7fefffffffffffff
//   fires iff any of: x=NaN, x=Inf, y=NaN, y=Inf, y=0, x=0.
//   Normal finite non-zero inputs have ax,ay in [1, 0x7fefffffffffffff], so
//   each of ax-1, ay-1 is at most 0x7feffffffffffffe < 0x7fefffffffffffff.
//   ay=0: ay-1 wraps to 0xffffffffffffffff >= 0x7fefffffffffffff.
//   ax=0: ax-1 wraps to 0xffffffffffffffff >= 0x7fefffffffffffff.
//   ax=Inf (0x7ff0000000000000): ax-1 = 0x7fefffffffffffff.
//   ax=NaN (> 0x7ff0000000000000): ax-1 >= 0x7ff0000000000000.
//
//   Work with ax = |x|, ay = |y|; apply sign of x at end.
//
//   Fast path (ax <= ay): n is 0 or 1, no division needed.
//     2*ax <= ay  -> n=0, result = x
//     2*ax >  ay  -> n=1, result = sign(x)*(ax - ay) [exact by Sterbenz]
//
//   Single-step path (exponent diff <= 52):
//     n = RneD(ax/ay): vroundsd imm=8 (SSE4.1) or trunc(q+0.5) fallback.
//     r = ax - n*ay via vfnmadd231sd.
//     Boundary correction for FP rounding overshoot, plus fallback tie correction.
//
//   Large-exponent path (d > 52): 24-bit chunk reduction, then single step.
//
// Special cases (IEEE 754):
//   remainder(x,   0) -> NaN, invalid
//   remainder(Inf, y) -> NaN, invalid  (when y is not NaN)
//   remainder(SNaN, y) -> NaN, invalid
//   remainder(x, SNaN) -> NaN, invalid
//   remainder(QNaN, y) -> QNaN         (quiet, no exception)
//   remainder(x, QNaN) -> QNaN         (quiet, no exception)
//   remainder(x, Inf) -> x             (finite x)
//   remainder(0,   y) -> 0             (sign preserved)

#include <stdint.h>
#include <math.h>
#include <float.h>

#include "libm_macros.h"
#include "libm_util_amd.h"
#include <libm/alm_special.h>
#include <libm/typehelper.h>
#include <libm/amd_funcs_internal.h>
#include <libm/compiler.h>

// RneD: round q to nearest integer, ties to even, independent of the current
// FP rounding mode.  q is always nonnegative at every call site here.

#ifdef __SSE4_1__

// With SSE4.1: vroundsd imm=8 -- one instruction, exact nearest-even.
#include <immintrin.h>
static inline double RneD(double q) {
    __m128d v = _mm_set_sd(q);
    return _mm_cvtsd_f64(_mm_round_sd(v, v, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
}

#else // __SSE4_1__

// Without SSE4.1: trunc(q + 0.5) for nonneg q equals round-half-away-from-zero.
// Called only with q in [0, 2^52), so q + 0.5 is always exact (no rounding from the
// addition), making this independent of the current FP rounding mode.
//
// Differs from nearest-even only at half-integer ties; corrected at each call site with:
//   if (r + r == -y && (int64_t)n & 1) r += y;
// (no-op in the SSE4.1 path where nearest-even always gives an even n).
static inline double RneD(double q) { return trunc(q + 0.5); }

#endif // __SSE4_1__

double ALM_PROTO_OPT(remainder)(double x, double y)
{
    uint64_t ix = asuint64(x);
    uint64_t iy = asuint64(y);
    uint64_t ax = ix & UINT64_C(0x7fffffffffffffff);
    uint64_t ay = iy & UINT64_C(0x7fffffffffffffff);
    double result;

    // Single branch for all special cases.  For normal finite non-zero
    // inputs ax,ay in [1, 0x7fefffffffffffff], so ax-1 and ay-1 are each at
    // most 0x7feffffffffffffe -- both strictly below 0x7fefffffffffffff.
    // The branch is taken only when at least one of x, y is zero, Inf, or NaN.
    if (unlikely((ax - 1) >= UINT64_C(0x7fefffffffffffff) ||
                 (ay - 1) >= UINT64_C(0x7fefffffffffffff))) {

        // Signaling NaN (SNaN): exponent=0x7ff, quiet bit (bit 51) clear,
        // mantissa non-zero. Uint64 [0x7ff0000000000001, 0x7ff7ffffffffffff].
        // Per IEEE 754, any SNaN operand raises FE_INVALID.
        // Signaling NaN: quiet the operand's own bits to preserve the payload.
        if ((ax - UINT64_C(0x7ff0000000000001)) < UINT64_C(0x0007ffffffffffff)) {
            result = __alm_handle_error(ix | UINT64_C(0x0008000000000000), AMD_F_INVALID);
        } else if ((ay - UINT64_C(0x7ff0000000000001)) < UINT64_C(0x0007ffffffffffff)) {
            result = __alm_handle_error(iy | UINT64_C(0x0008000000000000), AMD_F_INVALID);

        // Quiet NaN: propagate without raising an exception.
        } else if (ax > UINT64_C(0x7ff0000000000000)) {
            result = x;   // x=QNaN
        } else if (ay > UINT64_C(0x7ff0000000000000)) {
            result = y;   // y=QNaN

        // x=Inf is invalid (y is confirmed not NaN here).
        } else if (ax == UINT64_C(0x7ff0000000000000)) {
            result = __alm_handle_error(QNANBITPATT_DP64, AMD_F_INVALID);

        // y=0: invalid (x is confirmed finite non-zero here).
        } else if (ay == 0) {
            result = __alm_handle_error(QNANBITPATT_DP64, AMD_F_INVALID);

        } else {
            result = x;   // y=Inf or x=0
        }

    } else {
        // Reconstruct |x| and |y| as doubles.  IEEE 754 nonnegative doubles are
        // monotone in their integer representation, so ax <= ay iff |x| <= |y|.
        double adx = asdouble(ax);
        double ady = asdouble(ay);

        // Fast path: |x| <= |y|.  n is 0 or 1; no division needed.
        // 2*adx is exact (no overflow: ax < 0x7ff0000000000000).
        // adx - ady is exact by Sterbenz (ady/2 < adx <= ady for n=1 case).
        if (likely(ax <= ay)) {
            double ax2 = adx + adx;
            if (ax2 <= ady) {
                result = x;  // n=0 (includes tie 2|x|==|y|: rounds to 0)
            } else {
                double r = adx - ady;  // exact by Sterbenz, n=1
                result = (r == 0.0) ? copysign(0.0, x)
                       : (ix & UINT64_C(0x8000000000000000)) ? -r : r;
            }
        } else {
            // Raw biased exponent fields suffice for d: the 1023 bias cancels in the
            // subtraction, and for subnormal y (biased field 0) nsteps is off by at
            // most 2, which the final RneD step absorbs without loss of correctness.
            int32_t xe_d = (int32_t)(ax >> EXPSHIFTBITS_DP64);
            int32_t ye_d = (int32_t)(ay >> EXPSHIFTBITS_DP64);
            int32_t d    = xe_d - ye_d;

            if (likely(d <= EXPSHIFTBITS_DP64)) {
                double half_ady = ady * 0x1p-1;
                double n_d = RneD(adx / ady);
                double r   = fma(-n_d, ady, adx);

                if (unlikely(r > half_ady)) {
                    r -= ady;
                } else if (unlikely(r < -half_ady || (r + r == -ady && (int64_t)n_d & 1))) {
                    r += ady;
                }

                result = (r == 0.0) ? copysign(0.0, x)
                       : (ix & UINT64_C(0x8000000000000000)) ? -r : r;
            } else {
                // d > 52: multi-step 24-bit chunk reduction.
                // Each step reduces adx mod (w) where w = ady * 2^(24*nsteps).
                // We use (uint64_t)(adx/w) to get the exact integer quotient (fits in
                // 24 bits by construction).  scalbn is used instead of bit-manipulation
                // so that subnormal ady is handled correctly.
                int32_t nsteps = d / 24;
                // asdouble(ay + ...) is safe only for normal ady (subnormal ady has
                // biased exponent field 0 so adding overflows into the mantissa bits).
                double w = ay >= UINT64_C(0x0010000000000000)
                    ? asdouble(ay + ((uint64_t)(24 * nsteps) << EXPSHIFTBITS_DP64))
                    : scalbn(ady, 24 * nsteps);
                for (int32_t i = 0; i < nsteps; i++) {
                    uint64_t q = (uint64_t)(adx / w);
                    adx -= (double)q * w;
                    w *= 0x1p-24;  // 2^-24
                }
                double half_ady = ady * 0x1p-1;
                double n = RneD(adx / w);
                adx = fma(-n, w, adx);
                if (unlikely(adx > half_ady)) {
                    adx -= ady;
                } else if (unlikely(adx < -half_ady || (adx + adx == -ady && (int64_t)n & 1))) {
                    adx += ady;
                }

                result = (adx == 0.0) ? copysign(0.0, x)
                       : (ix & UINT64_C(0x8000000000000000)) ? -adx : adx;
            }
        }
    }

    return result;
}
