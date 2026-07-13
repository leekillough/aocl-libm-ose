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
//   Single combined entry test:
//     (ax - IMPBIT_DP64) > 0x7fdfffffffffffff ||
//     (ay - IMPBIT_DP64) > 0x7fdfffffffffffff
//   fires iff any of: x or y is zero, subnormal, Inf, or NaN.
//   Normal finite nonzero ax,ay lie in [IMPBIT_DP64, 0x7fefffffffffffff];
//   subtracting IMPBIT_DP64 maps them to [0, 0x7fdfffffffffffff].
//   Special inputs fall outside: subnormals underflow to 0 or wrap, Inf/NaN exceed
//   the upper bound.  All such inputs are dispatched to RemainderGeneral().
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
//     y is guaranteed normal here (entry test excludes subnormals), so w is
//     computed by adding to the biased exponent field of ay.
//
//   RemainderGeneral(): cold NOINLINE helper for special and subnormal inputs.
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
//
// In the chunk-reduction final step, q = adx/w < 2^24 (exact, representable),
// so q + 0.5 is always exact and this is rounding-mode independent.
//
// In the single-step path, q = adx/ady and d may be as large as 52, allowing
// q to reach just under 2^53.  For q in [2^52, 2^53), doubles have ULP = 1,
// so q + 0.5 is rounded and may be rounding-mode dependent; the boundary
// corrections at each call site absorb any resulting off-by-one.
//
// Differs from nearest-even only at half-integer ties; corrected at each call
// site with:
//   if (r + r == -ady && (int64_t)n & 1) r += ady;
// (no-op in the SSE4.1 path where nearest-even always gives an even n).
static inline double RneD(double q) { return trunc(q + 0.5); }

#endif // __SSE4_1__

// Cold helper for special and subnormal inputs.  Handles NaN, Inf, zero, and
// subnormal x or y (i.e., every case excluded from the hot path).
static double NOINLINE_COLD RemainderGeneral(double x, double y)
{
    uint64_t ix = asuint64(x);
    uint64_t iy = asuint64(y);
    uint64_t ax = ix & POS_BITSET_DP64;
    uint64_t ay = iy & POS_BITSET_DP64;
    double result = x;

    // Signaling NaN: exponent=0x7ff, quiet bit clear, mantissa nonzero.
    // Uint64 range: [PINFBITPATT_DP64+1, PINFBITPATT_DP64+QNAN_MASK_64-1].
    // Raise FE_INVALID and quiet the payload.
    if ((ax - 0x7ff0000000000001) < 0x0007ffffffffffff)
    {
        result = __alm_handle_error(ix | QNAN_MASK_64, AMD_F_INVALID);
    }
    else if ((ay - 0x7ff0000000000001) < 0x0007ffffffffffff)
    {
        result = __alm_handle_error(iy | QNAN_MASK_64, AMD_F_INVALID);
    }
    // Quiet NaN x: result = x already set.
    else if (ax > PINFBITPATT_DP64)
    {
    }
    else if (ay > PINFBITPATT_DP64)
    {
        result = y;
    }
    // x=Inf is invalid (y is confirmed not NaN here).
    else if (ax == PINFBITPATT_DP64)
    {
        result = __alm_handle_error(QNANBITPATT_DP64, AMD_F_INVALID);
    }
    // y=0: invalid (x is confirmed finite here).
    else if (ay == 0)
    {
        result = __alm_handle_error(QNANBITPATT_DP64, AMD_F_INVALID);
    }
    // y=Inf or x=0: result = x already set.
    else if (ay == PINFBITPATT_DP64 || ax == 0)
    {
    }
    else
    {
        // Subnormal x or y (both finite nonzero).
        double adx = asdouble(ax);
        double ady = asdouble(ay);

        if (ax <= ay)
        {
            // Fast path: n is 0 or 1.
            double ax2 = adx + adx;
            if (ax2 > ady)
            {
                double r = adx - ady;  // exact by Sterbenz, n=1
                result = (r == 0.0) ? copysign(0.0, x)
                       : (ix & SIGNBIT_DP64) ? -r : r;
            }
            // else n=0: result = x already set.
        }
        else
        {
            // Compute true (unbiased) exponents for subnormal-safe chunk reduction.
            // For normal values, the biased exponent field minus 1023 gives the true
            // exponent.  For subnormals, __builtin_clzll gives the bit position of the
            // leading 1, from which the true exponent follows.
#if defined(__GNUC__) || defined(__clang__)
            int32_t xe = (ax >= IMPBIT_DP64)
                       ? (int32_t)(ax >> EXPSHIFTBITS_DP64) - 1023
                       : (int32_t)(63 - __builtin_clzll(ax)) - 1074;
            int32_t ye = (ay >= IMPBIT_DP64)
                       ? (int32_t)(ay >> EXPSHIFTBITS_DP64) - 1023
                       : (int32_t)(63 - __builtin_clzll(ay)) - 1074;
#else
            int32_t xe = ilogb(adx);
            int32_t ye = ilogb(ady);
#endif
            int32_t d = xe - ye;
            double half_ady = ady * 0x1p-1;

            if (d <= EXPSHIFTBITS_DP64)
            {
                double n_d = RneD(adx / ady);
                double r   = fma(-n_d, ady, adx);

                if (unlikely(r > half_ady))
                    r -= ady;
                else if (unlikely((r < -half_ady) || ((r + r == -ady) && (((int64_t)n_d & 1) != 0))))
                    r += ady;

                result = (r == 0.0) ? copysign(0.0, x)
                       : (ix & SIGNBIT_DP64) ? -r : r;
            }
            else
            {
                // d > 52: multi-step 24-bit chunk reduction using true exponents.
                // scalbn handles subnormal ady correctly.
                int32_t nsteps = d / 24;
                double w = scalbn(ady, 24 * nsteps);

                for (int32_t i = 0; i < nsteps; i++)
                {
                    uint64_t q = (uint64_t)(adx / w);
                    adx -= (double)q * w;
                    w *= 0x1p-24;
                }

                double n = RneD(adx / w);
                adx = fma(-n, w, adx);

                if (unlikely(adx > half_ady))
                    adx -= ady;
                else if (unlikely((adx < -half_ady) || ((adx + adx == -ady) && (((int64_t)n & 1) != 0))))
                    adx += ady;

                result = (adx == 0.0) ? copysign(0.0, x)
                       : (ix & SIGNBIT_DP64) ? -adx : adx;
            }
        }
    }

    return result;
}

double ALM_PROTO_OPT(remainder)(double x, double y)
{
    uint64_t ix = asuint64(x);
    uint64_t iy = asuint64(y);
    uint64_t ax = ix & POS_BITSET_DP64;
    uint64_t ay = iy & POS_BITSET_DP64;
    double result = x;

    // Dispatch zero, subnormal, Inf, and NaN inputs to the cold helper.
    // Normal finite nonzero ax,ay lie in [IMPBIT_DP64, 0x7fefffffffffffff];
    // subtracting IMPBIT_DP64 maps them to [0, 0x7fdfffffffffffff].
    // Any value outside this range (zero/subnormal wrap below, Inf/NaN exceed above)
    // triggers the branch.
    if (unlikely((ax - IMPBIT_DP64) > 0x7fdfffffffffffff ||
                 (ay - IMPBIT_DP64) > 0x7fdfffffffffffff))
    {
        result = RemainderGeneral(x, y);
    }
    else
    {
        // Hot path: both x and y are normal finite nonzero doubles.
        double adx = asdouble(ax);
        double ady = asdouble(ay);

        // Fast path: |x| <= |y|.  n is 0 or 1; no division needed.
        // 2*adx is exact (no overflow: ax < PINFBITPATT_DP64).
        // adx - ady is exact by Sterbenz (ady/2 < adx <= ady for n=1 case).
        if (likely(ax <= ay))
        {
            double ax2 = adx + adx;
            if (ax2 > ady)
            {
                double r = adx - ady;  // exact by Sterbenz, n=1
                result = (r == 0.0) ? copysign(0.0, x)
                       : (ix & SIGNBIT_DP64) ? -r : r;
            }
            // else n=0 (includes tie 2|x|==|y|: rounds to 0): result = x already set.
        }
        else
        {
            // Raw biased exponent fields (bias cancels in the subtraction).
            int32_t xe_d = (int32_t)(ax >> EXPSHIFTBITS_DP64);
            int32_t ye_d = (int32_t)(ay >> EXPSHIFTBITS_DP64);
            int32_t d    = xe_d - ye_d;
            double half_ady = ady * 0x1p-1;

            if (likely(d <= EXPSHIFTBITS_DP64))
            {
                double n_d = RneD(adx / ady);
                double r   = fma(-n_d, ady, adx);

                if (unlikely(r > half_ady))
                    r -= ady;
                else if (unlikely((r < -half_ady) || ((r + r == -ady) && (((int64_t)n_d & 1) != 0))))
                    r += ady;

                result = (r == 0.0) ? copysign(0.0, x)
                       : (ix & SIGNBIT_DP64) ? -r : r;
            }
            else
            {
                // d > 52: multi-step 24-bit chunk reduction.
                // Each step reduces adx mod w where w = ady * 2^(24*(nsteps-i)).
                // y is normal here (entry test excludes subnormals), so w is formed by
                // adding directly to the biased exponent field of ay.  The resulting
                // biased exponent ye_d + 24*nsteps = xe_d - r (r = d mod 24 in [0,23])
                // is at most xe_d <= 2046, so no overflow into the Inf exponent.
                // Each quotient adx/w < 2^24 (adx < 2^(xe_d-1022), w >= ady*2^(xe_d-r-1),
                // so adx/w < 2^(r+1) <= 2^24), and fits exactly in uint64_t.
                int32_t nsteps = d / 24;
                double w = asdouble(ay + ((uint64_t)(24 * nsteps) << EXPSHIFTBITS_DP64));

                for (int32_t i = 0; i < nsteps; i++)
                {
                    uint64_t q = (uint64_t)(adx / w);
                    adx -= (double)q * w;
                    w *= 0x1p-24;
                }

                double n = RneD(adx / w);
                adx = fma(-n, w, adx);

                if (unlikely(adx > half_ady))
                    adx -= ady;
                else if (unlikely((adx < -half_ady) || ((adx + adx == -ady) && (((int64_t)n & 1) != 0))))
                    adx += ady;

                result = (adx == 0.0) ? copysign(0.0, x)
                       : (ix & SIGNBIT_DP64) ? -adx : adx;
            }
        }
    }

    return result;
}
