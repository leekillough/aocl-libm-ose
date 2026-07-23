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

/******************************************
 * Implementation Notes:
 *
 * Prototype:
 * double remainder(double x, double y)
 *
 * Algorithm:
 * remainder(x, y) = x - n*y, where n = roundTiesToEven(x/y).
 *
 * This is identical to the fmod integer algorithm except that fmod uses
 * n = trunc(x/y) while remainder uses n = roundTiesToEven(x/y).
 *
 * Main path (|x| >= |y|):
 *   Compute rem = Mx * 2^shift mod My via Rem128P (same as fmod), but
 *   Rem128P also returns the integer quotient.  The XOR of LSBs of all
 *   quotients across all reduction steps gives the parity of n.
 *   After the reduction:
 *     2*rem < My: keep rem; sign(result) = sign(x).
 *     2*rem > My: use rem = My-rem (round up); sign(result) = -sign(x).
 *     2*rem == My: round up iff n is odd (tie-to-even); tracked via n_quot.
 *   MAXSHIFT is 63 for divq/_udiv128 and 75 for __uint128_t (same as fmod).
 *
 * Small-x path (0 < |x| < |y|):
 *   n is 0 or 1.  The comparison 2|x| vs |y| is done in float arithmetic,
 *   which is exact here: doubling is exact for finite non-max doubles, and
 *   the Sterbenz subtraction |y|-|x| is exact when |y|/2 < |x| < |y|.
 *   (The rare case that doubling overflows for very large |x| near max
 *   finite still gives a correct comparison result, though FE_OVERFLOW
 *   may be raised; this affects only inputs where |x| has maximum exponent
 *   and |y| is slightly larger than |x| -- an extremely rare edge case.)
 *
 */

#include <stdint.h>
#include <math.h>

#include "libm_macros.h"
#include "libm_util_amd.h"
#include <libm/alm_special.h>
#include <libm/typehelper.h>
#include <libm/amd_funcs_internal.h>
#include <libm/compiler.h>

#if defined(__GNUC__) || defined(__clang__)

#define CLZ64(x) __builtin_clzll(x)

#elif defined(_MSC_VER)

#include <intrin.h>
static inline int alm_clz64(uint64_t x)
{
    unsigned long idx;
    _BitScanReverse64(&idx, x);
    return 63 - (int)idx;
}
#define CLZ64(x) alm_clz64(x)

#else

#error "Intrinsic for counting leading zeroes not found"

#endif

typedef struct
{
    uint64_t m;  // 53-bit significand (including implicit 1)
    int e;       // biased exponent
} F64ExpMan;

static inline F64ExpMan F64Extract(uint64_t fax)
{
    int lz;
    return unlikely(fax < POS_LNORMAL_F64) ?
        lz = CLZ64(fax),
        (F64ExpMan) {
            .m = fax << (lz - (64 - MANTLENGTH_DP64)),
            .e = (64 - MANTLENGTH_DP64) + 1 - lz
        } :
        (F64ExpMan) {
            .m = (fax & MANTBITS_DP64) | IMPBIT_DP64,
            .e = (int)(fax >> EXPSHIFTBITS_DP64)
        };
}

/* Rem128P: same as Rem128 but also stores the quotient in *quot.
   The XOR of (quot & 1) across all steps gives the parity of n. */
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))

/* hi = Mx >> (64-d) < 2^52 <= My for d in [0,63] when My has implicit bit */
#define MAXSHIFT 63

static inline uint64_t Rem128P(uint64_t Mx, uint64_t My, int d, uint64_t *quot)
{
    uint64_t rem;
    uint64_t hi = (d != 0) ? Mx >> (64 - d) : 0;
    uint64_t lo = Mx << d;
    __asm__("divq %[divisor]"
            : "=a"(*quot), "=d"(rem)
            : "0"(lo), "1"(hi), [divisor] "r"(My));
    return rem;
}

#elif defined(__SIZEOF_INT128__)

/* (Mx << d) < 2^128 for d <= 75 since Mx < 2^53; remainder fits in uint64_t */
#define MAXSHIFT (128 - MANTLENGTH_DP64)

static inline uint64_t Rem128P(uint64_t Mx, uint64_t My, int d, uint64_t *quot)
{
    __uint128_t dividend = (__uint128_t)Mx << d;
    *quot = (uint64_t)(dividend / My);
    return (uint64_t)(dividend % My);
}

#elif defined(_MSC_VER) && defined(_M_X64)

/* hi = Mx >> (64-d) < 2^52 <= My for d in [0,63] when My has implicit bit.
 * Note: _udiv128 is unresolved at link time with clang-cl + lld-link due to
 * a known open bug (https://github.com/llvm/llvm-project/issues/59168).
 * This branch is only reachable with MSVC link.exe. */
#define MAXSHIFT 63

#include <intrin.h>
static inline uint64_t Rem128P(uint64_t Mx, uint64_t My, int d, uint64_t *quot)
{
    uint64_t rem;
    uint64_t hi = (d != 0) ? Mx >> (64 - d) : 0;
    uint64_t lo = Mx << d;
    *quot = _udiv128(hi, lo, My, &rem);
    return rem;
}

#else

#error "128-bit integer division not available"

#endif

double ALM_PROTO_OPT(remainder)(double x, double y)
{
    uint64_t fax = asuint64(x) & POS_BITSET_DP64;
    uint64_t fay = asuint64(y) & POS_BITSET_DP64;
    double result = x;
    uint64_t sign = asuint64(x) & SIGNBIT_DP64;

    if (unlikely(((fay - 1) | fax) >= POS_INF_F64))
    {
        if (fay > POS_INF_F64)
        {   // |y| NaN
            result = x * y;
        }
        else if (fax > POS_INF_F64)
        {   // |x| NaN
            result = x + x;
        }
        else if ((fax == POS_INF_F64) || (fay == 0))
        {   // |x| == Inf || y == 0
            result = __alm_handle_error(INDEFBITPATT_DP64, AMD_F_INVALID);
        }
        else if (fax >= fay)
        {   // |x| >= |y|
            goto normal;
        }
    }
    else if (likely(fax >= fay))
    {   // |x| >= |y|
    normal: ;
        int xe = (int)(fax >> EXPSHIFTBITS_DP64);
        int ye = (int)(fay >> EXPSHIFTBITS_DP64);
        int shift = xe - ye;
        uint64_t rem;
        F64ExpMan fpy;
        uint64_t n_quot = 0;

        if (likely((ye != 0) && (xe != 0) && (shift <= MAXSHIFT)))
        {
            // Fast path: both normal, small shift
            rem = (fax & MANTBITS_DP64) | IMPBIT_DP64;
            fpy = (F64ExpMan) { .m = (fay & MANTBITS_DP64) | IMPBIT_DP64, .e = ye };
        }
        else
        {
            // Slow path: subnormals or large shift
            F64ExpMan fpx = F64Extract(fax);
            fpy = F64Extract(fay);
            rem = fpx.m;
            shift = fpx.e - fpy.e;
            while (shift > MAXSHIFT) {
                uint64_t q;
                rem = Rem128P(rem, fpy.m, MAXSHIFT, &q);
                n_quot ^= q;
                shift -= MAXSHIFT;
            }
        }
        {
            uint64_t q;
            rem = Rem128P(rem, fpy.m, shift, &q);
            n_quot ^= q;
        }

        // Nearest-even rounding: round up if 2*rem > My, or on a tie (2*rem == My)
        // when n is odd.  XOR of all quotient LSBs gives parity of the total n.
        uint64_t two_rem = rem + rem;
        if (two_rem > fpy.m || (two_rem == fpy.m && (n_quot & 1)))
        {
            rem = fpy.m - rem;
            sign ^= SIGNBIT_DP64;
        }

        if (likely(rem != 0)) {
            int k = CLZ64(rem) - (64 - MANTLENGTH_DP64);
            rem = (fpy.e > k) ? ((uint64_t)(fpy.e - k) << EXPSHIFTBITS_DP64)
                | ((rem << k) & MANTBITS_DP64) :
                likely(fpy.e > 0) ? rem << (fpy.e - 1) : rem >> (1 - fpy.e);
        }
        result = asdouble(rem | sign);
    }
    else
    {   // 0 < |x| < |y|: n = 0 or 1; float arithmetic is exact here (Sterbenz)
        double adx = asdouble(fax);
        double ady = asdouble(fay);
        double ax2 = adx + adx;       // 2*|x|, exact for finite non-max doubles
        if (ax2 > ady)
        {
            // n=1: |result| = |y|-|x|, sign = -sign(x)
            double r = ady - adx;     // exact by Sterbenz (adx > ady/2)
            result = asdouble(asuint64(r) | (sign ^ SIGNBIT_DP64));
        }
        // else 2|x| <= |y|: n=0 (includes tie 2|x|==|y|: rounds to 0), result=x
    }

    return result;
}
