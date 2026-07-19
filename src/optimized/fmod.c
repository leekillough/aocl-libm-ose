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

/******************************************
 * Implementation Notes:
 *
 * Prototype:
 * double fmod(double x, double y)
 *
 * Algorithm:
 * fmod(x, y) = x - n*y, where n = trunc(x/y).
 *
 * Integer fast path (both normal, exponent difference shift <= MAXSHIFT):
 *   Extract 53-bit double significands Mx, My directly from bit patterns.
 *   Compute rem = Mx * 2^shift mod My via Rem128, a 128-bit-by-64-bit integer
 *   division.  Pack rem back into a double.  No floating-point operations are
 *   performed, so no exceptions are raised.  MAXSHIFT is 63 for the divq and
 *   _udiv128 implementations (limited by divq overflow: hi < My requires
 *   shift <= 63 since My has its implicit bit set), and 75 for __uint128_t.
 *
 * Slow path (subnormals, or exponent difference > MAXSHIFT):
 *   Extract significands and exponents via F64Extract (handles subnormals with
 *   CLZ64).  Pure integer iterative reduction: rem = Rem128(rem, My, MAXSHIFT),
 *   repeated until shift <= MAXSHIFT, then one final Rem128(rem, My, shift).
 *   Invariant: rem < My < 2^53, so Rem128 never overflows.  No floating-point
 *   operations are performed, so no exceptions are raised.
 *
 */

#include <limits.h>
#include <stdint.h>

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
            .m = fax << (lz - (int)(sizeof(uint64_t) * CHAR_BIT - MANTLENGTH_DP64)),
            .e = (int)(sizeof(uint64_t) * CHAR_BIT - MANTLENGTH_DP64 + 1) - lz
        } :
        (F64ExpMan) {
            .m = (fax & MANTBITS_DP64) | IMPBIT_DP64,
            .e = (int)(fax >> EXPSHIFTBITS_DP64)
        };
}

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))

/* hi = Mx >> (64-d) < 2^52 <= My for d in [0,63] when My has implicit bit */
#define MAXSHIFT 63

static inline uint64_t Rem128(uint64_t Mx, uint64_t My, int d)
{
#if 0
    uint64_t lo = Mx << d;
    if (likely((d = -d) != 0)) {
        uint64_t quot;
        uint64_t rem;
        uint64_t hi = Mx >> (d + 64);   /* d negated above: d+64 == 64-original_d */
        __asm__("divq %[divisor]"
                : "=a"(quot), "=d"(rem)
                : "0"(lo), "1"(hi), [divisor] "r"(My));
        return rem;
    } else {
        return lo % My;
    }
#else
    uint64_t quot;
    uint64_t rem;
    uint64_t hi = (d != 0) ? Mx >> (64 - d) : 0;
    uint64_t lo = Mx << d;
    __asm__("divq %[divisor]"
            : "=a"(quot), "=d"(rem)
            : "0"(lo), "1"(hi), [divisor] "r"(My));
    return rem;
#endif
}

#elif defined(__SIZEOF_INT128__)

/* (Mx << d) < 2^128 for d <= 75 since Mx < 2^53; remainder fits in uint64_t */
#define MAXSHIFT (128 - MANTLENGTH_DP64)

static inline uint64_t Rem128(uint64_t Mx, uint64_t My, int d)
{
    return (uint64_t)(((__uint128_t)Mx << d) % My);
}

#elif defined(_MSC_VER) && defined(_M_X64)

/* hi = Mx >> (64-d) < 2^52 <= My for d in [0,63] when My has implicit bit.
 * Note: _udiv128 is unresolved at link time with clang-cl + lld-link due to
 * a known open bug (https://github.com/llvm/llvm-project/issues/59168).
 * This branch is only reachable with MSVC link.exe. */
#define MAXSHIFT 63

#include <intrin.h>
static inline uint64_t Rem128(uint64_t Mx, uint64_t My, int d)
{
    uint64_t rem;
    uint64_t hi = (d != 0) ? Mx >> (64 - d) : 0;
    uint64_t lo = Mx << d;
    _udiv128(hi, lo, My, &rem);
    return rem;
}

#else

#error "128-bit integer division not available"

#endif

double ALM_PROTO_OPT(fmod)(double x, double y)
{
    uint64_t  fax = asuint64(x) & POS_BITSET_DP64;
    uint64_t  fay = asuint64(y) & POS_BITSET_DP64;
    double result = x;

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
    else if (fax >= fay)
    {   // |x| >= |y|
    normal:
        int xe = (int)(fax >> EXPSHIFTBITS_DP64);
        int ye = (int)(fay >> EXPSHIFTBITS_DP64);
        int shift = xe - ye;
        const int maxshift = MAXSHIFT;
        uint64_t rem;
        F64ExpMan fpy;

        if (likely((xe != 0) && (ye != 0) && (shift <= maxshift)))
        {
            // Fast path: both normal, small shift
            rem = (fax & MANTBITS_DP64) | IMPBIT_DP64;

            //  My  = (fay & MANTBITS_DP64) | IMPBIT_DP64;
            fpy = (F64ExpMan) { .m = (fay & MANTBITS_DP64) | IMPBIT_DP64, .e = ye };
        }
        else
        {
            // Slow path: subnormals or large shift
            F64ExpMan fpx = F64Extract(fax);
            fpy = F64Extract(fay);
            rem = fpx.m;
            shift = fpx.e - fpy.e;
            while (shift > maxshift) {
                rem = Rem128(rem, fpy.m, maxshift);
                shift -= maxshift;
            }
        }
        rem = Rem128(rem, fpy.m, shift);
        if (likely(rem != 0)) {
            int k = CLZ64(rem) + MANTLENGTH_DP64 - (int)(sizeof(uint64_t) * CHAR_BIT);
            rem = (fpy.e > k) ? ((uint64_t)(fpy.e - k) << EXPSHIFTBITS_DP64)
                | ((rem << k) & MANTBITS_DP64) :
                likely(fpy.e > 0) ? rem << (fpy.e - 1) : rem >> (1 - fpy.e);
        }
        result = asdouble(rem | (asuint64(result) & SIGNBIT_DP64));
    }
    return result;
}
