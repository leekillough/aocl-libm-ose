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

#pragma once

#include <stdint.h>
#include <limits.h>

/*
 * Test data structs for lround/llround family.
 *
 * lround/llround round half-way cases away from zero (not round-to-even).
 * Exception flags are not tested because the existing lround.c reference
 * implementation raises AMD_F_NONE (no FE_INVALID) for out-of-range inputs.
 *
 * For f64 (double input) variants:
 *   in:      double input as uint64_t bit pattern
 *   out_l:   expected long result (LONG_MIN for out-of-range / NaN / Inf)
 *   out_ll:  expected long long result (LLONG_MIN for out-of-range / NaN / Inf)
 *
 * For f32 (float input) variants:
 *   in:      float input as uint32_t bit pattern
 *   out_l:   expected long result
 *   out_ll:  expected long long result
 */

struct lround_f64_data {
    uint64_t  in;
    long      out_l;
    long long out_ll;
};

struct lround_f32_data {
    uint32_t  in;
    long      out_l;
    long long out_ll;
};

/* Reuse the same bit-pattern constants as the lrint data header. */
/* +0.0 */ static const uint64_t R_D_POS_ZERO  = 0x0000000000000000ULL;
/* -0.0 */ static const uint64_t R_D_NEG_ZERO  = 0x8000000000000000ULL;
/* +0.5 */ static const uint64_t R_D_POS_HALF  = 0x3FE0000000000000ULL;
/* -0.5 */ static const uint64_t R_D_NEG_HALF  = 0xBFE0000000000000ULL;
/* +1.0 */ static const uint64_t R_D_ONE       = 0x3FF0000000000000ULL;
/* -1.0 */ static const uint64_t R_D_NEG_ONE   = 0xBFF0000000000000ULL;
/* +1.5 */ static const uint64_t R_D_POS_1P5   = 0x3FF8000000000000ULL;
/* -1.5 */ static const uint64_t R_D_NEG_1P5   = 0xBFF8000000000000ULL;
/* +2.0 */ static const uint64_t R_D_TWO       = 0x4000000000000000ULL;
/* +2.5 */ static const uint64_t R_D_POS_2P5   = 0x4004000000000000ULL;
/* -2.5 */ static const uint64_t R_D_NEG_2P5   = 0xC004000000000000ULL;
/* +3.0 */ static const uint64_t R_D_THREE     = 0x4008000000000000ULL;
/* +4.5 */ static const uint64_t R_D_POS_4P5   = 0x4012000000000000ULL;
/* -4.5 */ static const uint64_t R_D_NEG_4P5   = 0xC012000000000000ULL;
/* 2^52  */ static const uint64_t R_D_2P52     = 0x4330000000000000ULL;
/* -2^63 */ static const uint64_t R_D_NEG_2P63 = 0xC3E0000000000000ULL;
/* +Inf */ static const uint64_t R_D_POS_INF   = 0x7FF0000000000000ULL;
/* -Inf */ static const uint64_t R_D_NEG_INF   = 0xFFF0000000000000ULL;
/* +qNaN */ static const uint64_t R_D_QNAN     = 0x7FF8000000000000ULL;
/* sNaN */ static const uint64_t R_D_SNAN      = 0x7FF0000000000001ULL;

/* Float bit patterns */
/* +0.0f */ static const uint32_t R_F_POS_ZERO  = 0x00000000U;
/* -0.0f */ static const uint32_t R_F_NEG_ZERO  = 0x80000000U;
/* +0.5f */ static const uint32_t R_F_POS_HALF  = 0x3F000000U;
/* -0.5f */ static const uint32_t R_F_NEG_HALF  = 0xBF000000U;
/* +1.0f */ static const uint32_t R_F_ONE       = 0x3F800000U;
/* -1.0f */ static const uint32_t R_F_NEG_ONE   = 0xBF800000U;
/* +1.5f */ static const uint32_t R_F_POS_1P5   = 0x3FC00000U;
/* -1.5f */ static const uint32_t R_F_NEG_1P5   = 0xBFC00000U;
/* +2.5f */ static const uint32_t R_F_POS_2P5   = 0x40200000U;
/* -2.5f */ static const uint32_t R_F_NEG_2P5   = 0xC0200000U;
/* +4.5f */ static const uint32_t R_F_POS_4P5   = 0x40900000U;
/* -4.5f */ static const uint32_t R_F_NEG_4P5   = 0xC0900000U;
/* 2^23  */ static const uint32_t R_F_2P23      = 0x4B000000U;
/* -2^63 */ static const uint32_t R_F_NEG_2P63  = 0xDF000000U;
/* +Inf */ static const uint32_t R_F_POS_INF    = 0x7F800000U;
/* -Inf */ static const uint32_t R_F_NEG_INF    = 0xFF800000U;
/* +qNaN */ static const uint32_t R_F_QNAN      = 0x7FC00000U;
/* sNaN */ static const uint32_t R_F_SNAN       = 0x7F800001U;

/*
 * Test data for lround(double) and llround(double).
 *
 * lround/llround: half-integer rounds AWAY from zero.
 * +0.5 -> 1, -0.5 -> -1, +1.5 -> 2, -1.5 -> -2, +2.5 -> 3, -2.5 -> -3.
 */
static const struct lround_f64_data lround_f64_cases[] = {
    /* input               out_l        out_ll */
    /* zeros */
    { R_D_POS_ZERO,        0L,          0LL },
    { R_D_NEG_ZERO,        0L,          0LL },
    /* small exact values */
    { R_D_ONE,             1L,          1LL },
    { R_D_NEG_ONE,        -1L,         -1LL },
    { R_D_TWO,             2L,          2LL },
    { R_D_THREE,           3L,          3LL },
    /* half-integer rounds away from zero */
    { R_D_POS_HALF,        1L,          1LL },
    { R_D_NEG_HALF,       -1L,         -1LL },
    { R_D_POS_1P5,         2L,          2LL },
    { R_D_NEG_1P5,        -2L,         -2LL },
    { R_D_POS_2P5,         3L,          3LL },
    { R_D_NEG_2P5,        -3L,         -3LL },
    { R_D_POS_4P5,         5L,          5LL },
    { R_D_NEG_4P5,        -5L,         -5LL },
    /* large exact integral double */
    { R_D_2P52,            (long)4503599627370496LL, 4503599627370496LL },
    /* NaN -> implementation-defined (LONG_MIN / LLONG_MIN in this impl) */
    { R_D_QNAN,            LONG_MIN,    LLONG_MIN },
    { R_D_SNAN,            LONG_MIN,    LLONG_MIN },
    /* +/-Inf -> LONG_MIN / LLONG_MIN */
    { R_D_POS_INF,         LONG_MIN,    LLONG_MIN },
    { R_D_NEG_INF,         LONG_MIN,    LLONG_MIN },
    /* -2^63 is LLONG_MIN, valid for llround */
    { R_D_NEG_2P63,        LONG_MIN,    LLONG_MIN },
};

/*
 * Test data for lroundf(float) and llroundf(float).
 */
static const struct lround_f32_data lround_f32_cases[] = {
    /* input               out_l        out_ll */
    /* zeros */
    { R_F_POS_ZERO,        0L,          0LL },
    { R_F_NEG_ZERO,        0L,          0LL },
    /* small exact values */
    { R_F_ONE,             1L,          1LL },
    { R_F_NEG_ONE,        -1L,         -1LL },
    /* half-integer rounds away from zero */
    { R_F_POS_HALF,        1L,          1LL },
    { R_F_NEG_HALF,       -1L,         -1LL },
    { R_F_POS_1P5,         2L,          2LL },
    { R_F_NEG_1P5,        -2L,         -2LL },
    { R_F_POS_2P5,         3L,          3LL },
    { R_F_NEG_2P5,        -3L,         -3LL },
    { R_F_POS_4P5,         5L,          5LL },
    { R_F_NEG_4P5,        -5L,         -5LL },
    /* 2^23 is exactly representable as float and integral */
    { R_F_2P23,            (long)8388608L, 8388608LL },
    /* NaN -> LONG_MIN / LLONG_MIN */
    { R_F_QNAN,            LONG_MIN,    LLONG_MIN },
    { R_F_SNAN,            LONG_MIN,    LLONG_MIN },
    /* +/-Inf -> LONG_MIN / LLONG_MIN */
    { R_F_POS_INF,         LONG_MIN,    LLONG_MIN },
    { R_F_NEG_INF,         LONG_MIN,    LLONG_MIN },
    /* -2^63 is LLONG_MIN, valid for llroundf */
    { R_F_NEG_2P63,        LONG_MIN,    LLONG_MIN },
};
