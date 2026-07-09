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

#include <gtest.h>
#include <stdint.h>
#include <limits.h>
#include <cstring>

#include "test_lround_data.h"

/*
 * Forward declarations for the AOCL lround family.
 * These are declared extern "C" to match the C linkage of the library.
 *
 * Exception flags are not tested: the existing lround implementations
 * raise no FE_INVALID for NaN / Inf / overflow inputs.
 */
extern "C" {
    long      amd_lround(double x);
    long      amd_lroundf(float x);
    long long amd_llround(double x);
    long long amd_llroundf(float x);
}

/* Reinterpret uint64_t bits as double. */
static inline double bits_to_double(uint64_t bits)
{
    double v;
    std::memcpy(&v, &bits, sizeof(v));
    return v;
}

/* Reinterpret uint32_t bits as float. */
static inline float bits_to_float(uint32_t bits)
{
    float v;
    std::memcpy(&v, &bits, sizeof(v));
    return v;
}

/* ------------------------------------------------------------------ */
/* lround(double)                                                      */
/* ------------------------------------------------------------------ */

TEST(lround, f64_special_cases)
{
    for (size_t i = 0; i < sizeof(lround_f64_cases) / sizeof(lround_f64_cases[0]); ++i) {
        const struct lround_f64_data &tc = lround_f64_cases[i];
        double x = bits_to_double(tc.in);
        long result = amd_lround(x);
        EXPECT_EQ(result, tc.out_l)
            << "lround(" << x << "): got " << result << ", expected " << tc.out_l
            << " (case " << i << ")";
    }
}

/* ------------------------------------------------------------------ */
/* llround(double)                                                     */
/* ------------------------------------------------------------------ */

TEST(llround, f64_special_cases)
{
    for (size_t i = 0; i < sizeof(lround_f64_cases) / sizeof(lround_f64_cases[0]); ++i) {
        const struct lround_f64_data &tc = lround_f64_cases[i];
        double x = bits_to_double(tc.in);
        long long result = amd_llround(x);
        EXPECT_EQ(result, tc.out_ll)
            << "llround(" << x << "): got " << result << ", expected " << tc.out_ll
            << " (case " << i << ")";
    }
}

/* ------------------------------------------------------------------ */
/* lroundf(float)                                                      */
/* ------------------------------------------------------------------ */

TEST(lroundf, f32_special_cases)
{
    for (size_t i = 0; i < sizeof(lround_f32_cases) / sizeof(lround_f32_cases[0]); ++i) {
        const struct lround_f32_data &tc = lround_f32_cases[i];
        float x = bits_to_float(tc.in);
        long result = amd_lroundf(x);
        EXPECT_EQ(result, tc.out_l)
            << "lroundf(" << x << "): got " << result << ", expected " << tc.out_l
            << " (case " << i << ")";
    }
}

/* ------------------------------------------------------------------ */
/* llroundf(float)                                                     */
/* ------------------------------------------------------------------ */

TEST(llroundf, f32_special_cases)
{
    for (size_t i = 0; i < sizeof(lround_f32_cases) / sizeof(lround_f32_cases[0]); ++i) {
        const struct lround_f32_data &tc = lround_f32_cases[i];
        float x = bits_to_float(tc.in);
        long long result = amd_llroundf(x);
        EXPECT_EQ(result, tc.out_ll)
            << "llroundf(" << x << "): got " << result << ", expected " << tc.out_ll
            << " (case " << i << ")";
    }
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
