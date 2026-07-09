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
#include <fenv.h>
#include <stdint.h>
#include <limits.h>
#include <cstring>

#include "test_lrint_data.h"

/*
 * Forward declarations for the AOCL lrint family.
 * These are declared extern "C" to match the C linkage of the library.
 */
extern "C" {
    long      amd_lrint(double x);
    long      amd_lrintf(float x);
    long long amd_llrint(double x);
    long long amd_llrintf(float x);
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
/* lrint(double)                                                       */
/* ------------------------------------------------------------------ */

TEST(lrint, f64_special_cases)
{
    for (size_t i = 0; i < sizeof(lrint_f64_cases) / sizeof(lrint_f64_cases[0]); ++i) {
        const struct lrint_f64_data &tc = lrint_f64_cases[i];
        double x = bits_to_double(tc.in);

        feclearexcept(FE_ALL_EXCEPT);
        long result = amd_lrint(x);

        EXPECT_EQ(result, tc.out)
            << "lrint(" << x << "): got " << result << ", expected " << tc.out
            << " (case " << i << ")";

        if (tc.excepts != 0) {
            int raised = fetestexcept(tc.excepts);
            EXPECT_NE(raised, 0)
                << "lrint(" << x << "): expected FE_INVALID but it was not raised"
                << " (case " << i << ")";
        }
    }
}

/* ------------------------------------------------------------------ */
/* llrint(double)                                                      */
/* ------------------------------------------------------------------ */

TEST(llrint, f64_special_cases)
{
    for (size_t i = 0; i < sizeof(llrint_f64_cases) / sizeof(llrint_f64_cases[0]); ++i) {
        const struct llrint_f64_data &tc = llrint_f64_cases[i];
        double x = bits_to_double(tc.in);

        feclearexcept(FE_ALL_EXCEPT);
        long long result = amd_llrint(x);

        EXPECT_EQ(result, tc.out)
            << "llrint(" << x << "): got " << result << ", expected " << tc.out
            << " (case " << i << ")";

        if (tc.excepts != 0) {
            int raised = fetestexcept(tc.excepts);
            EXPECT_NE(raised, 0)
                << "llrint(" << x << "): expected FE_INVALID but it was not raised"
                << " (case " << i << ")";
        }
    }
}

/* ------------------------------------------------------------------ */
/* lrintf(float)                                                       */
/* ------------------------------------------------------------------ */

TEST(lrintf, f32_special_cases)
{
    for (size_t i = 0; i < sizeof(lrint_f32_cases) / sizeof(lrint_f32_cases[0]); ++i) {
        const struct lrint_f32_data &tc = lrint_f32_cases[i];
        float x = bits_to_float(tc.in);

        feclearexcept(FE_ALL_EXCEPT);
        long result = amd_lrintf(x);

        EXPECT_EQ(result, tc.out)
            << "lrintf(" << x << "): got " << result << ", expected " << tc.out
            << " (case " << i << ")";

        if (tc.excepts != 0) {
            int raised = fetestexcept(tc.excepts);
            EXPECT_NE(raised, 0)
                << "lrintf(" << x << "): expected FE_INVALID but it was not raised"
                << " (case " << i << ")";
        }
    }
}

/* ------------------------------------------------------------------ */
/* llrintf(float)                                                      */
/* ------------------------------------------------------------------ */

TEST(llrintf, f32_special_cases)
{
    for (size_t i = 0; i < sizeof(llrint_f32_cases) / sizeof(llrint_f32_cases[0]); ++i) {
        const struct llrint_f32_data &tc = llrint_f32_cases[i];
        float x = bits_to_float(tc.in);

        feclearexcept(FE_ALL_EXCEPT);
        long long result = amd_llrintf(x);

        EXPECT_EQ(result, tc.out)
            << "llrintf(" << x << "): got " << result << ", expected " << tc.out
            << " (case " << i << ")";

        if (tc.excepts != 0) {
            int raised = fetestexcept(tc.excepts);
            EXPECT_NE(raised, 0)
                << "llrintf(" << x << "): expected FE_INVALID but it was not raised"
                << " (case " << i << ")";
        }
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
