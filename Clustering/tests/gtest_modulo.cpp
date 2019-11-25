#include "Helper_Functions.h"
#include <gtest/gtest.h>

namespace {

using ::testing::TestWithParam;
using ::testing::Values;

TEST(ModuloTest, PositiveNos) {
    ASSERT_EQ(5, modulo(59, 6));
    ASSERT_EQ(0, modulo(78, 2));
    ASSERT_EQ(72671, modulo(72671, 923555));
    ASSERT_EQ(0, modulo(0, 3555));
}

TEST(ModuloTest, NegativeNos) {
    ASSERT_EQ(1, modulo(77, -2));
    ASSERT_EQ(15, modulo(9721,-46));
    ASSERT_EQ(134, modulo(-587,-721));
}

}

