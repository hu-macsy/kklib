#include <gtest/gtest.h>

#include "test.hpp"
#include <kklib/mpi_helper.hpp>
#include <kklib/util.hpp>

TEST(RandomEngine, integral)
{
    constexpr int min = 1;
    constexpr int max = 5;
    kklib::RandomEngine<int> engine{};
    EXPECT_NO_THROW(engine(min, max));

    unsigned int const result = engine(min, max);
    EXPECT_GE(result, min);
    EXPECT_LE(result, max);
}

TEST(RandomEngine, floating_point)
{
    constexpr auto min = 1.f;
    constexpr auto max = 5.f;
    kklib::RandomEngine<float> engine{};
    EXPECT_NO_THROW(engine(min, max));

    unsigned int const result = engine(min, max);
    EXPECT_GE(result, min);
    EXPECT_LE(result, max);
}

GTEST_API_ int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    return result;
}
