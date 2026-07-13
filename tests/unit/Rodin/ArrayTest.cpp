#include <random>
#include <gtest/gtest.h>

#include <Rodin/Array.h>

namespace Rodin::Tests::Unit
{
  // -----------------------------------------------------------------------------
  // IndexArrayEquality Tests (order sensitive)
  // -----------------------------------------------------------------------------
  /// @brief Verifies both empty arrays for index array equality by checking true predicates.
  TEST(IndexArrayEquality, BothEmptyArrays)
  {
    Rodin::IndexArray a(0), b(0);
    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  /// @brief Verifies single element equal for index array equality by checking true predicates.
  TEST(IndexArrayEquality, SingleElementEqual)
  {
    Rodin::IndexArray a(1), b(1);
    a << 42;
    b << 42;
    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  /// @brief Verifies single element different for index array equality by checking false predicates.
  TEST(IndexArrayEquality, SingleElementDifferent)
  {
    Rodin::IndexArray a(1), b(1);
    a << 42;
    b << 43;
    Rodin::IndexArrayEquality eq;
    EXPECT_FALSE(eq(a, b));
  }

  // ----- Tests for size 2 -----
  /// @brief Verifies size 2 equal for index array equality by checking true predicates.
  TEST(IndexArrayEquality, Size2Equal)
  {
    Rodin::IndexArray a(2), b(2);
    a << 7, 8;
    b << 7, 8;
    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  /// @brief Order matters here.
  TEST(IndexArrayEquality, Size2Different)
  {
    // Order matters here.
    Rodin::IndexArray a(2), b(2);
    a << 7, 8;
    b << 8, 7;
    Rodin::IndexArrayEquality eq;
    EXPECT_FALSE(eq(a, b));
  }

  // ----- Tests for size 3 -----
  /// @brief Verifies size 3 equal for index array equality by checking true predicates.
  TEST(IndexArrayEquality, Size3Equal)
  {
    Rodin::IndexArray a(3), b(3);
    a << 1, 2, 3;
    b << 1, 2, 3;
    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  /// @brief Even though a permutation might be considered equal in symmetric sense,.
  TEST(IndexArrayEquality, Size3Different)
  {
    // Even though a permutation might be considered equal in symmetric sense,
    // order-sensitive equality should fail.
    Rodin::IndexArray a(3), b(3);
    a << 1, 2, 3;
    b << 3, 2, 1;
    Rodin::IndexArrayEquality eq;
    EXPECT_FALSE(eq(a, b));
  }

  /// @brief Verifies different sizes for index array equality by checking false predicates.
  TEST(IndexArrayEquality, DifferentSizes)
  {
    Rodin::IndexArray a(4), b(5);
    a << 1, 2, 3, 4;
    b << 1, 2, 3, 4, 5;
    Rodin::IndexArrayEquality eq;
    EXPECT_FALSE(eq(a, b));
  }
}


