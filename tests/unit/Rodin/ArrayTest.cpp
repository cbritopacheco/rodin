#include <gtest/gtest.h>

#include <Rodin/Array.h>

namespace Rodin::Tests::Unit
{
  // -----------------------------------------------------------------------------
  // IndexArrayEquality Tests
  // -----------------------------------------------------------------------------
  TEST(IndexArrayEqualityTest, EqualArrays)
  {
    // Create two arrays with the same size and same elements.
    Rodin::IndexArray a(3), b(3);
    a << 1, 2, 3;
    b << 1, 2, 3;

    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  TEST(IndexArrayEqualityTest, DifferentArrays)
  {
    // Create two arrays with same size but differing in at least one element.
    Rodin::IndexArray a(3), b(3);
    a << 1, 2, 3;
    b << 1, 2, 4;

    Rodin::IndexArrayEquality eq;
    EXPECT_FALSE(eq(a, b));
  }

  TEST(IndexArrayEqualityTest, DifferentSizes)
  {
    // Create arrays of different sizes.
    Rodin::IndexArray a(3), b(4);
    a << 1, 2, 3;
    b << 1, 2, 3, 4;

    Rodin::IndexArrayEquality eq;
    EXPECT_FALSE(eq(a, b));
  }

  // Additional tests for arrays with sizes 6 through 10.
  TEST(IndexArrayEqualityLargeTest, EqualArraysLarge)
  {
    for (int n = 6; n <= 10; n++)
    {
      Rodin::IndexArray a(n), b(n);
      for (int i = 0; i < n; i++)
      {
        a(i) = i + 1;
        b(i) = i + 1;
      }
      Rodin::IndexArrayEquality eq;
      EXPECT_TRUE(eq(a, b)) << "Equality failed for array size " << n;
    }
  }

  TEST(IndexArrayEqualityLargeTest, DifferentArraysLarge)
  {
    for (int n = 6; n <= 10; n++)
    {
      Rodin::IndexArray a(n), b(n);
      for (int i = 0; i < n; i++)
      {
        a(i) = i + 1;
        b(i) = i + 1;
      }
      // Modify one element in b.
      b(0) = 999;
      Rodin::IndexArrayEquality eq;
      EXPECT_FALSE(eq(a, b)) << "Inequality not detected for array size " << n;
    }
  }

  // -----------------------------------------------------------------------------
  // IndexArraySymmetricEquality Tests
  // -----------------------------------------------------------------------------
  TEST(IndexArraySymmetricEqualityTest, SingleElement)
  {
    Rodin::IndexArray a(1), b(1);
    a << 5;
    b << 5;

    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));

    b << 6;
    EXPECT_FALSE(seq(a, b));
  }

  TEST(IndexArraySymmetricEqualityTest, TwoElementsOrderIndependent)
  {
    Rodin::IndexArray a(2), b(2);
    a << 1, 2;
    b << 1, 2;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));

    b << 2, 1;
    EXPECT_TRUE(seq(a, b));

    b << 2, 3;
    EXPECT_FALSE(seq(a, b));
  }

  TEST(IndexArraySymmetricEqualityTest, ThreeOrMoreElementsPermutation)
  {
    Rodin::IndexArray a(3), b(3);
    a << 3, 1, 2;
    b << 1, 2, 3;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));

    b << 1, 2, 4;
    EXPECT_FALSE(seq(a, b));
  }

  TEST(IndexArraySymmetricEqualityTest, DifferentSizes)
  {
    // Arrays with different sizes (non-empty) should return false.
    Rodin::IndexArray a(3), b(4);
    a << 1, 2, 3;
    b << 1, 2, 3, 4;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_FALSE(seq(a, b));
  }

  // Additional tests for arrays with sizes 6 through 10.
  TEST(IndexArraySymmetricEqualityLargeTest, PermutationEqualityLarge)
  {
    for (int n = 6; n <= 10; n++)
    {
      Rodin::IndexArray a(n), b(n);
      for (int i = 0; i < n; i++)
      {
        a(i) = i + 1;
      }
      // b is a permutation of a: reverse the order.
      for (int i = 0; i < n; i++)
      {
        b(i) = a(n - i - 1);
      }
      Rodin::IndexArraySymmetricEquality seq;
      EXPECT_TRUE(seq(a, b)) << "Symmetric equality failed for array size " << n;
    }
  }

  TEST(IndexArraySymmetricEqualityLargeTest, DifferentElementsLarge)
  {
    for (int n = 6; n <= 10; n++)
    {
      Rodin::IndexArray a(n), b(n);
      for (int i = 0; i < n; i++)
      {
        a(i) = i + 1;
        b(i) = i + 1;
      }
      // Modify one element in b so that b is no longer a permutation of a.
      b(0) = 999;
      Rodin::IndexArraySymmetricEquality seq;
      EXPECT_FALSE(seq(a, b)) << "Symmetric inequality failed for array size " << n;
    }
  }

  // -----------------------------------------------------------------------------
  // IndexArrayHash Tests
  // -----------------------------------------------------------------------------
  TEST(IndexArrayHashTest, SameArrayProducesSameHash)
  {
    Rodin::IndexArray a(3), b(3);
    a << 4, 5, 6;
    b << 4, 5, 6;
    Rodin::IndexArrayHash hashFunc;

    size_t hashA = hashFunc(a);
    size_t hashB = hashFunc(b);
    EXPECT_EQ(hashA, hashB);
  }

  TEST(IndexArrayHashTest, OrderMattersForHash)
  {
    Rodin::IndexArray a(2), b(2);
    a << 10, 20;
    b << 20, 10;
    Rodin::IndexArrayHash hashFunc;

    size_t hashA = hashFunc(a);
    size_t hashB = hashFunc(b);
    // Although hash collisions can occur, with these values it is very unlikely.
    EXPECT_NE(hashA, hashB);
  }

  // Additional tests for arrays with sizes 6 through 10.
  TEST(IndexArrayHashLargeTest, SameArrayProducesSameHashLarge)
  {
    for (int n = 6; n <= 10; n++)
    {
      Rodin::IndexArray a(n), b(n);
      for (int i = 0; i < n; i++)
      {
        a(i) = i + 10;
        b(i) = i + 10;
      }
      Rodin::IndexArrayHash hashFunc;
      size_t hashA = hashFunc(a);
      size_t hashB = hashFunc(b);
      EXPECT_EQ(hashA, hashB) << "Hash values differ for equal arrays of size " << n;
    }
  }

  TEST(IndexArrayHashLargeTest, OrderMattersLarge)
  {
    for (int n = 6; n <= 10; n++)
    {
      Rodin::IndexArray a(n), b(n);
      for (int i = 0; i < n; i++)
      {
        a(i) = i + 1;
        b(i) = i + 1;
      }
      // Reverse b.
      for (int i = 0; i < n / 2; i++)
      {
        std::swap(b(i), b(n - i - 1));
      }
      Rodin::IndexArrayHash hashFunc;
      size_t hashA = hashFunc(a);
      size_t hashB = hashFunc(b);
      if (!(a == b).all())
      {
        EXPECT_NE(hashA, hashB) << "Order does not affect hash for array size " << n;
      }
    }
  }

  // -----------------------------------------------------------------------------
  // IndexArraySymmetricHash Tests
  // -----------------------------------------------------------------------------
  TEST(IndexArraySymmetricHashTest, SingleElement)
  {
    Rodin::IndexArray a(1), b(1);
    a << 7;
    b << 7;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));
  }

  TEST(IndexArraySymmetricHashTest, TwoElementsOrderIndependent)
  {
    Rodin::IndexArray a(2), b(2);
    a << 100, 200;
    b << 200, 100;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));
  }

  TEST(IndexArraySymmetricHashTest, ThreeOrMoreElementsOrderIndependent)
  {
    Rodin::IndexArray a(4), b(4);
    a << 1, 2, 3, 4;
    b << 4, 3, 2, 1;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));

    b << 4, 3, 2, 5;
    EXPECT_NE(symHash(a), symHash(b));
  }

  // Additional tests for arrays with sizes 6 through 10.
  TEST(IndexArraySymmetricHashLargeTest, PermutationHashLarge)
  {
    for (int n = 6; n <= 10; n++)
    {
      Rodin::IndexArray a(n), b(n);
      for (int i = 0; i < n; i++)
      {
        a(i) = i + 100;
        b(i) = i + 100;
      }
      // Permute b by reversing its order.
      for (int i = 0; i < n / 2; i++)
      {
        std::swap(b(i), b(n - i - 1));
      }
      Rodin::IndexArraySymmetricHash symHash;
      EXPECT_EQ(symHash(a), symHash(b))
        << "Symmetric hash does not match for permutation of array size " << n;
    }
  }

  TEST(IndexArraySymmetricHashLargeTest, DifferentElementsHashLarge)
  {
    for (int n = 6; n <= 10; n++)
    {
      Rodin::IndexArray a(n), b(n);
      for (int i = 0; i < n; i++)
      {
        a(i) = i + 50;
        b(i) = i + 50;
      }
      b(0) = 999;
      Rodin::IndexArraySymmetricHash symHash;
      EXPECT_NE(symHash(a), symHash(b))
        << "Symmetric hash failed to detect difference for array size " << n;
    }
  }
}


