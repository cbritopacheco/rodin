#include <gtest/gtest.h>

#include <Rodin/Array.h>

namespace Rodin::Tests::Unit
{
  TEST(Rodin_IndexArray_EqualityTest, EqualArrays)
  {
    // Create two arrays with the same size and same elements.
    Rodin::IndexArray a(3), b(3);
    a << 1, 2, 3;
    b << 1, 2, 3;

    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  TEST(Rodin_IndexArray_EqualityTest, DifferentArrays)
  {
    // Create two arrays with same size but differing in at least one element.
    Rodin::IndexArray a(3), b(3);
    a << 1, 2, 3;
    b << 1, 2, 4;

    Rodin::IndexArrayEquality eq;
    EXPECT_FALSE(eq(a, b));
  }

  // -----------------------------------------------------------------------------
  // IndexArraySymmetricEquality Tests
  // -----------------------------------------------------------------------------
  TEST(Rodin_IndexArray_SymmetricEqualityTest, SingleElement)
  {
    // For size 1, order is trivial.
    Rodin::IndexArray a(1), b(1);
    a << 5;
    b << 5;

    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));

    b << 6;
    EXPECT_FALSE(seq(a, b));
  }

  TEST(Rodin_IndexArray_SymmetricEqualityTest, TwoElementsOrderIndependent)
  {
    // For size 2, order should not matter.
    Rodin::IndexArray a(2), b(2);
    a << 1, 2;
    b << 1, 2;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));

    b << 2, 1;
    EXPECT_TRUE(seq(a, b));

    // Change one element.
    b << 2, 3;
    EXPECT_FALSE(seq(a, b));
  }

  TEST(Rodin_IndexArray_SymmetricEqualityTest, ThreeOrMoreElementsPermutation)
  {
    // For size > 2, the function uses std::is_permutation.
    Rodin::IndexArray a(3), b(3);
    a << 3, 1, 2;
    b << 1, 2, 3;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));

    // Not a permutation.
    b << 1, 2, 4;
    EXPECT_FALSE(seq(a, b));
  }

  // -----------------------------------------------------------------------------
  // IndexArrayHash Tests
  // -----------------------------------------------------------------------------
  TEST(Rodin_IndexArray_HashTest, SameArrayProducesSameHash)
  {
    // Two arrays with the same ordering should have the same hash.
    Rodin::IndexArray a(3), b(3);
    a << 4, 5, 6;
    b << 4, 5, 6;
    Rodin::IndexArrayHash hashFunc;

    size_t hashA = hashFunc(a);
    size_t hashB = hashFunc(b);
    EXPECT_EQ(hashA, hashB);
  }

  TEST(Rodin_IndexArray_HashTest, OrderMattersForHash)
  {
    // For IndexArrayHash, order is important.
    Rodin::IndexArray a(2), b(2);
    a << 10, 20;
    b << 20, 10;
    Rodin::IndexArrayHash hashFunc;

    size_t hashA = hashFunc(a);
    size_t hashB = hashFunc(b);
    // Although hash collisions can occur in principle, with these small values
    // it is very likely that the hash values are different.
    EXPECT_NE(hashA, hashB);
  }

  // -----------------------------------------------------------------------------
  // IndexArraySymmetricHash Tests
  // -----------------------------------------------------------------------------
  TEST(Rodin_IndexArray_SymmetricHashTest, SingleElement)
  {
    Rodin::IndexArray a(1), b(1);
    a << 7;
    b << 7;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));
  }

  TEST(Rodin_IndexArray_SymmetricHashTest, TwoElementsOrderIndependent)
  {
    // For size 2, the symmetric hash should yield the same value regardless of order.
    Rodin::IndexArray a(2), b(2);
    a << 100, 200;
    b << 200, 100;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));
  }

  TEST(Rodin_IndexArray_SymmetricHashTest, ThreeOrMoreElementsOrderIndependent)
  {
    // For size > 2, symmetric hash should be independent of ordering.
    Rodin::IndexArray a(4), b(4);
    a << 1, 2, 3, 4;
    b << 4, 3, 2, 1;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));

    // Changing an element should change the hash.
    b << 4, 3, 2, 5;
    EXPECT_NE(symHash(a), symHash(b));
  }
}


