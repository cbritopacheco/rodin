#include <random>
#include <gtest/gtest.h>

#include <Rodin/Array.h>

namespace Rodin::Tests::Unit
{
  // -----------------------------------------------------------------------------
  // IndexArrayEquality Tests (order sensitive)
  // -----------------------------------------------------------------------------
  TEST(IndexArrayEquality, BothEmptyArrays)
  {
    Rodin::IndexArray a(0), b(0);
    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  TEST(IndexArrayEquality, SingleElementEqual)
  {
    Rodin::IndexArray a(1), b(1);
    a << 42;
    b << 42;
    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  TEST(IndexArrayEquality, SingleElementDifferent)
  {
    Rodin::IndexArray a(1), b(1);
    a << 42;
    b << 43;
    Rodin::IndexArrayEquality eq;
    EXPECT_FALSE(eq(a, b));
  }

  // ----- Tests for size 2 -----
  TEST(IndexArrayEquality, Size2Equal)
  {
    Rodin::IndexArray a(2), b(2);
    a << 7, 8;
    b << 7, 8;
    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

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
  TEST(IndexArrayEquality, Size3Equal)
  {
    Rodin::IndexArray a(3), b(3);
    a << 1, 2, 3;
    b << 1, 2, 3;
    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

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

  TEST(IndexArrayEquality, DifferentSizes)
  {
    Rodin::IndexArray a(4), b(5);
    a << 1, 2, 3, 4;
    b << 1, 2, 3, 4, 5;
    Rodin::IndexArrayEquality eq;
    EXPECT_FALSE(eq(a, b));
  }

  // -----------------------------------------------------------------------------
  // IndexArraySymmetricEquality Tests (order independent)
  // -----------------------------------------------------------------------------
  TEST(IndexArraySymmetricEquality, BothEmptyArrays)
  {
    Rodin::IndexArray a(0), b(0);
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));
  }

  TEST(IndexArraySymmetricEquality, SingleElementEqual)
  {
    Rodin::IndexArray a(1), b(1);
    a << 42;
    b << 42;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));
  }

  TEST(IndexArraySymmetricEquality, SingleElementDifferent)
  {
    Rodin::IndexArray a(1), b(1);
    a << 42;
    b << 43;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_FALSE(seq(a, b));
  }

  // ----- Tests for size 2 -----
  TEST(IndexArraySymmetricEquality, Size2PermutationEqual)
  {
    Rodin::IndexArray a(2), b(2);
    a << 7, 8;
    b << 8, 7;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));
  }

  TEST(IndexArraySymmetricEquality, Size2Different)
  {
    Rodin::IndexArray a(2), b(2);
    a << 7, 8;
    b << 7, 9;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_FALSE(seq(a, b));
  }

  // ----- Tests for size 3 -----
  TEST(IndexArraySymmetricEquality, Size3PermutationEqual)
  {
    Rodin::IndexArray a(3), b(3);
    a << 1, 2, 3;
    b << 3, 1, 2;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));
  }

  TEST(IndexArraySymmetricEquality, Size3NonPermutation)
  {
    Rodin::IndexArray a(3), b(3);
    a << 1, 2, 3;
    b << 3, 2, 4;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_FALSE(seq(a, b));
  }

  TEST(IndexArraySymmetricEquality, DifferentSizes)
  {
    Rodin::IndexArray a(4), b(5);
    a << 1, 2, 3, 4;
    b << 1, 2, 3, 4, 5;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_FALSE(seq(a, b));
  }

  // -----------------------------------------------------------------------------
  // IndexArrayHash Tests (order sensitive)
  // -----------------------------------------------------------------------------
  TEST(IndexArrayHash, BothEmptyArrays)
  {
    Rodin::IndexArray a(0), b(0);
    Rodin::IndexArrayHash hashFunc;
    // For empty arrays, no elements are processed, so seed remains 0.
    EXPECT_EQ(hashFunc(a), 0);
    EXPECT_EQ(hashFunc(a), hashFunc(b));
  }

  TEST(IndexArrayHash, SingleElementEqualHash)
  {
    Rodin::IndexArray a(1), b(1);
    a << 99;
    b << 99;
    Rodin::IndexArrayHash hashFunc;
    EXPECT_EQ(hashFunc(a), hashFunc(b));
  }

  // ----- Tests for size 2 -----
  TEST(IndexArrayHash, Size2EqualHash)
  {
    Rodin::IndexArray a(2), b(2);
    a << 11, 22;
    b << 11, 22;
    Rodin::IndexArrayHash hashFunc;
    EXPECT_EQ(hashFunc(a), hashFunc(b));
  }

  TEST(IndexArrayHash, Size2DifferentOrderHash)
  {
    // Order matters, so a permutation should yield a different hash.
    Rodin::IndexArray a(2), b(2);
    a << 11, 22;
    b << 22, 11;
    Rodin::IndexArrayHash hashFunc;
    EXPECT_NE(hashFunc(a), hashFunc(b));
  }

  // ----- Tests for size 3 -----
  TEST(IndexArrayHash, Size3EqualHash)
  {
    Rodin::IndexArray a(3), b(3);
    a << 5, 6, 7;
    b << 5, 6, 7;
    Rodin::IndexArrayHash hashFunc;
    EXPECT_EQ(hashFunc(a), hashFunc(b));
  }

  TEST(IndexArrayHash, Size3DifferentOrderHash)
  {
    Rodin::IndexArray a(3), b(3);
    a << 5, 6, 7;
    b << 7, 6, 5;
    Rodin::IndexArrayHash hashFunc;
    EXPECT_NE(hashFunc(a), hashFunc(b));
  }

  // -----------------------------------------------------------------------------
  // IndexArraySymmetricHash Tests (order independent)
  // -----------------------------------------------------------------------------
  TEST(IndexArraySymmetricHash, BothEmptyArrays)
  {
    Rodin::IndexArray a(0), b(0);
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), 0);
    EXPECT_EQ(symHash(a), symHash(b));
  }

  TEST(IndexArraySymmetricHash, SingleElementEqualHash)
  {
    Rodin::IndexArray a(1), b(1);
    a << 88;
    b << 88;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));
  }

  // ----- Tests for size 2 -----
  TEST(IndexArraySymmetricHash, Size2PermutationEqualHash)
  {
    Rodin::IndexArray a(2), b(2);
    a << 11, 22;
    b << 22, 11;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));
  }

  // ----- Tests for size 3 -----
  TEST(IndexArraySymmetricHash, Size3PermutationEqualHash)
  {
    Rodin::IndexArray a(3), b(3);
    a << 5, 6, 7;
    b << 7, 5, 6;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));
  }

  TEST(IndexArraySymmetricHash, Size3DifferentHash)
  {
    Rodin::IndexArray a(3), b(3);
    a << 5, 6, 7;
    b << 7, 6, 8;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_NE(symHash(a), symHash(b));
  }

  // -----------------------------------------------------------------------------
  // More Complex / Mixed Tests
  // -----------------------------------------------------------------------------
  TEST(ComplexTest, RandomShuffleSymmetricHashConsistency)
  {
    // Create a 10-element array with duplicate values.
    const int n = 10;
    Rodin::IndexArray original(n);
    for (int i = 0; i < n; i++)
    {
      original(i) = i % 4;  // values will repeat: 0,1,2,3,0,1,...
    }
    Rodin::IndexArraySymmetricHash symHash;

    // Create a permutation of the original array using std::shuffle.
    std::vector<Rodin::Index> vec(original.data(), original.data() + original.size());
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(vec.begin(), vec.end(), g);
    Rodin::IndexArray permuted(n);
    for (int i = 0; i < n; i++)
    {
      permuted(i) = vec[i];
    }

    EXPECT_EQ(symHash(original), symHash(permuted))
      << "Symmetric hash should be invariant under permutation.";

    // Alter one element in the permuted array and check that the hash changes.
    permuted(0) = permuted(0) + 1;
    EXPECT_NE(symHash(original), symHash(permuted))
      << "Changing an element should change the symmetric hash.";
  }

  TEST(ComplexTest, MixedEqualityAndHash)
  {
    // Create two arrays of size 6 with the same elements in different orders.
    const int n = 6;
    Rodin::IndexArray a(n), b(n);
    a << 3, 1, 4, 1, 5, 9;
    b << 1, 3, 1, 9, 4, 5; // permutation of a

    Rodin::IndexArrayEquality eq;
    Rodin::IndexArraySymmetricEquality seq;
    Rodin::IndexArrayHash hashFunc;
    Rodin::IndexArraySymmetricHash symHash;

    // Order-sensitive equality should fail while symmetric equality passes.
    EXPECT_FALSE(eq(a, b));
    EXPECT_TRUE(seq(a, b));

    // Order-sensitive hash is likely different; symmetric hash must match.
    EXPECT_NE(hashFunc(a), hashFunc(b));
    EXPECT_EQ(symHash(a), symHash(b));
  }
}


