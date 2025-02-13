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

  TEST(IndexArrayEquality, Size4Equal)
  {
    Rodin::IndexArray a(4), b(4);
    a << 1, 2, 3, 4;
    b << 1, 2, 3, 4;
    Rodin::IndexArrayEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  TEST(IndexArrayEquality, Size5Different)
  {
    Rodin::IndexArray a(5), b(5);
    a << 5, 6, 7, 8, 9;
    b << 5, 6, 7, 8, 10;
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

  TEST(IndexArraySymmetricEquality, TwoElementsOrderIndependent)
  {
    Rodin::IndexArray a(2), b(2);
    a << 10, 20;
    b << 20, 10;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));
  }

  TEST(IndexArraySymmetricEquality, TwoElementsDifferent)
  {
    Rodin::IndexArray a(2), b(2);
    a << 10, 20;
    b << 10, 30;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_FALSE(seq(a, b));
  }

  TEST(IndexArraySymmetricEquality, Size4PermutationEqual)
  {
    Rodin::IndexArray a(4), b(4);
    a << 1, 2, 3, 4;
    // b is a permutation of a (reversed order)
    b << 4, 3, 2, 1;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));
  }

  TEST(IndexArraySymmetricEquality, Size5PermutationWithDuplicatesEqual)
  {
    Rodin::IndexArray a(5), b(5);
    a << 2, 2, 3, 4, 5;
    // b is a permutation of a (order changed)
    b << 5, 4, 3, 2, 2;
    Rodin::IndexArraySymmetricEquality seq;
    EXPECT_TRUE(seq(a, b));
  }

  TEST(IndexArraySymmetricEquality, Size5Different)
  {
    Rodin::IndexArray a(5), b(5);
    a << 1, 2, 3, 4, 5;
    b << 5, 4, 3, 2, 6;
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
    // For empty arrays, the for_each loop does nothing so seed remains 0.
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

  TEST(IndexArrayHash, Size4EqualHash)
  {
    Rodin::IndexArray a(4), b(4);
    a << 1, 2, 3, 4;
    b << 1, 2, 3, 4;
    Rodin::IndexArrayHash hashFunc;
    EXPECT_EQ(hashFunc(a), hashFunc(b));
  }

  TEST(IndexArrayHash, Size5DifferentOrderHash)
  {
    // For IndexArrayHash the order matters.
    Rodin::IndexArray a(5), b(5);
    a << 10, 20, 30, 40, 50;
    b << 50, 40, 30, 20, 10;
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

  TEST(IndexArraySymmetricHash, TwoElementsOrderIndependentHash)
  {
    Rodin::IndexArray a(2), b(2);
    a << 11, 22;
    b << 22, 11;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));
  }

  TEST(IndexArraySymmetricHash, Size4PermutationEqualHash)
  {
    Rodin::IndexArray a(4), b(4);
    a << 1, 2, 3, 4;
    b << 4, 3, 2, 1;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));
  }

  TEST(IndexArraySymmetricHash, Size5PermutationWithDuplicatesEqualHash)
  {
    Rodin::IndexArray a(5), b(5);
    a << 2, 2, 3, 4, 5;
    b << 5, 4, 3, 2, 2;
    Rodin::IndexArraySymmetricHash symHash;
    EXPECT_EQ(symHash(a), symHash(b));
  }

  TEST(IndexArraySymmetricHash, Size5DifferentHash)
  {
    Rodin::IndexArray a(5), b(5);
    a << 1, 2, 3, 4, 5;
    b << 5, 4, 3, 2, 6;
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
    size_t originalHash = symHash(original);

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


