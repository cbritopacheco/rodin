#include <gtest/gtest.h>

#include <Rodin/FlatSet.h>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  //=== IndexSetEquality Tests ===============================================

  TEST(Rodin_FlatSet, IndexSetEquality_BothEmpty)
  {
    IndexSet a, b;
    IndexSetEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  TEST(Rodin_FlatSet, IndexSetEquality_DifferentSizes)
  {
    IndexSet a{1, 2, 3};
    IndexSet b{4, 5};
    IndexSetEquality eq;
    EXPECT_FALSE(eq(a, b));
  }

  TEST(Rodin_FlatSet, IndexSetEquality_SameSizeDifferentElements)
  {
    IndexSet a{1, 2, 3};
    IndexSet b{1, 2, 4};
    IndexSetEquality eq;
    EXPECT_FALSE(eq(a, b));
  }

  TEST(Rodin_FlatSet, IndexSetEquality_SameElements)
  {
    IndexSet a{1, 2, 3};
    IndexSet b{1, 2, 3};
    IndexSetEquality eq;
    EXPECT_TRUE(eq(a, b));
  }

  TEST(Rodin_FlatSet, IndexSetEquality_SingleElement)
  {
    IndexSet a{42};
    IndexSet b{42};
    IndexSet c{43};
    IndexSetEquality eq;
    
    EXPECT_TRUE(eq(a, b));
    EXPECT_FALSE(eq(a, c));
  }

  TEST(Rodin_FlatSet, IndexSetEquality_OrderMatters)
  {
    // Note: FlatSet maintains sorted order, so this tests the underlying comparison
    IndexSet a{1, 2, 3};
    IndexSet b{3, 1, 2};  // Will be sorted to {1, 2, 3}
    IndexSetEquality eq;
    EXPECT_TRUE(eq(a, b));  // Should be equal because FlatSet sorts elements
  }

  TEST(Rodin_FlatSet, IndexSetEquality_LargeSet)
  {
    IndexSet a, b;
    for (Index i = 0; i < 1000; ++i)
    {
      a.insert(i);
      b.insert(i);
    }
    
    IndexSetEquality eq;
    EXPECT_TRUE(eq(a, b));
    
    // Modify one element
    b.erase(500);
    b.insert(9999);
    EXPECT_FALSE(eq(a, b));
  }

  //=== IndexSetHash Tests ===================================================

  TEST(Rodin_FlatSet, IndexSetHash_EmptySet)
  {
    IndexSet empty;
    IndexSetHash hasher;
    size_t hash = hasher(empty);
    
    // Hash of empty set should be deterministic
    EXPECT_EQ(hash, hasher(empty));
  }

  TEST(Rodin_FlatSet, IndexSetHash_SingleElement)
  {
    IndexSet a{42};
    IndexSet b{42};
    IndexSet c{43};
    
    IndexSetHash hasher;
    
    EXPECT_EQ(hasher(a), hasher(b));  // Same elements should have same hash
    EXPECT_NE(hasher(a), hasher(c));  // Different elements should have different hash (usually)
  }

  TEST(Rodin_FlatSet, IndexSetHash_MultipleElements)
  {
    IndexSet a{1, 2, 3};
    IndexSet b{1, 2, 3};
    IndexSet c{1, 2, 4};
    
    IndexSetHash hasher;
    
    EXPECT_EQ(hasher(a), hasher(b));
    EXPECT_NE(hasher(a), hasher(c));
  }

  TEST(Rodin_FlatSet, IndexSetHash_OrderIndependent)
  {
    IndexSet a{1, 2, 3};
    IndexSet b{3, 1, 2};  // Will be sorted to {1, 2, 3}
    
    IndexSetHash hasher;
    EXPECT_EQ(hasher(a), hasher(b));  // Should be equal since FlatSet sorts
  }

  TEST(Rodin_FlatSet, IndexSetHash_Consistency)
  {
    IndexSet set{5, 10, 15, 20};
    IndexSetHash hasher;
    
    size_t hash1 = hasher(set);
    size_t hash2 = hasher(set);
    size_t hash3 = hasher(set);
    
    EXPECT_EQ(hash1, hash2);
    EXPECT_EQ(hash2, hash3);
  }

  TEST(Rodin_FlatSet, IndexSetHash_DifferentSets)
  {
    IndexSetHash hasher;
    
    IndexSet set1{1};
    IndexSet set2{1, 2};
    IndexSet set3{1, 2, 3};
    IndexSet set4{2, 3, 4};
    
    size_t hash1 = hasher(set1);
    size_t hash2 = hasher(set2);
    size_t hash3 = hasher(set3);
    size_t hash4 = hasher(set4);
    
    // These should all be different (though hash collisions are theoretically possible)
    EXPECT_NE(hash1, hash2);
    EXPECT_NE(hash1, hash3);
    EXPECT_NE(hash1, hash4);
    EXPECT_NE(hash2, hash3);
    EXPECT_NE(hash2, hash4);
    EXPECT_NE(hash3, hash4);
  }

  TEST(Rodin_FlatSet, IndexSetHash_LargeSet)
  {
    IndexSet largeSet;
    for (Index i = 0; i < 100; ++i)
    {
      largeSet.insert(i * 7);  // Insert some pattern
    }
    
    IndexSetHash hasher;
    size_t hash = hasher(largeSet);
    
    // Test consistency
    EXPECT_EQ(hash, hasher(largeSet));
    
    // Test that modification changes hash
    IndexSet modifiedSet = largeSet;
    modifiedSet.insert(999999);
    EXPECT_NE(hash, hasher(modifiedSet));
  }

  //=== Integration Tests ====================================================

  TEST(Rodin_FlatSet, EqualityAndHashConsistency)
  {
    IndexSet a{10, 20, 30};
    IndexSet b{10, 20, 30};
    IndexSet c{10, 20, 31};
    
    IndexSetEquality eq;
    IndexSetHash hasher;
    
    // If two sets are equal, their hashes should be equal
    EXPECT_TRUE(eq(a, b));
    EXPECT_EQ(hasher(a), hasher(b));
    
    // If two sets are not equal, their hashes should usually be different
    EXPECT_FALSE(eq(a, c));
    EXPECT_NE(hasher(a), hasher(c));
  }

  TEST(Rodin_FlatSet, FunctorCallOperator)
  {
    // Test that functors can be used with parentheses operator
    IndexSet set1{1, 2, 3};
    IndexSet set2{1, 2, 3};
    
    IndexSetEquality equality;
    IndexSetHash hash;
    
    bool isEqual = equality(set1, set2);
    size_t hashValue = hash(set1);
    
    EXPECT_TRUE(isEqual);
    EXPECT_GT(hashValue, 0);  // Hash should be non-zero for non-empty set
  }
}