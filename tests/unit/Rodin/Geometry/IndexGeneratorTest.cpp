/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <vector>
#include <set>

#include "Rodin/Geometry/IndexGenerator.h"

using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ============================================================================
  // EmptyIndexGenerator Tests
  // ============================================================================

  TEST(Geometry_EmptyIndexGenerator, DefaultConstruction)
  {
    EmptyIndexGenerator gen;
    EXPECT_TRUE(gen.end());
  }

  TEST(Geometry_EmptyIndexGenerator, CopyConstruction)
  {
    EmptyIndexGenerator gen1;
    EmptyIndexGenerator gen2(gen1);
    
    EXPECT_TRUE(gen2.end());
  }

  TEST(Geometry_EmptyIndexGenerator, MoveConstruction)
  {
    EmptyIndexGenerator gen1;
    EmptyIndexGenerator gen2(std::move(gen1));
    
    EXPECT_TRUE(gen2.end());
  }

  TEST(Geometry_EmptyIndexGenerator, AlwaysAtEnd)
  {
    EmptyIndexGenerator gen;
    EXPECT_TRUE(gen.end());
  }

  TEST(Geometry_EmptyIndexGenerator, Copy)
  {
    EmptyIndexGenerator gen;
    std::unique_ptr<EmptyIndexGenerator> copy(gen.copy());
    
    EXPECT_TRUE(copy->end());
  }

  TEST(Geometry_EmptyIndexGenerator, Move)
  {
    EmptyIndexGenerator gen;
    std::unique_ptr<EmptyIndexGenerator> moved(gen.move());
    
    EXPECT_TRUE(moved->end());
  }

  // ============================================================================
  // BoundedIndexGenerator Tests
  // ============================================================================

  TEST(Geometry_BoundedIndexGenerator, Construction)
  {
    BoundedIndexGenerator gen(0, 5);
    EXPECT_FALSE(gen.end());
  }

  TEST(Geometry_BoundedIndexGenerator, EmptyRange)
  {
    BoundedIndexGenerator gen(5, 5);
    EXPECT_TRUE(gen.end());
  }

  TEST(Geometry_BoundedIndexGenerator, CopyConstruction)
  {
    BoundedIndexGenerator gen1(0, 3);
    BoundedIndexGenerator gen2(gen1);
    
    // Copy constructor sets curr to end
    EXPECT_TRUE(gen2.end());
  }

  TEST(Geometry_BoundedIndexGenerator, MoveConstruction)
  {
    BoundedIndexGenerator gen1(0, 3);
    BoundedIndexGenerator gen2(std::move(gen1));
    
    EXPECT_FALSE(gen2.end());
  }

  TEST(Geometry_BoundedIndexGenerator, Dereference)
  {
    BoundedIndexGenerator gen(10, 15);
    EXPECT_EQ(*gen, 10);
  }

  TEST(Geometry_BoundedIndexGenerator, Increment)
  {
    BoundedIndexGenerator gen(0, 3);
    
    EXPECT_EQ(*gen, 0);
    ++gen;
    EXPECT_EQ(*gen, 1);
    ++gen;
    EXPECT_EQ(*gen, 2);
    ++gen;
    EXPECT_TRUE(gen.end());
  }

  TEST(Geometry_BoundedIndexGenerator, FullIteration)
  {
    BoundedIndexGenerator gen(5, 10);
    
    std::vector<Index> indices;
    while (!gen.end())
    {
      indices.push_back(*gen);
      ++gen;
    }
    
    EXPECT_EQ(indices.size(), 5);
    EXPECT_EQ(indices[0], 5);
    EXPECT_EQ(indices[1], 6);
    EXPECT_EQ(indices[2], 7);
    EXPECT_EQ(indices[3], 8);
    EXPECT_EQ(indices[4], 9);
  }

  TEST(Geometry_BoundedIndexGenerator, Copy)
  {
    BoundedIndexGenerator gen(0, 3);
    std::unique_ptr<BoundedIndexGenerator> copy(gen.copy());
    
    // Copy should be at end
    EXPECT_TRUE(copy->end());
  }

  TEST(Geometry_BoundedIndexGenerator, Move)
  {
    BoundedIndexGenerator gen(0, 3);
    std::unique_ptr<BoundedIndexGenerator> moved(gen.move());
    
    EXPECT_FALSE(moved->end());
    EXPECT_EQ(**moved, 0);
  }

  TEST(Geometry_BoundedIndexGenerator, SingleElement)
  {
    BoundedIndexGenerator gen(42, 43);
    
    EXPECT_FALSE(gen.end());
    EXPECT_EQ(*gen, 42);
    ++gen;
    EXPECT_TRUE(gen.end());
  }

  TEST(Geometry_BoundedIndexGenerator, LargeRange)
  {
    BoundedIndexGenerator gen(0, 1000);
    
    int count = 0;
    while (!gen.end())
    {
      ++gen;
      ++count;
    }
    
    EXPECT_EQ(count, 1000);
  }

  // ============================================================================
  // VectorIndexGenerator Tests
  // ============================================================================

  TEST(Geometry_VectorIndexGenerator, EmptyVector)
  {
    std::vector<Index> indices;
    VectorIndexGenerator gen(std::move(indices));
    
    EXPECT_TRUE(gen.end());
  }

  TEST(Geometry_VectorIndexGenerator, SingleElement)
  {
    std::vector<Index> indices{42};
    VectorIndexGenerator gen(std::move(indices));
    
    EXPECT_FALSE(gen.end());
    EXPECT_EQ(*gen, 42);
    ++gen;
    EXPECT_TRUE(gen.end());
  }

  TEST(Geometry_VectorIndexGenerator, MultipleElements)
  {
    std::vector<Index> indices{10, 20, 30, 40};
    VectorIndexGenerator gen(std::move(indices));
    
    EXPECT_EQ(*gen, 10);
    ++gen;
    EXPECT_EQ(*gen, 20);
    ++gen;
    EXPECT_EQ(*gen, 30);
    ++gen;
    EXPECT_EQ(*gen, 40);
    ++gen;
    EXPECT_TRUE(gen.end());
  }

  TEST(Geometry_VectorIndexGenerator, FullIteration)
  {
    std::vector<Index> original{1, 3, 5, 7, 9};
    VectorIndexGenerator gen(std::move(original));
    
    std::vector<Index> result;
    while (!gen.end())
    {
      result.push_back(*gen);
      ++gen;
    }
    
    EXPECT_EQ(result.size(), 5);
    EXPECT_EQ(result[0], 1);
    EXPECT_EQ(result[1], 3);
    EXPECT_EQ(result[2], 5);
    EXPECT_EQ(result[3], 7);
    EXPECT_EQ(result[4], 9);
  }

  TEST(Geometry_VectorIndexGenerator, Copy)
  {
    std::vector<Index> indices{1, 2, 3};
    VectorIndexGenerator gen(std::move(indices));
    
    std::unique_ptr<VectorIndexGenerator> copy(gen.copy());
    
    // Both should iterate the same way
    EXPECT_EQ(**copy, 1);
    ++(*copy);
    EXPECT_EQ(**copy, 2);
  }

  TEST(Geometry_VectorIndexGenerator, Move)
  {
    std::vector<Index> indices{5, 10, 15};
    VectorIndexGenerator gen(std::move(indices));
    
    std::unique_ptr<VectorIndexGenerator> moved(gen.move());
    
    EXPECT_FALSE(moved->end());
    EXPECT_EQ(**moved, 5);
  }

  // ============================================================================
  // SetIndexGenerator Tests
  // ============================================================================

  TEST(Geometry_SetIndexGenerator, EmptySet)
  {
    IndexSet indices;
    SetIndexGenerator gen(std::move(indices));
    
    EXPECT_TRUE(gen.end());
  }

  TEST(Geometry_SetIndexGenerator, SingleElement)
  {
    IndexSet indices{42};
    SetIndexGenerator gen(std::move(indices));
    
    EXPECT_FALSE(gen.end());
    EXPECT_EQ(*gen, 42);
    ++gen;
    EXPECT_TRUE(gen.end());
  }

  TEST(Geometry_SetIndexGenerator, MultipleElementsSorted)
  {
    IndexSet indices{30, 10, 20, 40};  // Set will sort these
    SetIndexGenerator gen(std::move(indices));
    
    // Sets are sorted, so we should get them in order
    std::vector<Index> result;
    while (!gen.end())
    {
      result.push_back(*gen);
      ++gen;
    }
    
    EXPECT_EQ(result.size(), 4);
    EXPECT_EQ(result[0], 10);
    EXPECT_EQ(result[1], 20);
    EXPECT_EQ(result[2], 30);
    EXPECT_EQ(result[3], 40);
  }

  TEST(Geometry_SetIndexGenerator, Copy)
  {
    IndexSet indices{1, 2, 3};
    SetIndexGenerator gen(indices);
    
    std::unique_ptr<SetIndexGenerator> copy(gen.copy());
    
    EXPECT_FALSE(copy->end());
    EXPECT_EQ(**copy, 1);
  }

  TEST(Geometry_SetIndexGenerator, Move)
  {
    IndexSet indices{5, 10, 15};
    SetIndexGenerator gen(std::move(indices));
    
    std::unique_ptr<SetIndexGenerator> moved(gen.move());
    
    EXPECT_FALSE(moved->end());
    EXPECT_EQ(**moved, 5);
  }

  TEST(Geometry_SetIndexGenerator, ConstructWithConstRef)
  {
    IndexSet indices{7, 14, 21};
    SetIndexGenerator gen(indices);  // const ref constructor
    
    EXPECT_FALSE(gen.end());
    EXPECT_EQ(*gen, 7);
  }
}
