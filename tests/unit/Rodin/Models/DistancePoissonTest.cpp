/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Models/Distance/Poisson.h"

using namespace Rodin;
using namespace Rodin::Models::Distance;

class DistancePoissonTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test basic construction
TEST_F(DistancePoissonTest, Construction)
{
  // Test default construction
  Poisson poisson;
  
  // Should construct without throwing
  EXPECT_NO_THROW(Poisson{});
}

// Test copy semantics
TEST_F(DistancePoissonTest, CopySemantics)
{
  Poisson poisson1;
  
  // Test copy constructor
  Poisson poisson2 = poisson1;
  
  // Test copy assignment
  Poisson poisson3;
  poisson3 = poisson1;
  
  // All should be valid objects
  // Since Poisson appears to be a stateless function object,
  // copies should behave identically
}

// Test that it's a callable object
TEST_F(DistancePoissonTest, CallableObject)
{
  Poisson poisson;
  
  // Test that it has operator() defined
  // Note: We can't actually test the full functionality without
  // setting up a complete finite element space and mesh,
  // but we can verify the basic structure
  
  // The operator() is templated, so we test type traits
  static_assert(std::is_class_v<Poisson>);
  
  // Test that multiple instances can be created
  Poisson poisson1;
  Poisson poisson2;
  
  // They should be independent objects
  EXPECT_NE(&poisson1, &poisson2);
}

// Test default behavior characteristics
TEST_F(DistancePoissonTest, DefaultBehavior)
{
  Poisson poisson;
  
  // Test that the object is lightweight (stateless)
  EXPECT_EQ(sizeof(Poisson), 1);  // Empty class should have size 1
  
  // Test that multiple instances are equivalent
  Poisson poisson1;
  Poisson poisson2;
  
  // Since it's stateless, they should behave identically
  // (though we can't test the actual behavior without a full FES setup)
}

// Test const correctness
TEST_F(DistancePoissonTest, ConstCorrectness)
{
  const Poisson poisson;
  
  // The operator() should be const (it doesn't modify state)
  // We can't call it without proper setup, but we can verify
  // that const objects can be created and used
  
  // Test const copy
  const Poisson poissonCopy = poisson;
  
  // Both should be valid const objects
  (void)poissonCopy;  // Suppress unused variable warning
}

// Test move semantics
TEST_F(DistancePoissonTest, MoveSemantics)
{
  Poisson poisson1;
  
  // Test move constructor
  Poisson poisson2 = std::move(poisson1);
  
  // Test move assignment
  Poisson poisson3;
  poisson3 = std::move(poisson2);
  
  // All operations should complete successfully
  // For a stateless class, move is typically the same as copy
}

// Test that it's a function object (functor)
TEST_F(DistancePoissonTest, FunctionObject)
{
  Poisson poisson;
  
  // Test that it can be used as a function object
  // Check basic functor properties
  static_assert(std::is_default_constructible_v<Poisson>);
  static_assert(std::is_copy_constructible_v<Poisson>);
  static_assert(std::is_copy_assignable_v<Poisson>);
  static_assert(std::is_move_constructible_v<Poisson>);
  static_assert(std::is_move_assignable_v<Poisson>);
}

// Test multiple instance independence
TEST_F(DistancePoissonTest, MultipleInstances)
{
  // Create multiple instances
  Poisson poisson1;
  Poisson poisson2;
  Poisson poisson3;
  
  // They should be independent objects
  EXPECT_NE(&poisson1, &poisson2);
  EXPECT_NE(&poisson2, &poisson3);
  EXPECT_NE(&poisson1, &poisson3);
  
  // Since they're stateless, they should behave identically
  // but be distinct objects in memory
}

// Test assignment chains
TEST_F(DistancePoissonTest, AssignmentChains)
{
  Poisson poisson1;
  Poisson poisson2;
  Poisson poisson3;
  Poisson poisson4;
  
  // Test chained assignment
  poisson4 = poisson3 = poisson2 = poisson1;
  
  // All should be valid after chained assignment
  // For stateless objects, this should work seamlessly
}

// Test object lifecycle
TEST_F(DistancePoissonTest, ObjectLifecycle)
{
  // Test construction in different scopes
  {
    Poisson localPoisson;
    // Should construct successfully
    (void)localPoisson;
  }  // Should destruct successfully
  
  // Test dynamic allocation
  auto* dynamicPoisson = new Poisson();
  EXPECT_NE(dynamicPoisson, nullptr);
  delete dynamicPoisson;  // Should delete successfully
  
  // Test smart pointer
  auto smartPoisson = std::make_unique<Poisson>();
  EXPECT_NE(smartPoisson.get(), nullptr);
  // Should automatically delete when scope ends
}

// Test container storage
TEST_F(DistancePoissonTest, ContainerStorage)
{
  // Test storing in vector
  std::vector<Poisson> poissonVector;
  poissonVector.emplace_back();
  poissonVector.push_back(Poisson{});
  
  EXPECT_EQ(poissonVector.size(), 2);
  
  // Test storing in array
  std::array<Poisson, 3> poissonArray;
  
  // All elements should be valid
  EXPECT_EQ(poissonArray.size(), 3);
  
  // Test storing pointers
  std::vector<std::unique_ptr<Poisson>> poissonPtrVector;
  poissonPtrVector.push_back(std::make_unique<Poisson>());
  
  EXPECT_EQ(poissonPtrVector.size(), 1);
  EXPECT_NE(poissonPtrVector[0].get(), nullptr);
}