/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2024.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "InterfaceTestUtils.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_Jump, ConstantFunction_JumpIsZero)
  {
    // Jump of a constant across any interior face should be zero
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    RealFunction f(42.0);
    auto jump_f = Jump(f);

    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(jump_f.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Jump, ConstantVectorFunction_JumpIsZero)
  {
    // Jump of a constant vector function should be a zero vector
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    VectorFunction f{1.0, 2.0};
    auto jump_f = Jump(f);

    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    auto result = jump_f.getValue(p);
    EXPECT_NEAR(result(0), 0.0, 1e-10);
    EXPECT_NEAR(result(1), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Jump, GridFunction_ContinuousFunction_JumpIsZero)
  {
    // Jump of a continuous P1 grid function across an interior face should be zero
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    GridFunction gf(fes);
    RealFunction f(5.0);
    gf.project(f);
    auto jump_gf = Jump(gf);

    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(jump_gf.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Jump, CopyConstruction)
  {
    RealFunction f(7.0);
    auto jump_f = Jump(f);
    auto jump_copy = jump_f;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(jump_copy.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Jump, MoveConstruction)
  {
    RealFunction f(7.0);
    auto jump_f = Jump(f);
    auto jump_moved = std::move(jump_f);

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(jump_moved.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Jump, PolymorphicCopy)
  {
    RealFunction f(3.0);
    auto jump_f = Jump(f);
    std::unique_ptr<decltype(jump_f)> copy(jump_f.copy());

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(copy->getValue(p), 0.0, 1e-10);
  }
}
