/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2024.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Helper to find an interior face and create a point on it.
   */
  static std::pair<bool, Point> findInteriorFacePointAvg(Mesh<Context::Local>& mesh)
  {
    const size_t d = mesh.getDimension();
    mesh.getConnectivity().compute(d - 1, d);
    for (auto it = mesh.getFace(); it; ++it)
    {
      if (!it->isBoundary())
      {
        Math::SpatialPoint rc(1);
        rc(0) = 0.5;
        return { true, Point(*it, std::move(rc)) };
      }
    }
    return { false, Point(*(mesh.getFace()), Math::SpatialPoint(1)) };
  }

  TEST(Rodin_Variational_Average, ConstantFunction_AverageEqualsConstant)
  {
    // Average of a constant function across any face should equal the constant
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    RealFunction f(42.0);
    auto avg_f = Average(f);

    auto [found, p] = findInteriorFacePointAvg(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(avg_f.getValue(p), 42.0, 1e-10);
  }

  TEST(Rodin_Variational_Average, ConstantVectorFunction_AverageEqualsConstant)
  {
    // Average of a constant vector function should equal the same vector
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    VectorFunction f{3.0, 7.0};
    auto avg_f = Average(f);

    auto [found, p] = findInteriorFacePointAvg(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    auto result = avg_f.getValue(p);
    EXPECT_NEAR(result(0), 3.0, 1e-10);
    EXPECT_NEAR(result(1), 7.0, 1e-10);
  }

  TEST(Rodin_Variational_Average, GridFunction_ContinuousFunction_AverageEqualsValue)
  {
    // Average of a continuous P1 grid function should equal the function value
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    GridFunction gf(fes);
    RealFunction f(5.0);
    gf.project(f);
    auto avg_gf = Average(gf);

    auto [found, p] = findInteriorFacePointAvg(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(avg_gf.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Average, CopyConstruction)
  {
    RealFunction f(9.0);
    auto avg_f = Average(f);
    auto avg_copy = avg_f;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    auto [found, p] = findInteriorFacePointAvg(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(avg_copy.getValue(p), 9.0, 1e-10);
  }

  TEST(Rodin_Variational_Average, MoveConstruction)
  {
    RealFunction f(9.0);
    auto avg_f = Average(f);
    auto avg_moved = std::move(avg_f);

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    auto [found, p] = findInteriorFacePointAvg(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(avg_moved.getValue(p), 9.0, 1e-10);
  }

  TEST(Rodin_Variational_Average, PolymorphicCopy)
  {
    RealFunction f(11.0);
    auto avg_f = Average(f);
    std::unique_ptr<decltype(avg_f)> copy(avg_f.copy());

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    auto [found, p] = findInteriorFacePointAvg(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(copy->getValue(p), 11.0, 1e-10);
  }

  TEST(Rodin_Variational_Average, JumpAndAverage_Consistency)
  {
    // For a constant function, Jump=0 and Average=constant
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    RealFunction f(10.0);
    auto avg_f = Average(f);
    auto jump_f = Jump(f);

    auto [found, p] = findInteriorFacePointAvg(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(jump_f.getValue(p), 0.0, 1e-10);
    EXPECT_NEAR(avg_f.getValue(p), 10.0, 1e-10);
  }
}
