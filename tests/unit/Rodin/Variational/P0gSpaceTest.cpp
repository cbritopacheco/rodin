/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Real and complex global constant spaces on every cell type. */

#include <gtest/gtest.h>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  class P0gSpaceTest : public ::testing::TestWithParam<Polytope::Type> {};

  LocalMesh makeP0gMesh(Polytope::Type geometry)
  {
    const size_t dim = Polytope::Traits(geometry).getDimension();
    if (dim == 1) return LocalMesh::UniformGrid(geometry, {3});
    if (dim == 2) return LocalMesh::UniformGrid(geometry, {3, 3});
    return LocalMesh::UniformGrid(geometry, {3, 3, 3});
  }

  TEST_P(P0gSpaceTest, ScalarsShareOneComplexOrRealDOF)
  {
    const auto geometry = GetParam();
    auto mesh = makeP0gMesh(geometry);
    const size_t dim = mesh.getDimension();
    P0g realSpace(mesh);
    ComplexP0g<LocalMesh> complexSpace(mesh);
    ASSERT_EQ(realSpace.getSize(), 1);
    ASSERT_EQ(complexSpace.getSize(), 1);
    GridFunction realField(realSpace);
    GridFunction complexField(complexSpace);
    realField[0] = 2.5;
    complexField[0] = Complex(2.5, -1.25);
    for (Index i = 0; i < mesh.getCellCount(); ++i)
    {
      EXPECT_EQ(realSpace.getDOFs(dim, i)[0], 0);
      EXPECT_EQ(complexSpace.getDOFs(dim, i)[0], 0);
      const auto cell = mesh.getCell(i);
      const Point point(*cell, Polytope::Traits(geometry).getCentroid());
      EXPECT_NEAR(realField(point), 2.5, 1e-12);
      EXPECT_NEAR(std::abs(complexField(point) - Complex(2.5, -1.25)),
        0, 1e-12);
    }
  }

  TEST_P(P0gSpaceTest, VectorsShareOneDOFPerComponent)
  {
    const auto geometry = GetParam();
    auto mesh = makeP0gMesh(geometry);
    const size_t dim = mesh.getDimension();
    P0g<Math::SpatialVector<Real>, LocalMesh> realSpace(mesh, 2);
    P0g<Math::SpatialVector<Complex>, LocalMesh> complexSpace(mesh, 2);
    ASSERT_EQ(realSpace.getSize(), 2);
    ASSERT_EQ(complexSpace.getSize(), 2);
    GridFunction realField(realSpace);
    GridFunction complexField(complexSpace);
    realField[0] = 1.25;
    realField[1] = -2.5;
    complexField[0] = Complex(1.25, 2);
    complexField[1] = Complex(-2.5, -3);
    for (Index i = 0; i < mesh.getCellCount(); ++i)
    {
      const auto& realDOFs = realSpace.getDOFs(dim, i);
      const auto& complexDOFs = complexSpace.getDOFs(dim, i);
      ASSERT_EQ(realDOFs.size(), 2);
      ASSERT_EQ(complexDOFs.size(), 2);
      EXPECT_EQ(realDOFs[0], 0);
      EXPECT_EQ(realDOFs[1], 1);
      EXPECT_EQ(complexDOFs[0], 0);
      EXPECT_EQ(complexDOFs[1], 1);
      const auto cell = mesh.getCell(i);
      const Point point(*cell, Polytope::Traits(geometry).getCentroid());
      const auto realValue = realField(point);
      const auto complexValue = complexField(point);
      EXPECT_NEAR(realValue(0), 1.25, 1e-12);
      EXPECT_NEAR(realValue(1), -2.5, 1e-12);
      EXPECT_NEAR(std::abs(complexValue(0) - Complex(1.25, 2)), 0, 1e-12);
      EXPECT_NEAR(std::abs(complexValue(1) - Complex(-2.5, -3)), 0, 1e-12);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, P0gSpaceTest,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge));
}
