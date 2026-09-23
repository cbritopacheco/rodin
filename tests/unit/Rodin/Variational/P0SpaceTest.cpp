/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief P0 value types and cellwise degrees of freedom. */

#include <gtest/gtest.h>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  class P0SpaceTest : public ::testing::TestWithParam<Polytope::Type> {};

  TEST_P(P0SpaceTest, ScalarRealAndComplexAreCellwiseConstant)
  {
    const auto geometry = GetParam();
    const size_t dim = Polytope::Traits(geometry).getDimension();
    LocalMesh mesh;
    if (dim == 1) mesh = LocalMesh::UniformGrid(geometry, {3});
    if (dim == 2) mesh = LocalMesh::UniformGrid(geometry, {3, 3});
    if (dim == 3) mesh = LocalMesh::UniformGrid(geometry, {3, 3, 3});

    P0 realSpace(mesh);
    ComplexP0<LocalMesh> complexSpace(mesh);
    ASSERT_EQ(realSpace.getSize(), mesh.getCellCount());
    ASSERT_EQ(complexSpace.getSize(), mesh.getCellCount());

    GridFunction realField(realSpace);
    GridFunction complexField(complexSpace);
    for (Index i = 0; i < mesh.getCellCount(); ++i)
    {
      realField[i] = Real(i + 1);
      complexField[i] = Complex(Real(i + 1), -Real(i + 2));
      EXPECT_EQ(realSpace.getDOFs(dim, i)[0], i);
      EXPECT_EQ(complexSpace.getDOFs(dim, i)[0], i);
      const auto cell = mesh.getCell(i);
      const Point point(*cell, Polytope::Traits(geometry).getCentroid());
      EXPECT_NEAR(realField(point), Real(i + 1), 1e-12);
      EXPECT_NEAR(std::abs(complexField(point)
        - Complex(Real(i + 1), -Real(i + 2))), 0, 1e-12);
    }
  }

  TEST_P(P0SpaceTest, RealAndComplexVectorsHaveIndependentCellDOFs)
  {
    const auto geometry = GetParam();
    const size_t dim = Polytope::Traits(geometry).getDimension();
    LocalMesh mesh;
    if (dim == 1) mesh = LocalMesh::UniformGrid(geometry, {3});
    if (dim == 2) mesh = LocalMesh::UniformGrid(geometry, {3, 3});
    if (dim == 3) mesh = LocalMesh::UniformGrid(geometry, {3, 3, 3});

    P0<Math::SpatialVector<Real>, LocalMesh> realSpace(mesh, 2);
    P0<Math::SpatialVector<Complex>, LocalMesh> complexSpace(mesh, 2);
    ASSERT_EQ(realSpace.getSize(), 2 * mesh.getCellCount());
    ASSERT_EQ(complexSpace.getSize(), 2 * mesh.getCellCount());
    GridFunction realField(realSpace);
    GridFunction complexField(complexSpace);

    for (Index i = 0; i < mesh.getCellCount(); ++i)
    {
      for (Index component = 0; component < 2; ++component)
      {
        const Index global = 2 * i + component;
        EXPECT_EQ(realSpace.getGlobalIndex({dim, i}, component), global);
        EXPECT_EQ(complexSpace.getGlobalIndex({dim, i}, component), global);
        realField[global] = Real(global + 1);
        complexField[global] = Complex(Real(global + 1), Real(component + 1));
      }
      const auto cell = mesh.getCell(i);
      const Point point(*cell, Polytope::Traits(geometry).getCentroid());
      const auto realValue = realField(point);
      const auto complexValue = complexField(point);
      for (Index component = 0; component < 2; ++component)
      {
        const Index global = 2 * i + component;
        EXPECT_NEAR(realValue(component), Real(global + 1), 1e-12);
        EXPECT_NEAR(std::abs(complexValue(component)
          - Complex(Real(global + 1), Real(component + 1))), 0, 1e-12);
      }
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, P0SpaceTest,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge));
}
