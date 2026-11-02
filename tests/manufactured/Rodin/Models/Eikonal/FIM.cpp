/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Eikonal/FIM.h"
#include "Rodin/Geometry/SubMesh.h"
#include "Rodin/Variational/MatrixFunction.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief Native multidimensional updates propagate an exact plane over multiple cells.
TEST(FIMManufactured, NativeGeometryPlanarWaves)
{
  for (auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {3})
      : dimension == 2         ? LocalMesh::UniformGrid(geometry, {3, 3})
                               : LocalMesh::UniformGrid(geometry, {3, 3, 3});
    mesh.scale(Real(0.5));
    mesh.getConnectivity().compute(dimension, 0);
    mesh.getConnectivity().compute(0, dimension);
    P1 space(mesh);
    GridFunction arrival(space);
    IndexVector seeds;
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      if (vertex->getCoordinates()(0) == 0)
        seeds.push_back(vertex->getIndex());
    }
    ASSERT_FALSE(seeds.empty());
    Eikonal::FIM solver(arrival, [](const Point&) { return Real(2); });
    solver.seed(seeds).solve();
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      const Real computed = arrival[space.getDOFs(0, vertex->getIndex())[0]];
      const Real exact = vertex->getCoordinates()(0) / 2;
      ASSERT_TRUE(std::isfinite(computed));
      EXPECT_NEAR(computed, exact, 2e-5);
    }
  }
}

/// @brief Rotated constant-tensor point-source errors decrease under native mesh refinement.
TEST(FIMManufactured, AnisotropicPointSourceRefinement)
{
  Math::SpatialMatrix<Real> tensor(2, 2);
  tensor(0, 0) = 4;
  tensor(0, 1) = tensor(1, 0) = 1;
  tensor(1, 1) = 1;
  MatrixFunction coefficient(tensor);
  for (auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral})
  {
    Real previousError = std::numeric_limits<Real>::infinity();
    for (size_t intervals : {2u, 4u, 8u})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      SCOPED_TRACE(intervals);
      auto mesh = LocalMesh::UniformGrid(geometry, {intervals + 1, intervals + 1});
      mesh.scale(Real(1) / intervals);
      mesh.getConnectivity().compute(2, 0);
      mesh.getConnectivity().compute(0, 2);
      P1 space(mesh);
      GridFunction arrival(space);
      IndexVector seeds;
      for (auto vertex = mesh.getVertex(); vertex; ++vertex)
      {
        if (vertex->getCoordinates().squaredNorm() == 0)
          seeds.push_back(vertex->getIndex());
      }
      ASSERT_EQ(seeds.size(), 1);
      Eikonal::FIM solver(arrival, coefficient);
      solver.setLocalTolerance(1e-6).seed(seeds).solve();
      Real meanError = 0;
      for (auto vertex = mesh.getVertex(); vertex; ++vertex)
      {
        const auto& point = vertex->getCoordinates();
        // M^{-1} = [1,-1; -1,4]/3, independently of the solver's factorization.
        const Real exact = std::sqrt(
          (point(0) * point(0) - 2 * point(0) * point(1) + 4 * point(1) * point(1)) / 3);
        const Real computed = arrival[space.getDOFs(0, vertex->getIndex())[0]];
        ASSERT_TRUE(std::isfinite(computed));
        EXPECT_GE(computed, exact - 1e-8);
        meanError += std::abs(computed - exact);
      }
      meanError /= mesh.getVertexCount();
      RecordProperty("mean_error_" + std::to_string(static_cast<int>(geometry)) + "_" +
          std::to_string(intervals),
        meanError);
      // Integrate the field error over the fixed physical domain. A changing
      // set of nodal samples overweights the exact source-cell stencil on
      // very coarse meshes and need not give a decreasing nodal average.
      Real integratedError = 0;
      const auto& quadrature = QF::QuadratureFormula::get(6, geometry);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        for (size_t sample = 0; sample < quadrature.getSize(); ++sample)
        {
          const Point point(*cell, quadrature.getPoint(sample));
          const auto& coordinates = point.getPhysicalCoordinates();
          const Real exact = std::sqrt(
            (coordinates(0) * coordinates(0) - 2 * coordinates(0) * coordinates(1) +
              4 * coordinates(1) * coordinates(1)) /
            3);
          integratedError += quadrature.getWeight(sample) * point.getDistortion() *
            std::abs(arrival(point) - exact);
        }
      }
      RecordProperty("integrated_error_" + std::to_string(static_cast<int>(geometry)) +
          "_" + std::to_string(intervals),
        integratedError);
      EXPECT_LT(integratedError, previousError);
      previousError = integratedError;
    }
  }
}
