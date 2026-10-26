/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Copied Krylov solvers reproduce affine harmonic fields on every cell geometry. */
#include <gtest/gtest.h>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <utility>
#include "Rodin/Assembly.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include RODIN_TEST_SOLVER_HEADER

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  class HarmonicSolverOwnership : public ::testing::TestWithParam<Polytope::Type>
  {};

  /**
   * @brief Solves two affine Dirichlet problems after the original solver is destroyed.
   *
   * The manufactured field is harmonic and belongs to every scalar P1 cell
   * space. Every cell vertex and centroid is checked, so the oracle includes
   * interior unknowns and interpolation between degrees of freedom.
   */
  TEST_P(HarmonicSolverOwnership, RepeatedAffineDirichletSolves)
  {
    const auto geometry = GetParam();
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    constexpr size_t PointsPerAxis = 4;
    auto mesh =
      LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, PointsPerAxis));
    mesh.scale(Real(1) / Real(PointsPerAxis - 1));
    mesh.getConnectivity().compute(dimension - 1, dimension);
    P1 space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    Problem problem(u, v);
    using System = std::remove_reference_t<decltype(problem.getLinearSystem())>;
    using Algorithm = Rodin::Solver::RODIN_TEST_SOLVER<System>;
    constexpr Real SolverTolerance = 1e-12;
    constexpr size_t IterationLimit = 500;
    constexpr size_t RestartDimension = 32;
    constexpr Real FieldTolerance = 1e-9;
    std::unique_ptr<Algorithm> survivor;
    {
      Algorithm original(problem);
      const auto configure = [&](auto& solver) {
        solver.setTolerance(SolverTolerance).setMaxIterations(IterationLimit);
        if constexpr (requires { solver.setS(size_t(2)); })
          solver.setS(2);
        else
          solver.setRestart(RestartDimension);
      };
      configure(original);
      Algorithm copied(original);
      Algorithm moved(std::move(copied));
      survivor.reset(moved.copy());
    }
    for (const Real phase : {Real(1), Real(2)})
    {
      SCOPED_TRACE(phase);
      const auto exact = RealFunction([phase](const Point& point) {
        Real value = 1;
        for (size_t i = 0; i < point.getPhysicalCoordinates().size(); ++i)
          value += phase * Real(i + 2) * point(i);
        return value;
      });
      problem = Integral(Grad(u), Grad(v)) + DirichletBC(u, exact);
      survivor->solve();
      ASSERT_TRUE(survivor->success());
      Real maximumError = 0;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const Polytope::Traits traits(cell->getGeometry());
        const Point centroid(*cell, traits.getCentroid());
        maximumError =
          std::max(maximumError, std::abs(u.getSolution()(centroid) - exact(centroid)));
        for (size_t local = 0; local < traits.getVertexCount(); ++local)
        {
          const Point vertex(*cell, traits.getVertex(local));
          maximumError =
            std::max(maximumError, std::abs(u.getSolution()(vertex) - exact(vertex)));
        }
      }
      EXPECT_LT(maximumError, FieldTolerance);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, HarmonicSolverOwnership,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge));
}
