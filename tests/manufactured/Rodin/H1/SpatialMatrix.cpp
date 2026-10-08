/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Manufactured matrix reaction-diffusion solutions on every cell geometry. */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
namespace
{
  template <class FES>
  void checkSolve(const FES& fes, bool diffusion)
  {
    auto exact = MatrixFunction(fes.getRows(), fes.getColumns(), [&](const Point& p) {
      Math::SpatialMatrix<Real> value(fes.getRows(), fes.getColumns());
      for (size_t r = 0; r < fes.getRows(); ++r)
      {
        for (size_t c = 0; c < fes.getColumns(); ++c)
        {
          value(r, c) = 1 + 7 * r + c;
          if (diffusion)
            for (size_t k = 0; k < p.getCoordinates().size(); ++k)
              value(r, c) += (r + c + k + 1) * p.getCoordinates()[k];
        }
      }
      return value;
    });
    exact.setOrder(diffusion ? 1 : 0);
    TrialFunction u(fes);
    TestFunction v(fes);
    auto mass = Integral(u, v);
    auto load = Integral(exact, v);
    auto stiffness = Integral(Grad(u), Grad(v));
    // Collapsed pyramid bases require more than the polynomial degree rule.
    mass.setOrder(12);
    load.setOrder(12);
    stiffness.setOrder(12);
    Problem problem(u, v);
    if (diffusion)
      problem = stiffness + mass - load + DirichletBC(u, exact);
    else
      problem = mass - load;
    Solver::CG solver(problem);
    solver.setTolerance(1e-12).setMaxIterations(2000).solve();
    GridFunction expected(fes);
    expected = exact;
    ASSERT_EQ(u.getSolution().getData().size(), expected.getData().size());
    EXPECT_LE((u.getSolution().getData() - expected.getData()).norm(),
      1e-7 * std::max(Real(1), expected.getData().norm()));
  }
}
TEST(SpatialMatrixManufactured, EverySpaceAndGeometry)
{
  using G = Polytope::Type;
  for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
         G::Pyramid, G::Wedge, G::Hexahedron})
  {
    const size_t D = Polytope::Traits(geometry).getDimension();
    auto mesh = D == 1 ? LocalMesh::UniformGrid(geometry, {4})
      : D == 2         ? LocalMesh::UniformGrid(geometry, {3, 3})
                       : LocalMesh::UniformGrid(geometry, {3, 3, 3});
    for (size_t d = 1; d <= D; ++d)
    {
      for (size_t lower = 0; lower < d; ++lower)
      {
        mesh.getConnectivity().compute(d, lower);
        mesh.getConnectivity().compute(lower, d);
      }
    }
    for (const auto& shape : {std::pair<size_t, size_t>{2, 3}, {3, 3}})
    {
      const auto [rows, cols] = shape;
      checkSolve(P0(mesh, rows, cols), false);
      checkSolve(P0g(mesh, rows, cols), false);
      checkSolve(P1(mesh, rows, cols), true);
      checkSolve(H1(std::integral_constant<size_t, 1>{}, mesh, rows, cols), true);
      checkSolve(H1(std::integral_constant<size_t, 2>{}, mesh, rows, cols), true);
      checkSolve(H1(std::integral_constant<size_t, 3>{}, mesh, rows, cols), true);
    }
  }
}

TEST(SpatialMatrixManufactured, QuadraticReactionDiffusion)
{
  auto check = []<class FES>(const FES& fes) {
    const size_t D = fes.getMesh().getDimension();
    auto exact = MatrixFunction(2, 3, [](const Point& p) {
      Math::SpatialMatrix<Real> value(2, 3);
      const Real square = p.getCoordinates().squaredNorm();
      for (size_t r = 0; r < 2; ++r)
      {
        for (size_t c = 0; c < 3; ++c)
          value(r, c) = 1 + r + c + (1 + 3 * r + c) * square;
      }
      return value;
    });
    auto forcing = MatrixFunction(2, 3, [&](const Point& p) {
      auto value = exact.getValue(p);
      for (size_t r = 0; r < 2; ++r)
      {
        for (size_t c = 0; c < 3; ++c)
          value(r, c) -= 2 * D * (1 + 3 * r + c);
      }
      return value;
    });
    exact.setOrder(2);
    forcing.setOrder(2);
    TrialFunction u(fes);
    TestFunction v(fes);
    auto mass = Integral(u, v);
    auto stiffness = Integral(Grad(u), Grad(v));
    auto load = Integral(forcing, v);
    mass.setOrder(12);
    stiffness.setOrder(12);
    load.setOrder(12);
    Problem problem(u, v);
    problem = stiffness + mass - load + DirichletBC(u, exact);
    Solver::CG solver(problem);
    solver.setTolerance(1e-12).setMaxIterations(2000).solve();
    GridFunction expected(fes);
    expected = exact;
    EXPECT_LE((u.getSolution().getData() - expected.getData()).norm(),
      1e-7 * std::max(Real(1), expected.getData().norm()));
  };
  using G = Polytope::Type;
  for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
         G::Pyramid, G::Wedge, G::Hexahedron})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t D = Polytope::Traits(geometry).getDimension();
    auto mesh = D == 1 ? LocalMesh::UniformGrid(geometry, {2})
      : D == 2         ? LocalMesh::UniformGrid(geometry, {2, 2})
                       : LocalMesh::UniformGrid(geometry, {2, 2, 2});
    for (size_t d = 1; d <= D; ++d)
    {
      for (size_t lower = 0; lower < d; ++lower)
      {
        mesh.getConnectivity().compute(d, lower);
        mesh.getConnectivity().compute(lower, d);
      }
    }
    check(H1(std::integral_constant<size_t, 2>{}, mesh, 2, 3));
    check(H1(std::integral_constant<size_t, 3>{}, mesh, 2, 3));
  }
}
