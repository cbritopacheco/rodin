/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Kokkos.h>
#include <Rodin/Kokkos/Solver.h>

#include "../../../unit/Rodin/Kokkos/TestContext.h"

namespace Rodin::Tests::Manufactured
{
  using namespace Rodin::Variational;
  using Polytope = Geometry::Polytope;
  using Mesh = Geometry::Mesh<::Rodin::Tests::Unit::TestContext>;

  /**
   * @brief Solves an advection-diffusion-reaction equation with exact linear field.
   *
   * For @f$ u=1+x_0 @f$ and @f$ \beta=e_0 @f$, the source of
   * @f$ -\Delta u+\beta\cdot\nabla u+u=f @f$ is @f$ f=2+x_0 @f$.
   * On segments the advection coefficient is zero and @f$ f=1+x_0 @f$.
   * The field is reproduced by P1 on every supported reference geometry.
   * @param geometry Reference geometry tested.
   */
  void checkBiCGSTABManufactured(Polytope::Type geometry)
  {
    using HostMesh = Geometry::Mesh<Context::Local>;
    HostMesh source;
    if (geometry == Polytope::Type::Segment)
      source = HostMesh::UniformGrid(geometry, {9});
    else if (geometry == Polytope::Type::Triangle ||
      geometry == Polytope::Type::Quadrilateral)
      source = HostMesh::UniformGrid(geometry, {5, 5});
    else
      source = HostMesh::UniformGrid(geometry, {4, 4, 4});
    const size_t dimension = source.getDimension();
    source.getConnectivity().compute(dimension - 1, dimension);
    const Mesh mesh(::Rodin::Tests::Unit::TestContext{}, source);
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    RealFunction exact(KOKKOS_LAMBDA(const Rodin::Kokkos::Geometry::PointView& p)->Real {
      return 1 + p.x();
    });
    const bool advective = mesh.getSpaceDimension() > 1;
    RealFunction force(KOKKOS_LAMBDA(const Rodin::Kokkos::Geometry::PointView& p)->Real {
      return (advective ? 2 : 1) + p.x();
    });
    const auto solve = [&](const auto& beta) {
      Problem problem(u, v);
      if constexpr (std::is_same_v<std::decay_t<decltype(beta)>, std::nullptr_t>)
      {
        problem = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(force, v) +
          DirichletBC(u, exact);
      }
      else
      {
        problem = Integral(Grad(u), Grad(v)) + Integral(Dot(beta, Grad(u)), v) +
          Integral(u, v) - Integral(force, v) + DirichletBC(u, exact);
      }
      Solver::BiCGSTAB solver(problem);
      solver.setTolerance(1e-12).setMaxIterations(2000).solve();
      ASSERT_TRUE(solver.success());
      GridFunction expected(fes);
      expected = exact;
      const auto coefficients = expected.getData().toHost();
      const auto actual = u.getSolution().getData().toHost();
      EXPECT_LE((actual - coefficients).norm() / coefficients.norm(), 1e-10);
      const auto& system = problem.getLinearSystem();
      const auto a = system.getOperator().toHost();
      const auto b = system.getVector().toHost();
      const auto x = system.getSolution().toHost();
      const Real residual = (b - a * x).norm() / b.norm();
      EXPECT_LE(residual, 1e-12);
      EXPECT_NEAR(solver.getError(), residual, 1e-13);
    };
    if (mesh.getSpaceDimension() == 1)
      solve(nullptr);
    else if (mesh.getSpaceDimension() == 2)
      solve(VectorFunction{Real(1), Real(0)});
    else
      solve(VectorFunction{Real(1), Real(0), Real(0)});
  }

  class KokkosBiCGSTABManufactured : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(KokkosBiCGSTABManufactured, LinearAdvectionDiffusionReaction)
  {
    checkBiCGSTABManufactured(GetParam());
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, KokkosBiCGSTABManufactured,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge));
}
