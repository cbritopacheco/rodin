/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

/**
 * @brief Manufactured solution tests for the upwind DG advection problem.
 *
 * Solves the steady linear advection equation
 * @f[
 *   \boldsymbol{\beta} \cdot \nabla u = f \quad \text{in } \Omega = [0,1]^2,
 *   \qquad u = g \quad \text{on } \partial\Omega_{\mathrm{in}},
 * @f]
 * using the DG IBP upwind formulation:
 * @f[
 *   -\int_\Omega u\,(\boldsymbol{\beta}\cdot\nabla v)\,dx
 *   + \int_{\partial\Omega} (\boldsymbol{\beta}\cdot\mathbf{n})^+\, u\, v \, ds
 *   = \int_\Omega f\, v\, dx
 *   - \int_{\partial\Omega} (\boldsymbol{\beta}\cdot\mathbf{n})^-\, g\, v \, ds.
 * @f]
 */
namespace Rodin::Tests::Manufactured::DGAdvection
{
  template <size_t M>
  class Manufactured_DGAdvection_Test : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      Mesh<Context::Local> getMesh()
      {
        Mesh mesh;
        mesh = mesh.UniformGrid(GetParam(), { M, M });
        mesh.scale(1.0 / (M - 1));
        mesh.getConnectivity().compute(1, 2);
        return mesh;
      }
  };

  using Manufactured_DGAdvection_Test_16x16 =
    Manufactured_DGAdvection_Test<16>;
  using Manufactured_DGAdvection_Test_32x32 =
    Manufactured_DGAdvection_Test<32>;
  using Manufactured_DGAdvection_Test_64x64 =
    Manufactured_DGAdvection_Test<64>;

  /**
   * @brief Advection with β=(1,0) and manufactured solution u = x(1−x)y(1−y).
   *
   * For β=(1,0) the advection PDE reduces to @f$ \partial_x u = f @f$.
   * With @f$ u = x(1-x)y(1-y) @f$ the source is @f$ f = (1-2x)y(1-y) @f$
   * and the inflow condition at @f$ x=0 @f$ is @f$ g=0 @f$ (since @f$
   * u|_{x=0}=0 @f$).  The Max/Min upwinding enforces @f$ g=0 @f$ on the
   * inflow without explicit attribute restrictions.
   */
  TEST_P(Manufactured_DGAdvection_Test_32x32, DGAdvection_HorizontalFlow)
  {
    Mesh mesh = this->getMesh();

    // Advection velocity β = (1, 0)
    Math::Vector<Real> betaVec(2);
    betaVec << 1.0, 0.0;
    VectorFunction beta(betaVec);

    // Manufactured solution and source
    auto f = (1.0 - 2.0 * F::x) * F::y * (1.0 - F::y);
    auto g = Zero();

    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);
    BoundaryNormal n(mesh);

    Problem advection(u, v);
    advection = -Integral(u, Dot(beta, Grad(v)))
              + BoundaryIntegral(Max(Dot(beta, n), Zero()) * u, v)
              - Integral(f, v)
              - BoundaryIntegral(Min(Dot(beta, n), Zero()) * g, v);

    SparseLU(advection).solve();

    auto solution = F::x * (1.0 - F::x) * F::y * (1.0 - F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real error = Integral(diff).compute();
    EXPECT_LT(error, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Advection with β=(1,1)/√2 and manufactured solution u = sin(πx)sin(πy).
   *
   * Uses a diagonal flow direction so that all four boundary segments are
   * either inflow or outflow, fully exercising the Max/Min flux splitting.
   */
  TEST_P(Manufactured_DGAdvection_Test_32x32, DGAdvection_DiagonalFlow)
  {
    Mesh mesh = this->getMesh();

    const Real pi    = Math::Constants::pi();
    const Real inv_r = 1.0 / std::sqrt(2.0);

    // Advection velocity β = (1/√2, 1/√2)
    Math::Vector<Real> betaVec(2);
    betaVec << inv_r, inv_r;
    VectorFunction beta(betaVec);

    // Manufactured solution: u = sin(πx)sin(πy)
    // β·∇u = (1/√2)(π cos(πx)sin(πy) + π sin(πx)cos(πy))
    auto solution = sin(pi * F::x) * sin(pi * F::y);
    auto f        = inv_r * pi * (cos(pi * F::x) * sin(pi * F::y)
                                + sin(pi * F::x) * cos(pi * F::y));
    auto g        = solution; // inflow BC = exact solution trace

    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);
    BoundaryNormal n(mesh);

    Problem advection(u, v);
    advection = -Integral(u, Dot(beta, Grad(v)))
              + BoundaryIntegral(Max(Dot(beta, n), Zero()) * u, v)
              - Integral(f, v)
              - BoundaryIntegral(Min(Dot(beta, n), Zero()) * g, v);

    SparseLU(advection).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real error = Integral(diff).compute();
    EXPECT_LT(error, RODIN_FUZZY_CONSTANT);
  }

  // Register parameterised tests for triangles
  INSTANTIATE_TEST_SUITE_P(Triangle, Manufactured_DGAdvection_Test_16x16,
    ::testing::Values(Polytope::Type::Triangle));

  INSTANTIATE_TEST_SUITE_P(Triangle, Manufactured_DGAdvection_Test_32x32,
    ::testing::Values(Polytope::Type::Triangle));

} // namespace Rodin::Tests::Manufactured::DGAdvection
