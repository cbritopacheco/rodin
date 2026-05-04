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
 * @brief Manufactured solution tests for the P0 upwind DG advection problem.
 *
 * Solves the steady linear advection equation
 * @f[
 *   \boldsymbol{\beta} \cdot \nabla u = f \quad \text{in } \Omega = [0,1]^2,
 *   \qquad u = g \quad \text{on } \partial\Omega_{\mathrm{in}},
 * @f]
 * using a P0 (piecewise constant) DG upwind scheme.  Since @f$ \nabla v = 0
 * @f$ for P0, the scheme is assembled purely through interface integrals:
 * @f[
 *   \sum_{E} \int_E \Bigl[(\boldsymbol{\beta}\cdot\mathbf{n})\,\{\!\{u\}\!\}
 *   + \tfrac{1}{2}|\boldsymbol{\beta}\cdot\mathbf{n}|\,[\![u]\!]\Bigr][\![v]\!]\,ds
 *   + \int_{\partial\Omega} (\boldsymbol{\beta}\cdot\mathbf{n})^+\, u\, v\, ds
 *   = \int_\Omega f\, v\, dx
 *   - \int_{\partial\Omega} (\boldsymbol{\beta}\cdot\mathbf{n})^-\, g\, v\, ds.
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
        // Assign attribute 1 to every cell for consistent FaceNormal orientation.
        for (auto it = mesh.getCell(); it; ++it)
          mesh.setAttribute({ mesh.getDimension(), it->getIndex() }, 1);
        return mesh;
      }
  };

  using Manufactured_DGAdvection_Test_16x16 =
    Manufactured_DGAdvection_Test<16>;
  using Manufactured_DGAdvection_Test_32x32 =
    Manufactured_DGAdvection_Test<32>;

  /**
   * @brief Advection with beta=(1,0) and manufactured solution u = x(1-x)y(1-y).
   *
   * For beta=(1,0) the advection PDE reduces to d_x u = f.
   * With u = x(1-x)y(1-y) the source is f = (1-2x)y(1-y)
   * and the inflow boundary condition at x=0 is g=0.
   */
  TEST_P(Manufactured_DGAdvection_Test_32x32, DGAdvection_HorizontalFlow)
  {
    Mesh mesh = this->getMesh();

    // Advection velocity beta = (1, 0)
    Math::Vector<Real> betaVec(2);
    betaVec << 1.0, 0.0;
    VectorFunction beta(betaVec);

    // Manufactured solution and source
    auto f = (1.0 - 2.0 * F::x) * F::y * (1.0 - F::y);
    auto g = Zero();

    P0 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);

    FaceNormal     fn(mesh);
    BoundaryNormal bn(mesh);
    auto fn1    = fn.traceOf(1);
    auto betaN  = Dot(beta, fn1);
    auto betaBN = Dot(beta, bn);

    Problem advection(u, v);
    advection = InterfaceIntegral(betaN * Average(u), Jump(v))
              + InterfaceIntegral(Abs(betaN) / 2.0 * Jump(u), Jump(v))
              + BoundaryIntegral(Max(betaBN, Zero()) * u, v)
              - Integral(f, v)
              - BoundaryIntegral(Min(betaBN, Zero()) * g, v);

    SparseLU(advection).solve();

    auto solution = F::x * (1.0 - F::x) * F::y * (1.0 - F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real error = Integral(diff).compute();
    EXPECT_LT(error, 0.01);
  }

  /**
   * @brief Advection with beta=(1,1)/sqrt(2) and manufactured solution u = sin(pi*x)sin(pi*y).
   *
   * Uses a diagonal flow direction to exercise inflow/outflow on all four
   * boundary segments via the Max/Min flux splitting.
   */
  TEST_P(Manufactured_DGAdvection_Test_32x32, DGAdvection_DiagonalFlow)
  {
    Mesh mesh = this->getMesh();

    const Real pi    = Math::Constants::pi();
    const Real inv_r = 1.0 / std::sqrt(2.0);

    // Advection velocity beta = (1/sqrt(2), 1/sqrt(2))
    Math::Vector<Real> betaVec(2);
    betaVec << inv_r, inv_r;
    VectorFunction beta(betaVec);

    // Manufactured solution: u = sin(pi*x)sin(pi*y)
    auto solution = sin(pi * F::x) * sin(pi * F::y);
    auto f        = inv_r * pi * (cos(pi * F::x) * sin(pi * F::y)
                                + sin(pi * F::x) * cos(pi * F::y));
    auto g        = solution;

    P0 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);

    FaceNormal     fn(mesh);
    BoundaryNormal bn(mesh);
    auto fn1    = fn.traceOf(1);
    auto betaN  = Dot(beta, fn1);
    auto betaBN = Dot(beta, bn);

    Problem advection(u, v);
    advection = InterfaceIntegral(betaN * Average(u), Jump(v))
              + InterfaceIntegral(Abs(betaN) / 2.0 * Jump(u), Jump(v))
              + BoundaryIntegral(Max(betaBN, Zero()) * u, v)
              - Integral(f, v)
              - BoundaryIntegral(Min(betaBN, Zero()) * g, v);

    SparseLU(advection).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real error = Integral(diff).compute();
    EXPECT_LT(error, 0.01);
  }

  // Register parameterised tests for triangles
  INSTANTIATE_TEST_SUITE_P(Triangle, Manufactured_DGAdvection_Test_16x16,
    ::testing::Values(Polytope::Type::Triangle));

  INSTANTIATE_TEST_SUITE_P(Triangle, Manufactured_DGAdvection_Test_32x32,
    ::testing::Values(Polytope::Type::Triangle));

} // namespace Rodin::Tests::Manufactured::DGAdvection
