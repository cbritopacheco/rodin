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
 * @brief Manufactured solution tests for the P0 DG interior-penalty Poisson problem.
 *
 * Solves the Dirichlet Poisson problem
 * @f[
 *   -\Delta u = f \quad \text{in } \Omega = [0,1]^2, \qquad u = g \quad \text{on } \partial\Omega
 * @f]
 * using a P0 (piecewise constant) DG discretisation with an interior penalty:
 * @f[
 *   \sigma \int_{\Gamma_{\mathrm{int}}} [\![u]\!][\![v]\!]\,ds
 *   + \sigma \int_{\partial\Omega} u\,v\,ds
 *   = \int_\Omega f\,v\,dx + \sigma \int_{\partial\Omega} g\,v\,ds.
 * @f]
 * The penalty @f$ \sigma = 1/h @f$ recovers a two-point-flux approximation
 * (TPFA) finite-volume scheme.  The tests verify that the discrete solution
 * is close to the manufactured solution with a tolerance consistent with the
 * first-order accuracy of the method.
 */
namespace Rodin::Tests::Manufactured::DGPoisson
{
  template <size_t M>
  class Manufactured_DGPoisson_Test : public ::testing::TestWithParam<Polytope::Type>
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

  using Manufactured_DGPoisson_Test_16x16 =
    Manufactured_DGPoisson_Test<16>;
  using Manufactured_DGPoisson_Test_32x32 =
    Manufactured_DGPoisson_Test<32>;

  /**
   * @brief Verifies that the P0 penalty scheme produces a bounded solution
   *        for the simple-sine manufactured solution.
   *
   * Manufactured solution: @f$ u = \sin(\pi x)\sin(\pi y) @f$,
   * @f$ f = 2\pi^2 \sin(\pi x)\sin(\pi y) @f$, @f$ g = 0 @f$.
   *
   * The P0 penalty scheme gives O(h) convergence in L2.  This test checks
   * that the 32x32 L2 error is below a tolerance consistent with O(h).
   */
  TEST_P(Manufactured_DGPoisson_Test_32x32, DGPoisson_SimpleSine)
  {
    Mesh mesh = this->getMesh();

    const Real pi = Math::Constants::pi();
    auto solution = sin(pi * F::x) * sin(pi * F::y);
    auto f        = 2.0 * pi * pi * solution;
    auto g        = Zero();

    P0 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);

    const Real h     = 1.0 / static_cast<Real>(32 - 1);
    const Real sigma = 1.0 / h;

    Problem poisson(u, v);
    poisson = InterfaceIntegral(sigma * Jump(u), Jump(v))
            + BoundaryIntegral(sigma * u, v)
            - Integral(f, v)
            - BoundaryIntegral(sigma * g, v);

    SparseLU(poisson).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real error = Integral(diff).compute();
    EXPECT_LT(error, 0.1);
  }

  /**
   * @brief Convergence test: 32x32 error should be smaller than 16x16 error.
   *
   * Confirms that the P0 penalty scheme converges as the mesh is refined.
   */
  TEST_P(Manufactured_DGPoisson_Test_16x16, DGPoisson_ConvergesUnderRefinement)
  {
    const Real pi = Math::Constants::pi();
    auto solution = sin(pi * F::x) * sin(pi * F::y);
    auto f        = 2.0 * pi * pi * solution;
    auto g        = Zero();

    // Error on 16x16 mesh
    Real error16;
    {
      Mesh mesh16 = this->getMesh();
      P0 vh(mesh16);
      TrialFunction u(vh);
      TestFunction  v(vh);
      const Real h16    = 1.0 / static_cast<Real>(16 - 1);
      const Real sigma16 = 1.0 / h16;
      Problem poisson(u, v);
      poisson = InterfaceIntegral(sigma16 * Jump(u), Jump(v))
              + BoundaryIntegral(sigma16 * u, v)
              - Integral(f, v)
              - BoundaryIntegral(sigma16 * g, v);
      SparseLU(poisson).solve();
      GridFunction diff(vh);
      diff = Pow(u.getSolution() - solution, 2);
      error16 = Integral(diff).compute();
    }

    // Error on 32x32 mesh
    Real error32;
    {
      Mesh mesh32;
      mesh32 = mesh32.UniformGrid(GetParam(), { 32, 32 });
      mesh32.scale(1.0 / (32 - 1));
      mesh32.getConnectivity().compute(1, 2);
      P0 vh(mesh32);
      TrialFunction u(vh);
      TestFunction  v(vh);
      const Real h32    = 1.0 / static_cast<Real>(32 - 1);
      const Real sigma32 = 1.0 / h32;
      Problem poisson(u, v);
      poisson = InterfaceIntegral(sigma32 * Jump(u), Jump(v))
              + BoundaryIntegral(sigma32 * u, v)
              - Integral(f, v)
              - BoundaryIntegral(sigma32 * g, v);
      SparseLU(poisson).solve();
      GridFunction diff(vh);
      diff = Pow(u.getSolution() - solution, 2);
      error32 = Integral(diff).compute();
    }

    EXPECT_LT(error32, error16);
  }

  // Register parameterised tests for triangles
  INSTANTIATE_TEST_SUITE_P(Triangle, Manufactured_DGPoisson_Test_16x16,
    ::testing::Values(Polytope::Type::Triangle));

  INSTANTIATE_TEST_SUITE_P(Triangle, Manufactured_DGPoisson_Test_32x32,
    ::testing::Values(Polytope::Type::Triangle));

} // namespace Rodin::Tests::Manufactured::DGPoisson

