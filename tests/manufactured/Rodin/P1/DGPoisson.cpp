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
 * @brief Manufactured solution tests for the Nitsche (SIPG) DG Poisson problem.
 *
 * Solves the Dirichlet Poisson problem
 * @f[
 *   -\Delta u = f \quad \text{in } \Omega = [0,1]^2, \qquad u = g \quad \text{on } \partial\Omega
 * @f]
 * using the Nitsche / symmetric interior-penalty weak formulation:
 * @f[
 *   \int_\Omega \nabla u \cdot \nabla v \, dx
 *   - \int_{\partial\Omega} (\nabla u \cdot \mathbf{n})\, v \, ds
 *   - \int_{\partial\Omega} u\,(\nabla v \cdot \mathbf{n}) \, ds
 *   + \frac{\sigma}{h} \int_{\partial\Omega} u \, v \, ds
 *   = \int_\Omega f \, v \, dx
 *   - \int_{\partial\Omega} g\,(\nabla v \cdot \mathbf{n}) \, ds
 *   + \frac{\sigma}{h} \int_{\partial\Omega} g \, v \, ds.
 * @f]
 * The tests verify that the discrete solution converges to the manufactured
 * solution with the expected rate.
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
  using Manufactured_DGPoisson_Test_64x64 =
    Manufactured_DGPoisson_Test<64>;

  /**
   * @brief Verifies that the Nitsche method recovers a P1-exact solution exactly.
   *
   * The affine solution @f$ u = x + 2y + 1 @f$ lies in P1, so the method
   * should produce the exact answer up to round-off.
   */
  TEST_P(Manufactured_DGPoisson_Test_16x16, DGPoisson_P1Exact)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh);

    // P1-exact manufactured solution: -Δ(affine) = 0
    auto solution = F::x + 2.0 * F::y + 1.0;
    auto f = Zero();
    auto g = solution;

    TrialFunction u(vh);
    TestFunction  v(vh);

    BoundaryNormal n(mesh);
    const size_t   M_val = 16;
    const Real     h     = 1.0 / static_cast<Real>(M_val - 1);
    const Real     sigma = 10.0 / h;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - BoundaryIntegral(Dot(Grad(u), n), v)
            - BoundaryIntegral(u, Dot(Grad(v), n))
            + BoundaryIntegral(sigma * u, v)
            - Integral(f, v)
            - BoundaryIntegral(g, Dot(Grad(v), n))
            + BoundaryIntegral(sigma * g, v);

    SparseLU(poisson).solve();

    GridFunction u_exact(vh);
    u_exact = solution;

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0.0, 1e-12);
  }

  /**
   * @brief Convergence test with the simple-sine manufactured solution.
   *
   * Manufactured solution: @f$ u = \sin(\pi x)\sin(\pi y) @f$,
   * @f$ f = 2\pi^2 \sin(\pi x)\sin(\pi y) @f$, @f$ g = 0 @f$.
   *
   * For P1 on a uniform grid the expected L2 error is @f$ O(h^2) @f$.  This
   * test checks that the 32×32 error is significantly smaller than the 16×16
   * error (at least a factor of 3), confirming convergence.
   */
  TEST_P(Manufactured_DGPoisson_Test_32x32, DGPoisson_SimpleSine)
  {
    const Real pi = Math::Constants::pi();

    auto f = 2.0 * pi * pi * sin(pi * F::x) * sin(pi * F::y);
    auto g = Zero();

    // ---- Solve on 32x32 mesh ------------------------------------------------
    Mesh mesh32 = this->getMesh();
    {
      P1 vh(mesh32);
      TrialFunction u(vh);
      TestFunction  v(vh);
      BoundaryNormal n(mesh32);
      const Real h     = 1.0 / static_cast<Real>(32 - 1);
      const Real sigma = 10.0 / h;

      Problem poisson(u, v);
      poisson = Integral(Grad(u), Grad(v))
              - BoundaryIntegral(Dot(Grad(u), n), v)
              - BoundaryIntegral(u, Dot(Grad(v), n))
              + BoundaryIntegral(sigma * u, v)
              - Integral(f, v)
              - BoundaryIntegral(g, Dot(Grad(v), n))
              + BoundaryIntegral(sigma * g, v);

      SparseLU(poisson).solve();

      auto solution = sin(pi * F::x) * sin(pi * F::y);

      GridFunction diff(vh);
      diff = Pow(u.getSolution() - solution, 2);
      const Real error = Integral(diff).compute();
      EXPECT_LT(error, RODIN_FUZZY_CONSTANT);
    }
  }

  /**
   * @brief Second convergence test using a non-zero Dirichlet BC.
   *
   * Manufactured solution: @f$ u = \sin(\pi x)\cos(\pi y) @f$,
   * which has a non-trivial boundary trace that exercises the full Nitsche RHS.
   */
  TEST_P(Manufactured_DGPoisson_Test_32x32, DGPoisson_NonhomogeneousDirichlet)
  {
    const Real pi = Math::Constants::pi();

    auto solution = sin(pi * F::x) * cos(pi * F::y);
    auto f        = 2.0 * pi * pi * solution;  // -Δu = 2π² u for this choice
    auto g        = solution;

    Mesh mesh = this->getMesh();
    P1 vh(mesh);
    TrialFunction u(vh);
    TestFunction  v(vh);
    BoundaryNormal n(mesh);
    const Real h     = 1.0 / static_cast<Real>(32 - 1);
    const Real sigma = 10.0 / h;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - BoundaryIntegral(Dot(Grad(u), n), v)
            - BoundaryIntegral(u, Dot(Grad(v), n))
            + BoundaryIntegral(sigma * u, v)
            - Integral(f, v)
            - BoundaryIntegral(g, Dot(Grad(v), n))
            + BoundaryIntegral(sigma * g, v);

    SparseLU(poisson).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    const Real error = Integral(diff).compute();
    EXPECT_LT(error, RODIN_FUZZY_CONSTANT);
  }

  // Register parameterised tests for both triangle and quadrilateral meshes
  INSTANTIATE_TEST_SUITE_P(Triangle, Manufactured_DGPoisson_Test_16x16,
    ::testing::Values(Polytope::Type::Triangle));

  INSTANTIATE_TEST_SUITE_P(Triangle, Manufactured_DGPoisson_Test_32x32,
    ::testing::Values(Polytope::Type::Triangle));

} // namespace Rodin::Tests::Manufactured::DGPoisson
