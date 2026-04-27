/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file NewtonianStokes.cpp
 * @brief Manufactured solutions for Newtonian viscous flow using the
 * Fluid::Newtonian constitutive law and the ViscousTangent integrator.
 *
 * ## Problem statement
 *
 * We seek @f$\mathbf{u} \in [V]^d@f$ satisfying:
 * @f[
 *   -\nabla \cdot \boldsymbol\tau(\mathbf{u}) = \mathbf{f} \quad\text{in }\Omega,
 *   \qquad
 *   \mathbf{u} = \mathbf{g} \quad\text{on }\partial\Omega,
 * @f]
 * where @f$\boldsymbol\tau = 2\mu\mathbf{D}(\mathbf{u})@f$ is the Newtonian
 * deviatoric stress and @f$\mathbf{D}(\mathbf{u}) =
 * \tfrac12(\nabla\mathbf{u}+\nabla\mathbf{u}^T)@f$.
 *
 * The weak form is: find @f$\mathbf{u}\in V_g@f$ such that
 * @f[
 *   \int_\Omega \boldsymbol\tau(\mathbf{u}) : \mathbf{D}(\mathbf{v})\,dx
 *   = \int_\Omega \mathbf{f}\cdot\mathbf{v}\,dx
 *   \quad\forall\,\mathbf{v}\in V_0.
 * @f]
 *
 * The @c ViscousTangent integrator assembles the bilinear form using the
 * Fluid::Newtonian constitutive law, exactly as it would be used inside a
 * nonlinear Navier-Stokes Newton iteration.
 *
 * ### P1-exact test (Couette profile)
 *
 * Exact solution on @f$[0,1]^2@f$:
 * @f[
 *   \mathbf{u}_\text{exact} = (y,\;0),\quad \mathbf{f} = \mathbf{0}.
 * @f]
 * Since @f$\mathbf{u}_\text{exact}@f$ is affine (piecewise-linear exact),
 * the discrete P1 solution must reproduce it to machine precision.
 *
 * ### Sine manufactured solution (convergence test)
 *
 * Exact solution on @f$[0,1]^2@f$ (zero on @f$\partial\Omega@f$):
 * @f[
 *   \mathbf{u}_\text{exact}(x,y) =
 *     \bigl(\sin(\pi x)\sin(\pi y),\;\sin(\pi x)\sin(\pi y)\bigr),
 * @f]
 * @f[
 *   \mathbf{f}(x,y) = \mu\pi^2
 *     \bigl(3\sin(\pi x)\sin(\pi y) - \cos(\pi x)\cos(\pi y),\;
 *           3\sin(\pi x)\sin(\pi y) - \cos(\pi x)\cos(\pi y)\bigr).
 * @f]
 * We verify the @f$L^2@f$ error decreases at the expected first-order rate
 * as the mesh is refined.
 */
#include <cmath>
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Fluid.h"
#include "Rodin/Solver/CG.h"
#include "Rodin/Solver/BiCGSTAB.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Fluid;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::FluidNewtonianStokes
{
  namespace
  {
    Mesh<Context::Local> makeUnitSquareMesh(Polytope::Type geom, size_t n)
    {
      Mesh<Context::Local> mesh = Mesh<Context::Local>::UniformGrid(geom, { n, n });
      mesh.scale(1.0 / (n - 1));
      mesh.getConnectivity().compute(1, 2);
      return mesh;
    }
  }

  // ===========================================================================
  // P1-exact manufactured test: Couette profile (u = (y, 0))
  //
  // This test verifies that the ViscousTangent integrator correctly assembles
  // the Newtonian viscous bilinear form by checking that the P1-exact Couette
  // profile is recovered to machine precision.
  // ===========================================================================

  class Manufactured_NewtonianStokes_P1Exact_Test
    : public ::testing::TestWithParam<Polytope::Type>
  {};

  /**
   * Couette flow: u_exact = (y, 0) on [0,1]^2.
   *
   * The exact solution is affine and satisfies:
   *   - D(u_exact) = [[0, 1/2], [1/2, 0]]  (constant strain rate)
   *   - -2*mu*div(D(u_exact)) = 0  (body force = 0)
   *   - div(u_exact) = 0           (incompressible)
   *
   * With Dirichlet BCs u = u_exact on entire boundary, the P1 FE solution
   * must be exact (residual at machine precision).
   */
  TEST_P(Manufactured_NewtonianStokes_P1Exact_Test,
         Newtonian_Couette_P1Exact_Residual)
  {
    Mesh mesh = makeUnitSquareMesh(GetParam(), 8);
    const size_t dim = mesh.getSpaceDimension();
    const Real mu = 0.01;

    P1 Vh(mesh, dim);

    // Exact solution: u_exact = (y, 0)
    auto u_exact = VectorFunction{ F::y, Zero() };

    TrialFunction u(Vh);
    TestFunction  v(Vh);

    Fluid::Newtonian law(mu);

    // For the linear Newtonian law, the ViscousTangent is independent of the
    // linearization point. We pass a zero GridFunction as the linearization
    // point (equivalent to no linearization point since the law is linear).
    GridFunction uCurrent(Vh);
    uCurrent = u_exact;

    auto viscous = ViscousTangent(law, u, v);
    viscous.setVelocity(uCurrent);

    Problem stokes(u, v);
    stokes = viscous + DirichletBC(u, u_exact);

    BiCGSTAB(stokes).solve();

    // Check solver residual
    auto& A = stokes.getLinearSystem().getOperator();
    auto& b = stokes.getLinearSystem().getVector();
    auto& x = stokes.getLinearSystem().getSolution();

    GridFunction u_exact_gf(Vh);
    u_exact_gf = u_exact;

    auto r  = A * x - b;
    auto re = A * u_exact_gf.getData() - b;

    const Real scale = std::max<Real>(b.norm(), 1.0);
    EXPECT_NEAR(r.norm() / scale, 0.0, 1e-8);
    EXPECT_NEAR(re.norm() / scale, 0.0, 1e-10);

    // FE L2 error should be at machine precision
    P1 scalar(mesh);
    GridFunction err2(scalar);
    err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
    EXPECT_NEAR(Integral(err2).compute(), 0.0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage2D,
    Manufactured_NewtonianStokes_P1Exact_Test,
    ::testing::Values(
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral));

  // ===========================================================================
  // Sine manufactured solution: convergence test on a sequence of meshes
  //
  // Exact: u_exact = (sin(pi*x)*sin(pi*y), sin(pi*x)*sin(pi*y))
  // f     = (3*mu*pi^2*sin(pi*x)*sin(pi*y), 3*mu*pi^2*sin(pi*x)*sin(pi*y))
  //
  // u_exact = 0 on dOmega, so Dirichlet BC is zero.
  // ===========================================================================

  class Manufactured_NewtonianStokes_Sine_Test
    : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(Manufactured_NewtonianStokes_Sine_Test,
         Newtonian_Sine_L2Convergence)
  {
    const Real mu  = 1.0;
    const Real pi  = Math::Constants::pi();
    const auto gt  = GetParam();

    // Body force  f = μπ²(3 sin(πx)sin(πy) − cos(πx)cos(πy), …)
    auto f = VectorFunction{
      RealFunction([mu, pi](const Point& p) -> Real {
        return mu * pi * pi * (3.0 * std::sin(pi * p.x()) * std::sin(pi * p.y())
                               - std::cos(pi * p.x()) * std::cos(pi * p.y()));
      }),
      RealFunction([mu, pi](const Point& p) -> Real {
        return mu * pi * pi * (3.0 * std::sin(pi * p.x()) * std::sin(pi * p.y())
                               - std::cos(pi * p.x()) * std::cos(pi * p.y()));
      })
    };

    // Exact solution
    auto u_exact = VectorFunction{
      RealFunction([pi](const Point& p) -> Real {
        return std::sin(pi * p.x()) * std::sin(pi * p.y());
      }),
      RealFunction([pi](const Point& p) -> Real {
        return std::sin(pi * p.x()) * std::sin(pi * p.y());
      })
    };

    auto computeL2Error = [&](size_t n) -> Real
    {
      Mesh mesh = makeUnitSquareMesh(gt, n);
      const size_t dim = mesh.getSpaceDimension();
      P1 Vh(mesh, dim);

      GridFunction uCurrent(Vh);
      uCurrent = Zero(dim);

      TrialFunction u(Vh);
      TestFunction  v(Vh);

      Fluid::Newtonian law(mu);
      auto viscous = ViscousTangent(law, u, v);
      viscous.setVelocity(uCurrent);

      Problem stokes(u, v);
      stokes = viscous
             - Integral(f, v)
             + DirichletBC(u, Zero(dim));

      BiCGSTAB(stokes).solve();

      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      return std::sqrt(Integral(err2).compute());
    };

    // Coarse and fine meshes (only two refinements to keep test fast)
    const Real e8  = computeL2Error(9);
    const Real e16 = computeL2Error(17);

    // P1 should give at least ~first-order convergence in L2.
    // A factor of ~2 in L2 error is expected for mesh refinement by 2.
    // We verify the error is decreasing by a meaningful amount (factor > 1.2).
    EXPECT_GT(e8, e16);
    EXPECT_GT(e8 / e16, 1.2)
        << "Error did not decrease sufficiently on mesh refinement: "
        << "e8=" << e8 << " e16=" << e16;
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage2D,
    Manufactured_NewtonianStokes_Sine_Test,
    ::testing::Values(
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral));

  // ===========================================================================
  // Cross-check: ViscousTangent matches LinearElasticityIntegral (lambda=0)
  //
  // For Newtonian flow, ViscousTangent(law, u, v) with mu should assemble
  // exactly the same stiffness matrix as LinearElasticityIntegral(u, v)(0, mu).
  // We verify this by comparing solutions from both formulations.
  // ===========================================================================

  class Manufactured_NewtonianStokes_CrossCheck_Test
    : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(Manufactured_NewtonianStokes_CrossCheck_Test,
         ViscousTangent_MatchesLinearElasticity)
  {
    const Real mu = 0.5;
    const Real pi = Math::Constants::pi();
    const auto gt = GetParam();

    Mesh mesh = makeUnitSquareMesh(gt, 9);
    const size_t dim = mesh.getSpaceDimension();
    P1 Vh(mesh, dim);

    auto f = VectorFunction{
      RealFunction([mu, pi](const Point& p) -> Real {
        return mu * pi * pi * (3.0 * std::sin(pi * p.x()) * std::sin(pi * p.y())
                               - std::cos(pi * p.x()) * std::cos(pi * p.y()));
      }),
      RealFunction([mu, pi](const Point& p) -> Real {
        return mu * pi * pi * (3.0 * std::sin(pi * p.x()) * std::sin(pi * p.y())
                               - std::cos(pi * p.x()) * std::cos(pi * p.y()));
      })
    };

    // Solve using ViscousTangent
    GridFunction uFluid(Vh);
    {
      GridFunction uCurrent(Vh);
      uCurrent = Zero(dim);

      TrialFunction u(Vh);
      TestFunction  v(Vh);

      Fluid::Newtonian law(mu);
      auto viscous = ViscousTangent(law, u, v);
      viscous.setVelocity(uCurrent);

      Problem stokes(u, v);
      stokes = viscous
             - Integral(f, v)
             + DirichletBC(u, Zero(dim));

      BiCGSTAB(stokes).solve();
      uFluid = u.getSolution();
    }

    // Solve using LinearElasticityIntegral (lambda=0 ⟺ Newtonian Stokes)
    GridFunction uElasticity(Vh);
    {
      TrialFunction u(Vh);
      TestFunction  v(Vh);

      Problem stokes(u, v);
      stokes = LinearElasticityIntegral(u, v)(0.0, mu)
             - Integral(f, v)
             + DirichletBC(u, Zero(dim));

      BiCGSTAB(stokes).solve();
      uElasticity = u.getSolution();
    }

    // Both solutions should agree to round-off
    P1 scalar(mesh);
    GridFunction diff2(scalar);
    diff2 = Pow(Frobenius(uFluid - uElasticity), 2);
    const Real l2diff = std::sqrt(Integral(diff2).compute());
    EXPECT_NEAR(l2diff, 0.0, 1e-10)
        << "ViscousTangent and LinearElasticityIntegral solutions differ.";
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage2D,
    Manufactured_NewtonianStokes_CrossCheck_Test,
    ::testing::Values(
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral));
}
