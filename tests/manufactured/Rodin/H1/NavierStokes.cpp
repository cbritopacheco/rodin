/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <limits>
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Solver/GMRES.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::NavierStokes
{
  template <size_t M>
  class Manufactured_NavierStokes_Test : public ::testing::TestWithParam<Polytope::Type>
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

  using Manufactured_NavierStokes_Test_16x16 =
    Manufactured_NavierStokes_Test<16>;

  TEST_P(Manufactured_NavierStokes_Test_16x16, NavierStokes_Picard_TaylorGreen)
  {
    constexpr Real nu = 1.0;
    // Keep the same Picard guard style as Stokes.cpp/NavierStokes_Picard_TaylorGreen.
    constexpr size_t maxIts = 12;
    constexpr Real tol = 1e-10;
    const auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);
    P0g p0g(mesh);
    P1 sh(mesh);

    VectorFunction u_exact{
      sin(pi * F::x) * cos(pi * F::y),
      -cos(pi * F::x) * sin(pi * F::y)
    };
    auto p_exact = -0.25 * (cos(2 * pi * F::x) + cos(2 * pi * F::y));

    VectorFunction f{
      2 * nu * pi * pi * sin(pi * F::x) * cos(pi * F::y) + pi * sin(2 * pi * F::x),
      -2 * nu * pi * pi * cos(pi * F::x) * sin(pi * F::y) + pi * sin(2 * pi * F::y)
    };

    TrialFunction u(uh);
    TrialFunction p(ph);
    TrialFunction lambda(p0g);
    TestFunction  v(uh);
    TestFunction  q(ph);
    TestFunction  mu(p0g);

    GridFunction w(uh);
    w = VectorFunction{ Zero(), Zero() };
    GridFunction u_prev(uh);
    u_prev = VectorFunction{ Zero(), Zero() };
    GridFunction p_prev(ph);
    p_prev = 0.0;

    Real updateNorm = std::numeric_limits<Real>::infinity();
    size_t picardIts = 0;
    for (size_t k = 0; k < maxIts; ++k)
    {
      Problem ns(u, p, lambda, v, q, mu);
      const auto conv_u = Mult(Jacobian(u), w);
      ns = nu * Integral(Jacobian(u), Jacobian(v))
         + Integral(Dot(conv_u, v))
         - Integral(p, Div(v))
         + Integral(Div(u), q)
         + Integral(lambda, q)
         + Integral(p, mu)
         - Integral(f, v)
         + DirichletBC(u, u_exact);

      SparseLU solver(ns);
      solver.solve();

      GridFunction u_new(uh);
      u_new.setData(u.getSolution().getData());
      GridFunction p_new(ph);
      p_new.setData(p.getSolution().getData());

      GridFunction du(sh);
      du = Pow(Frobenius(u_new - u_prev), 2);
      updateNorm = std::sqrt(Integral(du).compute());

      w.setData(u_new.getData());
      u_prev.setData(u_new.getData());
      p_prev.setData(p_new.getData());

      picardIts = k + 1;
      if (updateNorm < tol)
        break;
    }

    EXPECT_LT(picardIts, maxIts);
    EXPECT_LT(updateNorm, tol);

    GridFunction diff_u(sh);
    diff_u = Pow(Frobenius(u_prev - u_exact), 2);
    const Real error_u = Integral(diff_u).compute();

    const Real vol = mesh.getMeasure(mesh.getDimension());
    GridFunction mean(sh);
    mean = p_prev - p_exact;
    const Real mean_diff = Integral(mean).compute() / vol;

    GridFunction diff_p(sh);
    diff_p = Pow((p_prev - mean_diff) - p_exact, 2);
    const Real error_p = Integral(diff_p).compute();

    // Thresholds tuned to this manufactured setup on both triangle and quadrilateral meshes.
    EXPECT_NEAR(error_u, 0, 1e-8);
    EXPECT_NEAR(error_p, 0, 2e-5);

    GridFunction w_final(uh);
    w_final.setData(u_prev.getData());

    Problem ns_final(u, p, lambda, v, q, mu);
    const auto conv_u_final = Mult(Jacobian(u), w_final);
    ns_final = nu * Integral(Jacobian(u), Jacobian(v))
             + Integral(Dot(conv_u_final, v))
             - Integral(p, Div(v))
             + Integral(Div(u), q)
             + Integral(lambda, q)
             + Integral(p, mu)
             - Integral(f, v)
             + DirichletBC(u, u_exact);

    SparseLU solver_final(ns_final);
    solver_final.solve();

    auto& ls = ns_final.getLinearSystem();
    auto& A = ls.getOperator();
    auto& b = ls.getVector();
    auto& x = ls.getSolution();
    auto r = (A * x - b).eval();
    const Real scale = std::max<Real>(b.norm(), 1);
    EXPECT_NEAR(r.norm() / scale, 0, 1e-12);

    GridFunction u_final(uh);
    u_final.setData(u.getSolution().getData());
    GridFunction p_final(ph);
    p_final.setData(p.getSolution().getData());

    GridFunction du_fp(sh);
    du_fp = Pow(Frobenius(u_final - u_prev), 2);
    const Real fixedPointVelocityDefect = std::sqrt(Integral(du_fp).compute());
    EXPECT_LT(fixedPointVelocityDefect, 1e-8);

    GridFunction mean_fp(sh);
    mean_fp = p_final - p_prev;
    const Real mean_fp_diff = Integral(mean_fp).compute() / vol;

    GridFunction dp_fp(sh);
    dp_fp = Pow((p_final - mean_fp_diff) - p_prev, 2);
    const Real fixedPointPressureDefect = std::sqrt(Integral(dp_fp).compute());
    EXPECT_LT(fixedPointPressureDefect, 1e-6);
  }

  TEST_P(Manufactured_NavierStokes_Test_16x16, NavierStokes_Picard_PolynomialVortex)
  {
    constexpr size_t maxPicardIters = 12;
    constexpr Real picardTol = 1e-11;

    Mesh mesh = this->getMesh();

    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);
    P0g p0g(mesh);
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);

    VectorFunction u_exact{ F::y, -F::x };
    RealFunction p_exact = 0.0;
    VectorFunction f{ -F::x, -F::y };

    TrialFunction u(uh);
    TrialFunction p(ph);
    TrialFunction lambda(p0g);
    TestFunction  v(uh);
    TestFunction  q(ph);
    TestFunction  mu(p0g);

    GridFunction u_picard(uh);
    u_picard = Math::Vector<Real>{{ 0.0, 0.0 }};

    bool converged = false;
    for (size_t k = 0; k < maxPicardIters; ++k)
    {
      Problem ns(u, p, lambda, v, q, mu);
      ns = Integral(Jacobian(u), Jacobian(v))
         + Integral(Dot(Mult(Jacobian(u), u_picard), v))
         - Integral(p, Div(v))
         + Integral(Div(u), q)
         + Integral(lambda, q)
         + Integral(p, mu)
         - Integral(f, v)
         + DirichletBC(u, u_exact);

      SparseLU solver(ns);
      solver.solve();

      GridFunction diff_iter(sh);
      diff_iter = Pow(Frobenius(u.getSolution() - u_picard), 2);
      const Real picardDelta = Integral(diff_iter).compute();

      u_picard = u.getSolution();

      if (picardDelta < picardTol)
      {
        converged = true;
        break;
      }
    }

    EXPECT_TRUE(converged);

    GridFunction diff_u(sh);
    diff_u = Pow(Frobenius(u.getSolution() - u_exact), 2);
    const Real error_u = Integral(diff_u).compute();

    GridFunction diff_p(sh);
    diff_p = Pow(p.getSolution() - p_exact, 2);
    const Real error_p = Integral(diff_p).compute();

    EXPECT_NEAR(error_u, 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(error_p, 0, 1e-8);
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    Manufactured_NavierStokes_Test_16x16,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

}
