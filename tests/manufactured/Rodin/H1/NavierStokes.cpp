/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Solver/GMRES.h"

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
    constexpr size_t maxPicardIters = 12;
    constexpr Real picardTol = 1e-10;
    const auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);
    P0g p0g(mesh);
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);

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

    GridFunction u_picard(uh);
    u_picard = Math::Vector<Real>{{ 0.0, 0.0 }};

    bool converged = false;
    for (size_t k = 0; k < maxPicardIters; ++k)
    {
      Problem ns(u, p, lambda, v, q, mu);
      ns = nu * Integral(Jacobian(u), Jacobian(v))
         + Integral(Dot(Mult(Jacobian(u), u_picard), v))
         + 0.5 * Integral(Div(u_picard) * Dot(u, v))
         - Integral(p, Div(v))
         + Integral(Div(u), q)
         + Integral(lambda, q)
         + Integral(p, mu)
         - Integral(f, v)
         + DirichletBC(u, u_exact);

      GMRES gmres(ns);
      gmres.setTolerance(1e-12);
      gmres.setMaxIterations(1000);
      gmres.solve();

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

    const Real vol = mesh.getMeasure(mesh.getDimension());
    GridFunction mean(sh);
    mean = p.getSolution() - p_exact;
    const Real mean_diff = Integral(mean).compute() / vol;

    GridFunction diff_p(sh);
    diff_p = Pow((p.getSolution() - mean_diff) - p_exact, 2);
    const Real error_p = Integral(diff_p).compute();

    EXPECT_NEAR(error_u, 0, 5e-7);
    EXPECT_NEAR(error_p, 0, 5e-4);
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

      GMRES gmres(ns);
      gmres.setTolerance(1e-12);
      gmres.setMaxIterations(1000);
      gmres.solve();

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
    EXPECT_NEAR(error_p, 0, 1e-4);
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    Manufactured_NavierStokes_Test_16x16,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

}
