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
#include "Rodin/Solver/DGMRES.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::NavierStokes3D
{
  template <size_t NX, size_t NY, size_t NZ>
  class Manufactured_NavierStokes3D_Test : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      Mesh<Context::Local> getMesh()
      {
        Mesh mesh;
        mesh = mesh.UniformGrid(GetParam(), { NX, NY, NZ });
        mesh.scale(1.0 / (NX - 1));
        mesh.getConnectivity().compute(2, 3);
        mesh.getConnectivity().compute(3, 2);
        mesh.getConnectivity().compute(2, 1);
        mesh.getConnectivity().compute(1, 0);
        return mesh;
      }
  };

  using Manufactured_NavierStokes3D_Test_5 =
    Manufactured_NavierStokes3D_Test<5, 5, 5>;

  TEST_P(Manufactured_NavierStokes3D_Test_5, NavierStokes3D_Picard_TaylorGreen)
  {
    constexpr Real nu = 1.0;
    constexpr size_t maxPicardIters = 8;
    constexpr Real picardTol = 1e-8;
    const auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    P1 ph(mesh);
    P0g p0g(mesh);
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);

    VectorFunction u_exact{
      Sin(pi * F::x) * Cos(pi * F::y) * Cos(pi * F::z),
      -Cos(pi * F::x) * Sin(pi * F::y) * Cos(pi * F::z),
      Zero()
    };
    RealFunction p_exact = 0.0;

    VectorFunction f{
      3 * nu * pi * pi * Sin(pi * F::x) * Cos(pi * F::y) * Cos(pi * F::z)
      + 0.5 * pi * Sin(2 * pi * F::x) * Pow(Cos(pi * F::z), 2),

      -3 * nu * pi * pi * Cos(pi * F::x) * Sin(pi * F::y) * Cos(pi * F::z)
      + 0.5 * pi * Sin(2 * pi * F::y) * Pow(Cos(pi * F::z), 2),

      Zero()
    };

    TrialFunction u(uh);
    TrialFunction p(ph);
    TrialFunction lambda(p0g);
    TestFunction  v(uh);
    TestFunction  q(ph);
    TestFunction  mu(p0g);

    GridFunction u_picard(uh);
    u_picard = Math::Vector<Real>{{ 0.0, 0.0, 0.0 }};

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

      DGMRES gmres(ns);
      gmres.setTolerance(1e-10);
      gmres.setRestart(100);
      gmres.setMaxIterations(2000);
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

    EXPECT_NEAR(error_u, 0, 1e-4);
    EXPECT_NEAR(error_p, 0, 2e-2);
  }

  TEST_P(Manufactured_NavierStokes3D_Test_5, NavierStokes3D_Picard_PolynomialVortex)
  {
    constexpr size_t maxPicardIters = 8;
    constexpr Real picardTol = 1e-9;

    Mesh mesh = this->getMesh();

    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    P1 ph(mesh);
    P0g p0g(mesh);
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);

    VectorFunction u_exact{ F::y, -F::x, Zero() };
    RealFunction p_exact = 0.0;
    VectorFunction f{ -F::x, -F::y, Zero() };

    TrialFunction u(uh);
    TrialFunction p(ph);
    TrialFunction lambda(p0g);
    TestFunction  v(uh);
    TestFunction  q(ph);
    TestFunction  mu(p0g);

    GridFunction u_picard(uh);
    u_picard = Math::Vector<Real>{{ 0.0, 0.0, 0.0 }};

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

      DGMRES gmres(ns);
      gmres.setTolerance(1e-10);
      gmres.setRestart(100);
      gmres.setMaxIterations(2000);
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

    EXPECT_NEAR(error_u, 0, 1e-8);
    EXPECT_NEAR(error_p, 0, 1e-3);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Manufactured_NavierStokes3D_Test_5,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron
    )
  );
}
