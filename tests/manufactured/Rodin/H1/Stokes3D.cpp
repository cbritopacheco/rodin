/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Solver/ForwardDecls.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Solver/DGMRES.h"
#include "Rodin/Solver/IDRS.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::Stokes3D
{
  template <size_t NX, size_t NY, size_t NZ>
  class Manufactured_Stokes3D_Test : public ::testing::TestWithParam<Polytope::Type>
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

  using Manufactured_Stokes3D_Test_12 =
    Manufactured_Stokes3D_Test<12, 12, 12>;

  TEST_P(Manufactured_Stokes3D_Test_12, Stokes3D_P1ExactResidual)
  {
    Mesh mesh = this->getMesh();

    H1  uh(std::integral_constant<size_t, 1>{}, mesh, mesh.getSpaceDimension());
    P1  ph(mesh);
    P0g p0g(mesh);

    VectorFunction u_exact{ F::y, -F::y, Zero() }; // divergence-free affine (continuous)
    RealFunction   p_exact = 0.0;
    VectorFunction f{ Zero(), Zero(), Zero() };

    TrialFunction u(uh);
    TrialFunction p(ph);
    TestFunction  v(uh);
    TestFunction  q(ph);

    TrialFunction lambda(p0g);
    TestFunction  mu(p0g);

    // NOTE: keep the same sign convention as your 2D fixed test:
    // - Integral(Div(u), q) to match - Integral(p, Div(v)) for the usual block structure.
    Problem stokes(u, p, lambda, v, q, mu);
    stokes = Integral(Jacobian(u), Jacobian(v))
           - Integral(p, Div(v))
           - Integral(Div(u), q)
           + Integral(lambda, q)
           + Integral(p, mu)
           - Integral(f, v)
           + DirichletBC(u, u_exact);

    stokes.assemble();

    DGMRES gmres(stokes);
    gmres.setTolerance(1e-10);
    gmres.setRestart(50);
    gmres.setMaxIterations(500);
    gmres.solve();

    GridFunction u_exact_coeffs(uh);
    u_exact_coeffs = u_exact;

    GridFunction p_exact_coeffs(ph);
    p_exact_coeffs = p_exact;

    auto& ls = stokes.getLinearSystem();
    auto& A  = ls.getOperator();
    auto& b  = ls.getVector();
    auto& x  = ls.getSolution();

    const auto uSize = u_exact_coeffs.getData().size();
    const auto pSize = p_exact_coeffs.getData().size();
    const auto lambdaIndex = static_cast<Eigen::Index>(uSize + pSize);

    // Build x_exact = [u_exact, p_exact, lambda_exact], with lambda chosen consistently
    // with the assembled discrete system (1D LS minimizer on the pressure block).
    auto x_exact = x;
    x_exact.head(uSize) = u_exact_coeffs.getData();
    x_exact.segment(uSize, pSize) = p_exact_coeffs.getData();

    // Compute lambda_exact via assembled operator:
    // re0_p = pressure-block residual with lambda=0
    // g     = A_pλ column (pressure rows)
    // lambda* = -(g·re0_p)/(g·g)
    x_exact[lambdaIndex] = 0.0;

    auto re0   = (A * x_exact - b).eval();
    auto re0_p = re0.segment(uSize, pSize);

    decltype(x_exact) eL = x_exact;
    eL.setZero();
    eL[lambdaIndex] = 1.0;

    auto g_full = (A * eL).eval();
    auto g      = g_full.segment(uSize, pSize);

    Real lambda_exact = 0.0;
    const Real gg = g.squaredNorm();
    if (gg > 0)
      lambda_exact = - g.dot(re0_p) / gg;

    x_exact[lambdaIndex] = lambda_exact;

    auto r  = (A * x - b).eval();
    auto re = (A * x_exact - b).eval();

    const Real scale = std::max<Real>(b.norm(), 1);
    EXPECT_NEAR(r.norm() / scale, 0, 1e-8);

    // Block checks (do NOT expect full re-norm to be ~0 with P0g gauge)
    auto re_u  = re.head(uSize);
    auto re_p  = re.segment(uSize, pSize);
    Real re_mu = re[lambdaIndex];

    EXPECT_NEAR(re_u.norm(), 0, 1e-10);  // assembly + Dirichlet consistency
    EXPECT_NEAR(re_mu,       0, 1e-10);  // mean(p)=0 constraint row (p_exact=0)

    // Orthogonality condition at LS minimizer: re_p ⟂ span{g}
    const Real gnorm = std::max<Real>(g.norm(), 1);
    EXPECT_NEAR(g.dot(re_p) / gnorm, 0, 1e-10);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Manufactured_Stokes3D_Test_12,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge
      )
  );
}
