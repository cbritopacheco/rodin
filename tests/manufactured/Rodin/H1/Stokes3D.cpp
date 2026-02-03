/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Solver/GMRES.h"

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
      return mesh;
    }
  };

  using Manufactured_Stokes3D_Test_12 =
    Manufactured_Stokes3D_Test<12, 12, 12>;

  // 3D manufactured divergence-free velocity:
  // u = ( sin(pi x) cos(pi y) cos(pi z),
  //       -cos(pi x) sin(pi y) cos(pi z),
  //       0 )
  // ∇·u = 0
  // p = sin(pi x) sin(pi y) sin(pi z)
  // f = -Δu + ∇p
  TEST_P(Manufactured_Stokes3D_Test_12, Stokes3D_SimpleSine)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);

    VectorFunction u_exact{
      sin(pi * F::x) * cos(pi * F::y) * cos(pi * F::z),
      -cos(pi * F::x) * sin(pi * F::y) * cos(pi * F::z),
      Zero()
    };

    auto p_exact = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);

    // f = -Δu + ∇p
    // Δ u1 = -3 pi^2 sin(pi x) cos(pi y) cos(pi z)
    // Δ u2 = 3 pi^2 cos(pi x) sin(pi y) cos(pi z)
    VectorFunction f{
      3 * pi * pi * sin(pi * F::x) * cos(pi * F::y) * cos(pi * F::z)
      + pi * cos(pi * F::x) * sin(pi * F::y) * sin(pi * F::z),
      -3 * pi * pi * cos(pi * F::x) * sin(pi * F::y) * cos(pi * F::z)
      + pi * sin(pi * F::x) * cos(pi * F::y) * sin(pi * F::z),
      pi * cos(pi * F::x) * sin(pi * F::y) * cos(pi * F::z)
    };

    TrialFunction u(uh);
    TrialFunction p(ph);
    TestFunction  v(uh);
    TestFunction  q(ph);

    Problem stokes(u, p, v, q);
    stokes = Integral(Jacobian(u), Jacobian(v))
           - Integral(p, Div(v))
           + Integral(Div(u), q)
           - Integral(f, v)
           + DirichletBC(u, u_exact);

    GMRES gmres(stokes);
    gmres.setTolerance(1e-10);
    gmres.setMaxIterations(500);
    gmres.solve();

    // L2 errors
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction diff_u(sh);
    diff_u = Pow(Frobenius(u.getSolution() - u_exact), 2);
    Real error_u = Integral(diff_u).compute();

    // mean-corrected pressure error
    Real vol = mesh.getMeasure(mesh.getDimension());
    GridFunction mean(sh);
    mean = p.getSolution() - p_exact;
    Real mean_diff = Integral(mean).compute() / vol;

    GridFunction diff_p(sh);
    diff_p = Pow((p.getSolution() - mean_diff) - p_exact, 2);
    Real error_p = Integral(diff_p).compute();

    EXPECT_NEAR(error_u, 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(0.5 * error_p, 0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Manufactured_Stokes3D_Test_12,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge)
  );
}
