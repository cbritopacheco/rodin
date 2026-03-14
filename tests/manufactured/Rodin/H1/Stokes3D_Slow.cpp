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
#include "Rodin/Solver/SparseLU.h"

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

  TEST_P(Manufactured_Stokes3D_Test_12, Stokes3D_Trigonometric)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    H1 uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
    P1 ph(mesh);

    P0g p0g(mesh);

    VectorFunction u_exact{
      Sin(pi * F::x) * Cos(pi * F::y) * Cos(pi * F::z),
      -Cos(pi * F::x) * Sin(pi * F::y) * Cos(pi * F::z),
      Zero()
    };
    auto p_exact = Cos(2 * pi * F::x) * Cos(2 * pi * F::y) * Cos(2 * pi * F::z);

    VectorFunction f{
      3 * pi * pi * Sin(pi * F::x) * Cos(pi * F::y) * Cos(pi * F::z)
      - 2 * pi * Sin(2 * pi * F::x) * Cos(2 * pi * F::y) * Cos(2 * pi * F::z),

      -3 * pi * pi * Cos(pi * F::x) * Sin(pi * F::y) * Cos(pi * F::z)
      - 2 * pi * Cos(2 * pi * F::x) * Sin(2 * pi * F::y) * Cos(2 * pi * F::z),

      -2 * pi * Cos(2 * pi * F::x) * Cos(2 * pi * F::y) * Sin(2 * pi * F::z)
    };


    TrialFunction u(uh);
    TrialFunction p(ph);
    TestFunction  v(uh);
    TestFunction  q(ph);

    TrialFunction lambda(p0g);
    TestFunction  mu(p0g);

    Problem stokes(u, p, v, q, lambda, mu);
    stokes = Integral(Jacobian(u), Jacobian(v))
           - Integral(p, Div(v))
           + Integral(Div(u), q)
           + Integral(lambda, q)
           + Integral(p, mu)
           - Integral(f, v)
           + DirichletBC(u, u_exact);

    stokes.assemble();

    SparseLU solver(stokes);
    solver.solve();

    // L2 errors
    H1 sh(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction diff_u(sh);
    diff_u = Pow(Frobenius(u.getSolution() - u_exact), 2);
    Real error_u = Integral(diff_u).compute();

    GridFunction diff_p(sh);
    diff_p = Pow(p.getSolution() - p_exact, 2);
    Real error_p = Integral(diff_p).compute();

    EXPECT_NEAR(error_u, 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(error_p, 0, 2e-3);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Manufactured_Stokes3D_Test_12,
    ::testing::Values(
      // Polytope::Type::Tetrahedron,
      // Polytope::Type::Hexahedron,
      Polytope::Type::Wedge
      )
  );
}
