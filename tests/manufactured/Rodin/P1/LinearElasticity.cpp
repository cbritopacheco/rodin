/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"
#include "Rodin/Test/Random.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Manufactured::LinearElasticity
{
  template <size_t M>
  class Manufactured_LinearElasticity_Test : public ::testing::TestWithParam<Polytope::Type>
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

  using Manufactured_LinearElasticity_Test_16x16 =
    Rodin::Tests::Manufactured::LinearElasticity::Manufactured_LinearElasticity_Test<16>;
  using Manufactured_LinearElasticity_Test_32x32 =
    Rodin::Tests::Manufactured::LinearElasticity::Manufactured_LinearElasticity_Test<32>;
  using Manufactured_LinearElasticity_Test_64x64 =
    Rodin::Tests::Manufactured::LinearElasticity::Manufactured_LinearElasticity_Test<64>;


  TEST_P(Manufactured_LinearElasticity_Test_16x16, AffineExact)
  {
    Mesh mesh = this->getMesh();
    const Real lambda = 1.0, mu = 1.0;
    const size_t dim = mesh.getSpaceDimension();
    P1 vh(mesh, dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    auto u_exact = VectorFunction{ F::x, F::y };

    Problem elasticity(u, v);
    elasticity = Integral(lambda * Div(u), Div(v))
               + Integral(
                   mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()))
               + DirichletBC(u, u_exact);
    CG(elasticity).solve();

    mesh.save("LinearElasticity.mesh");
    u.getSolution().save("LinearElasticity.gf");

    P1 scalar(mesh);
    GridFunction err2(scalar);
    err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
    const Real L2error = Integral(err2).compute();
    EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    Manufactured_LinearElasticity_Test_16x16,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}
