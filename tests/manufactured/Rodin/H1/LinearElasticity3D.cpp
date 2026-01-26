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
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::H1LinearElasticity3D
{
  template <size_t M>
  class Manufactured_LinearElasticity3D_H1_Test : public ::testing::Test
  {
  protected:
    void SetUp() override
    {
      m_mesh = Mesh().UniformGrid(Polytope::Type::Tetrahedron, { M, M, M });
      m_mesh.scale(1.0 / (M - 1));
      m_mesh.getConnectivity().compute(2, 3);
      m_mesh.getConnectivity().compute(3, 2);
      m_mesh.getConnectivity().compute(2, 1);
      m_mesh.getConnectivity().compute(1, 0);
    }

    const auto& getMesh() const { return m_mesh; }

    static VectorFunction polynomialSolution()
    {
      return VectorFunction{
        F::x * (1 - F::x),
        F::y * (1 - F::y),
        F::z * (1 - F::z)
      };
    }

    static VectorFunction polynomialForcing(Real lambda, Real mu)
    {
      // For u_i = coord(1 - coord):
      // Δu_i = -2, div u = 3 - 2(x+y+z), grad(div u) = (-2,-2,-2)
      // f = -mu Δu - (lambda + mu) grad(div u)
      return VectorFunction{
        2 * mu + 2 * (lambda + mu),
        2 * mu + 2 * (lambda + mu),
        2 * mu + 2 * (lambda + mu)
      };
    }

  private:
    Mesh<Context::Local> m_mesh;
  };

  using Manufactured_LinearElasticity3D_H1_Test_8 =
    Manufactured_LinearElasticity3D_H1_Test<8>;

  TEST_F(Manufactured_LinearElasticity3D_H1_Test_8, LinearElasticity_Polynomial_H1_Order2)
  {
    constexpr auto order = std::integral_constant<size_t, 2>{};
    const Real lambda = 2.0;
    const Real mu = 1.0;

    const auto& mesh = this->getMesh();

    H1 vh(order, mesh, mesh.getSpaceDimension());

    TrialFunction u(vh);
    TestFunction  v(vh);

    const auto u_exact = polynomialSolution();
    const auto f = polynomialForcing(lambda, mu);

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    GridFunction diff(vh);
    diff = Pow(Frobenius(u.getSolution() - u_exact), 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  TEST_F(Manufactured_LinearElasticity3D_H1_Test_8, LinearElasticity_Polynomial_H1_Order3)
  {
    constexpr auto order = std::integral_constant<size_t, 3>{};
    const Real lambda = 1.5;
    const Real mu = 0.5;

    const auto& mesh = this->getMesh();

    H1 vh(order, mesh, mesh.getSpaceDimension());

    TrialFunction u(vh);
    TestFunction  v(vh);

    const auto u_exact = polynomialSolution();
    const auto f = polynomialForcing(lambda, mu);

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    GridFunction diff(vh);
    diff = Pow(Frobenius(u.getSolution() - u_exact), 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }
}
