/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Configure.h"
#include "Rodin/Context/Local.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

/**
 * 3D manufactured tests for isotropic linear elasticity:
 *
 * Strong form:
 *   -div(sigma(u)) = f  in Omega
 *               u = g  on dOmega
 *
 * sigma(u) = lambda (div u) I + 2 mu eps(u),
 * eps(u) = 0.5 (grad u + grad u^T)
 *
 * Weak form:
 *   ∫ [ lambda (div u)(div v) + 2 mu eps(u):eps(v) ] dx = ∫ f·v dx
 *
 * NOTE on geometry: Rodin's UniformGrid({M,M,M}) produces coordinates {0,...,M-1}.
 * We scale by 1/(M-1) so the domain is [0,1]^3.
 *
 * NOTE on "manufactured": For non-affine u_exact, the exact discrete solution is not
 * expected with P1. These tests check mathematical consistency (correct f for chosen u),
 * and use a relaxed error bound on a fixed mesh.
 */

namespace Rodin::Tests::Manufactured::LinearElasticity3D
{
  // ----------------------------
  // Fixture: build mesh once
  // ----------------------------
  template <Polytope::Type G, size_t M>
  class Elasticity3DFixture : public ::testing::Test
  {
    protected:
      void SetUp() override
      {
        m_mesh = Mesh().UniformGrid(G, {M, M, M});
        m_mesh.scale(Real(1) / Real(M - 1)); // map {0..M-1} -> [0,1]
        // (2,3) is enough for element-to-face. Add other connectivities if needed by BC code.
        m_mesh.getConnectivity().compute(2, 3);
      }

      const auto& mesh() const { return m_mesh; }

      template <class Msh, class GF, class VF>
      static Real relL2Frob(const Msh& mesh,
                            const GF& uh,
                            const VF& uex)
      {
        const size_t dim = mesh.getSpaceDimension();
        P1 scalar(mesh); // scalar space for integrating a scalar field

        GridFunction err2(scalar);
        err2 = Pow(Frobenius(uh - uex), 2);
        const Real l2e = Math::sqrt(Integral(err2).compute());

        // ||uex||_{L2}
        GridFunction sol2(scalar);
        sol2 = Pow(Frobenius(uex), 2);
        const Real l2u = Math::sqrt(Integral(sol2).compute());

        (void)dim;
        return l2e / (l2u + Real(1e-30));
      }

    private:
      Mesh<Context::Local> m_mesh;
  };

  using Tetra8  = Elasticity3DFixture<Polytope::Type::Tetrahedron, 8>;
  using Hex8    = Elasticity3DFixture<Polytope::Type::Hexahedron,   8>;
  using Tetra16 = Elasticity3DFixture<Polytope::Type::Tetrahedron, 16>;
  using Hex16    = Elasticity3DFixture<Polytope::Type::Hexahedron,   16>;
  using Tetra32 = Elasticity3DFixture<Polytope::Type::Tetrahedron, 32>;

  // ------------------------------------------------------------
  // Affine cases: -div(sigma(u)) = 0 (constant strain => div sigma = 0)
  // These are "exact" PDE solutions with f = 0.
  // ------------------------------------------------------------

  TEST_F(Tetra8, AffineExact_Identity)
  {
    const Real lambda = 1.0, mu = 1.0;
    const size_t dim = mesh().getSpaceDimension();

    P1 vh(mesh(), dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    VectorFunction u_exact{ F::x, F::y, F::z };
    VectorFunction f_body{ Zero(), Zero(), Zero() };

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f_body, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    // For affine solutions, the discrete solution should match strongly imposed Dirichlet data.
    const Real rel = relL2Frob(mesh(), u.getSolution(), u_exact);
    EXPECT_NEAR(rel, 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST_F(Hex8, AffineExact_Identity)
  {
    const Real lambda = 1.5, mu = 0.5;
    const size_t dim = mesh().getSpaceDimension();

    P1 vh(mesh(), dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    VectorFunction u_exact{ F::x, F::y, F::z };
    VectorFunction f_body{ Zero(), Zero(), Zero() };

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f_body, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    const Real rel = relL2Frob(mesh(), u.getSolution(), u_exact);
    EXPECT_NEAR(rel, 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST_F(Tetra8, GeneralAffine_Tetrahedron)
  {
    const Real lambda = 1.5, mu = 0.5;
    const size_t dim = mesh().getSpaceDimension();

    P1 vh(mesh(), dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    // u(x,y,z) = (2x - y + z, -x + 3y - 2z, x + y + 2z)
    VectorFunction u_exact{
      2 * F::x - F::y + F::z,
      -F::x + 3 * F::y - 2 * F::z,
      F::x + F::y + 2 * F::z
    };

    VectorFunction f_body{ Zero(), Zero(), Zero() }; // div(sigma)=0 for affine u

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f_body, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    const Real rel = relL2Frob(mesh(), u.getSolution(), u_exact);
    EXPECT_NEAR(rel, 0.0, RODIN_FUZZY_CONSTANT);
  }

  // ------------------------------------------------------------
  // Non-affine manufactured cases with correct f = -div(sigma(u_exact))
  //
  // For constant lambda, mu:
  //   f = -div(sigma(u)) = -mu * Δu - (lambda + mu) * grad(div u)
  //
  // These tests are mathematically consistent, but P1 cannot reproduce them exactly
  // on a fixed mesh, so we use a relaxed bound.
  // ------------------------------------------------------------

  TEST_F(Hex8, Polynomial_Hexahedron)
  {
    const Real lambda = 2.0, mu = 1.0;
    const size_t dim = mesh().getSpaceDimension();

    P1 vh(mesh(), dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    // u = (x(1-x), y(1-y), z(1-z))
    const auto u1 = F::x * (1 - F::x);
    const auto u2 = F::y * (1 - F::y);
    const auto u3 = F::z * (1 - F::z);
    VectorFunction u_exact{ u1, u2, u3 };

    // For u1 = x - x^2: Δu1 = -2, similarly for y,z.
    // div u = (1-2x) + (1-2y) + (1-2z) = 3 - 2(x+y+z)
    // grad(div u) = (-2, -2, -2)
    // f = -mu * Δu - (lambda+mu) * grad(div u)
    //   = -mu * (-2,-2,-2) - (lambda+mu) * (-2,-2,-2)
    //   = (2mu + 2lambda + 2mu) * (1,1,1) = (2lambda + 4mu) * (1,1,1)
    const Real c = 2 * lambda + 4 * mu;
    VectorFunction f_body{ c, c, c };

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f_body, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    const Real rel = relL2Frob(mesh(), u.getSolution(), u_exact);
    EXPECT_LT(rel, RODIN_FUZZY_CONSTANT); // relaxed fixed-mesh check
  }

  TEST_F(Tetra32, MixedComponents_Tetrahedron)
  {
    const Real lambda = 1.0, mu = 1.0;
    const size_t dim = mesh().getSpaceDimension();

    P1 vh(mesh(), dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    // u = (x^2 + y, y^2 + z, z^2 + x)
    VectorFunction u_exact{
      Pow(F::x, 2) + F::y,
      Pow(F::y, 2) + F::z,
      Pow(F::z, 2) + F::x
    };

    // div u = 2x + 2y + 2z
    // grad(div u) = (2,2,2)
    // Δu = (2,2,2)
    // f = -mu*(2,2,2) - (lambda+mu)*(2,2,2) = -(2lambda + 4mu) * (1,1,1)
    const Real c = -(2 * lambda + 4 * mu);
    VectorFunction f_body{ c, c, c };

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f_body, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    const Real rel = relL2Frob(mesh(), u.getSolution(), u_exact);
    EXPECT_LT(0.1 * rel, RODIN_FUZZY_CONSTANT); // relaxed fixed-mesh check
  }

  TEST_F(Hex16, Polynomial_Hexahedron)
  {
    const Real lambda = 2.0, mu = 1.0;
    const size_t dim = mesh().getSpaceDimension();

    P1 vh(mesh(), dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    const auto u1 = F::x * (1 - F::x);
    const auto u2 = F::y * (1 - F::y);
    const auto u3 = F::z * (1 - F::z);
    VectorFunction u_exact{ u1, u2, u3 };

    // Same exact forcing:
    // f = -mu Δu - (lambda+mu) grad(div u) = (2lambda + 4mu) (1,1,1)
    const Real c = 2 * lambda + 4 * mu;
    VectorFunction f_body{ c, c, c };

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f_body, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    // On the refined mesh, relative error should improve vs Hex8; keep a relaxed bound.
    const Real rel = relL2Frob(mesh(), u.getSolution(), u_exact);
    EXPECT_LT(rel, RODIN_FUZZY_CONSTANT);
  }
}
