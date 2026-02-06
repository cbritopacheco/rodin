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

/**
 * @brief Manufactured solutions for the linear elasticity problem.
 *
 * The strong form is:
 * @f[
 *   -\nabla\cdot\sigma(u) = 0 \quad\text{in }\Omega,
 *   \quad u = g \quad\text{on }\partial\Omega,
 * @f]
 * where
 * @f[
 *   \sigma(u) = \lambda\,\mathrm{tr}(\varepsilon(u))\,I + 2\mu\,\varepsilon(u),
 *   \quad \varepsilon(u)=\tfrac12(\nabla u + (\nabla u)^T),
 * @f]
 * with Lamé parameters @f$\lambda,\mu>0@f$.  The weak form: Find
 * @f$ u\in [V]^2 @f$ such that
 * @f[
 *   \int_\Omega \bigl[\lambda\,\mathrm{div}(u)\,\mathrm{div}(v)
 *       + 2\mu\,\varepsilon(u):\varepsilon(v)\bigr]\,dx = 0,
 *   \quad \forall v\in [V]^2,
 *   \quad u = g \text{ on }\partial\Omega.
 * @f]
 */
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

  /**
   * @brief Affine manufactured solution
   *
   * Domain and boundary:
   * @f[ \Omega = [0,1]\times[0,1], \quad \Gamma = \partial\Omega. @f]
   *
   * Lamé parameters:
   * @f[ \lambda = 1.0, \quad \mu = 1.0. @f]
   *
   * Manufactured solution (displacement):
   * @f[ u_{\rm exact}(x,y) = \begin{pmatrix} x \\ y \end{pmatrix}. @f]
   *
   * Boundary data:
   * @f[ g(x,y)=u_{\rm exact}(x,y). @f]
   *
   * No body forces (zero RHS).
   *
   * Checks that P1 elements reproduce affine fields exactly:
   * @f[ \|u_h - u_{\rm exact}\|_{L^2}=0. @f]
   */
  TEST_P(Manufactured_LinearElasticity_Test_16x16, AffineExact)
  {
    Mesh mesh = this->getMesh();
    const Real lambda = 1.0, mu = 1.0;
    const size_t dim = mesh.getSpaceDimension();
    P1 vh(mesh, dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    auto u_exact = VectorFunction{ F::x, F::y };

    {
      Problem elasticity(u, v);
      elasticity = Integral(lambda * Div(u), Div(v))
                 + Integral(
                     mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()))
                 + DirichletBC(u, u_exact);
      CG(elasticity).solve();

      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      const Real L2error = Integral(err2).compute();
      EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Problem elasticity(u, v);
      elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
                 + DirichletBC(u, u_exact);
      CG(elasticity).solve();

      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      const Real L2error = Integral(err2).compute();
      EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  /**
   * @brief  A general affine manufactured solution exactly representable in P1.
   *
   * Domain:
   * @f[
   *   \Omega = [0,1]\times[0,1],
   * \quad
   *   \Gamma = \partial\Omega.
   * @f]
   *
   * Lamé parameters:
   * @f[ \lambda = 1.5,\quad \mu = 0.5. @f]
   *
   * Exact displacement (affine):
   * @f[
   *   u_{\rm exact}(x,y)
   *   = \begin{pmatrix}
   *       2x - y + 1 \\[4pt]
   *       -x + 3y - 2
   *     \end{pmatrix}.
   * @f]
   *
   * Since @f$u_{\rm exact}\in [\mathbb{P}_1]^2@f$, the stress field
   * @f$\sigma(u_{\rm exact})@f$ is constant and @f$f = -\nabla\!\cdot\sigma = 0.@f$
   *
   * Boundary data:
   * @f[ g(x,y) = u_{\rm exact}(x,y). @f]
   *
   * We verify that the discrete @f$ L^2 @f$–error vanishes exactly on a
   * @f$ 16 \times 16 @f$ mesh.
   */
  TEST_P(Manufactured_LinearElasticity_Test_16x16, AffineGeneral)
  {
    Mesh mesh = this->getMesh();
    const Real lambda = 1.5, mu = 0.5;
    const size_t dim = mesh.getSpaceDimension();
    P1 vh(mesh, dim);

    // manufactured exact solution
    VectorFunction u_exact
    {
      2 * F::x - F::y + 1,
      -F::x + 3 * F::y - 2
    };

    // zero body force since divergence of constant stress = 0
    VectorFunction f = VectorFunction{ Zero(), Zero() };

    // assemble and solve
    TrialFunction u(vh);
    TestFunction  v(vh);

    {
      Problem elasticity(u, v);
      elasticity = Integral(lambda * Div(u), Div(v))
                 + Integral(
                     mu * (Jacobian(u) + Jacobian(u).T()),
                     0.5 * (Jacobian(v) + Jacobian(v).T()))
                 - Integral(f, v)
                 + DirichletBC(u, u_exact);
      CG(elasticity).solve();

      // compute L2 error
      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      const Real L2error = Integral(err2).compute();
      EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      Problem elasticity(u, v);
      elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
                 + DirichletBC(u, u_exact);
      CG(elasticity).solve();

      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      const Real L2error = Integral(err2).compute();
      EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  /**
   * @brief  A more complex manufactured solution mixing sine and cosine modes.
   *
   * Domain:
   * @f[
   *   \Omega = [0,1]\times[0,1],
   * \quad
   *   \Gamma = \partial\Omega.
   * @f]
   *
   * Exact displacement field:
   * @f[
   *   u_{\rm exact}(x,y)
   *   = \begin{pmatrix}
   *       \sin(\pi x)\,\cos(\pi y) \\[4pt]
   *       \cos(\pi x)\,\sin(\pi y)
   *     \end{pmatrix}.
   * @f]
   *
   * Lamé parameters:
   * @f[ \lambda = 1.0,\quad \mu = 1.0. @f]
   *
   * Body force @f$f = -\nabla\!\cdot\sigma(u_{\rm exact})@f$ can be shown to be
   * @f[
   *   f(x,y)
   *   = (2\lambda + 4\mu)\,\pi^2
   *     \begin{pmatrix}
   *       \sin(\pi x)\,\cos(\pi y) \\[4pt]
   *       \cos(\pi x)\,\sin(\pi y)
   *     \end{pmatrix}.
   * @f]
   *
   * Boundary data:
   * @f[ g(x,y) = u_{\rm exact}(x,y). @f]
   *
   * We verify that the discrete @f$L^2@f$–error of the computed solution
   * vanishes (up to RODIN_FUZZY_CONSTANT) on a @f$16\times16@f$ mesh.
   */
  TEST_P(Manufactured_LinearElasticity_Test_16x16, MixedTrigonometric)
  {
    // mesh + function space
    Mesh mesh = this->getMesh();
    const Real lambda = 1.0, mu = 1.0;
    const size_t dim = mesh.getSpaceDimension();
    P1 vh(mesh, dim);

    // exact solution
    auto pi = Rodin::Math::Constants::pi();
    VectorFunction u_exact
    {
      sin(pi * F::x) * cos(pi * F::y),
      cos(pi * F::x) * sin(pi * F::y)
    };

    // manufactured body force
    const Real coeff = (2 * lambda + 4 * mu) * pi * pi;
    VectorFunction f
    {
      coeff * sin(pi * F::x) * cos(pi * F::y),
      coeff * cos(pi * F::x) * sin(pi * F::y)
    };

    {
      // variational problem
      TrialFunction u(vh);
      TestFunction  v(vh);
      Problem elasticity(u, v);
      elasticity = Integral(lambda * Div(u), Div(v))
                 + Integral(
                     mu * (Jacobian(u) + Jacobian(u).T()),
                     0.5 * (Jacobian(v) + Jacobian(v).T()))
                 - Integral(f, v)
                 + DirichletBC(u, u_exact);
      CG(elasticity).solve();

      // compute L^2 error
      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      const Real L2error = Integral(err2).compute();
      EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      // variational problem
      TrialFunction u(vh);
      TestFunction  v(vh);
      Problem elasticity(u, v);
      elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
                 - Integral(f, v)
                 + DirichletBC(u, u_exact);
      CG(elasticity).solve();

      // compute L^2 error
      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      const Real L2error = Integral(err2).compute();
      EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  /**
   * @brief  A highly nontrivial manufactured solution combining polynomial and trigonometric components.
   *
   * Domain:
   * @f[
   *   \Omega = [0,1]\times[0,1],
   * \quad
   *   \Gamma = \partial\Omega.
   * @f]
   *
   * Exact displacement field:
   * @f[
   *   u_{\rm exact}(x,y)
   *   = \begin{pmatrix}
   *       x^2(1-x)\,\sin(\pi y) \\[4pt]
   *       y^2(1-y)\,\cos(\pi x)
   *     \end{pmatrix}.
   * @f]
   *
   * Lamé parameters:
   * @f[ \lambda = 2.5,\quad \mu = 0.75. @f]
   *
   * We compute the body force 
   * @f[
   *   f(x,y) = -\nabla\!\cdot\sigma\bigl(u_{\rm exact}(x,y)\bigr),
   * @f]
   * which works out to the vector
   * @f[
   *   f(x,y)
   *   = \begin{pmatrix}
   *       -\bigl[\lambda\,\partial_x(\nabla\!\cdot u_{\rm exact})
   *         + 2\mu\,\partial_x\varepsilon_{11}(u_{\rm exact})
   *         + 2\mu\,\partial_y\varepsilon_{12}(u_{\rm exact})\bigr] \\[6pt]
   *       -\bigl[\lambda\,\partial_y(\nabla\!\cdot u_{\rm exact})
   *         + 2\mu\,\partial_x\varepsilon_{12}(u_{\rm exact})
   *         + 2\mu\,\partial_y\varepsilon_{22}(u_{\rm exact})\bigr]
   *     \end{pmatrix},
   * @f]
   * where each derivative can be expanded analytically.
   *
   * Boundary data:
   * @f[
   *   g(x,y) = u_{\rm exact}(x,y).
   * @f]
   *
   * We verify that the discrete @f$L^2@f$–error vanishes (up to RODIN_FUZZY_CONSTANT)
   * on a @f$32\times32@f$ mesh of either triangles or quads.
   */
  TEST_P(Manufactured_LinearElasticity_Test_32x32, ComplexPolyTrig)
  {
    Mesh mesh = this->getMesh();
    const Real lambda = 2.5, mu = 0.75;
    const size_t dim = mesh.getSpaceDimension();
    P1 vh(mesh, dim);

    // exact solution
    auto pi = Rodin::Math::Constants::pi();
    VectorFunction u_exact
    {
      // p(x) * sin(pi y)
      F::x*F::x*(1 - F::x)*sin(pi * F::y),
      // q(y) * cos(pi x)
      F::y*F::y*(1 - F::y)*cos(pi * F::x)
    };

    // helper polynomials
    auto p   =  F::x*F::x*(1 - F::x);
    auto dp  = 2*F::x - 3*F::x*F::x;   // p'(x)
    auto ddp = 2    - 6*F::x;          // p''(x)
    auto q   =  F::y*F::y*(1 - F::y);
    auto dq  = 2*F::y - 3*F::y*F::y;   // q'(y)
    auto ddq = 2    - 6*F::y;          // q''(y)

    // manufactured body force f = -div σ(u_exact)
    VectorFunction f
    {
      -(
        (lambda + 2*mu)*ddp * sin(pi * F::y)
        - mu * pi*pi * p   * sin(pi * F::y)
        - (lambda + mu) * pi * dq * sin(pi * F::x)
      ),
      -(
        (lambda + 2*mu)*ddq * cos(pi * F::x)
        - mu * pi*pi * q   * cos(pi * F::x)
        + (lambda + mu) * pi * dp * cos(pi * F::y)
      )
    };

    {
      TrialFunction u(vh);
      TestFunction  v(vh);
      Problem elasticity(u, v);
      elasticity = Integral(lambda * Div(u), Div(v))
                 + Integral(
                     mu * (Jacobian(u) + Jacobian(u).T()),
                     0.5 * (Jacobian(v) + Jacobian(v).T()))
                 - Integral(f, v)
                 + DirichletBC(u, u_exact);
      CG(elasticity).solve();

      // compute L2 error
      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      const Real L2error = Integral(err2).compute();
      EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      TrialFunction u(vh);
      TestFunction  v(vh);
      Problem elasticity(u, v);
      elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
                 - Integral(f, v)
                 + DirichletBC(u, u_exact);
      CG(elasticity).solve();

      // compute L^2 error
      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      const Real L2error = Integral(err2).compute();
      EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  /**
   * @brief  A general affine manufactured solution exactly representable in P1.
   *
   * Domain:
   * @f[
   *   \Omega = [0,1]\times[0,1],
   * \quad
   *   \Gamma = \partial\Omega.
   * @f]
   *
   * Lamé parameters:
   * @f[ \lambda = 1.5,\quad \mu = 0.5. @f]
   *
   * Exact displacement (affine):
   * @f[
   *   u_{\rm exact}(x,y)
   *   = \begin{pmatrix}
   *       2x - y + 1 \\[4pt]
   *       -x + 3y - 2
   *     \end{pmatrix}.
   * @f]
   *
   * Since @f$u_{\rm exact}\in [\mathbb{P}_1]^2@f$, the stress field
   * @f$\sigma(u_{\rm exact})@f$ is constant and @f$f = -\nabla\!\cdot\sigma = 0.@f$
   *
   * Boundary data:
   * @f[ g(x,y) = u_{\rm exact}(x,y). @f]
   *
   * We verify that the discrete @f$L^2@f$–error vanishes exactly on a
   * @f$16\times16@f$ mesh.
   */
  TEST_P(Manufactured_LinearElasticity_Test_16x16, AffineGeneralDisplaced)
  {
    Mesh mesh = this->getMesh();
    const Real lambda = 1.5, mu = 0.5;
    const size_t dim = mesh.getSpaceDimension();
    P1 vh(mesh, dim);

    // manufactured exact solution
    VectorFunction u_exact
    {
      2 * F::x - F::y + 1,
      -F::x + 3 * F::y - 2
    };

    // zero body force since divergence of constant stress = 0
    VectorFunction f = VectorFunction{ Zero(), Zero() };

    {
      // assemble and solve
      TrialFunction u(vh);
      TestFunction  v(vh);
      Problem elasticity(u, v);
      elasticity = Integral(lambda * Div(u), Div(v))
                 + Integral(
                     mu * (Jacobian(u) + Jacobian(u).T()),
                     0.5 * (Jacobian(v) + Jacobian(v).T()))
                 - Integral(f, v)
                 + DirichletBC(u, u_exact);
      CG(elasticity).solve();

      // compute L2 error
      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      const Real L2error = Integral(err2).compute();
      EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      TrialFunction u(vh);
      TestFunction  v(vh);
      Problem elasticity(u, v);
      elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
                 - Integral(f, v)
                 + DirichletBC(u, u_exact);
      CG(elasticity).solve();

      // compute L^2 error
      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
      const Real L2error = Integral(err2).compute();
      EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
    }
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    Manufactured_LinearElasticity_Test_16x16,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

  INSTANTIATE_TEST_SUITE_P(
    MeshParams32x32,
    Manufactured_LinearElasticity_Test_32x32,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}
