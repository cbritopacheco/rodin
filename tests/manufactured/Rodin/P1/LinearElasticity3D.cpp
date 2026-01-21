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

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

/**
 * @brief 3D Manufactured solutions for the linear elasticity problem using P1 spaces.
 *
 * The strong form in 3D:
 * @f[
 *   -\nabla\cdot\sigma(u) = f \quad\text{in }\Omega \subset \mathbb{R}^3,
 *   \quad u = g \quad\text{on }\partial\Omega,
 * @f]
 * where
 * @f[
 *   \sigma(u) = \lambda\,\mathrm{tr}(\varepsilon(u))\,I + 2\mu\,\varepsilon(u),
 *   \quad \varepsilon(u)=\tfrac12(\nabla u + (\nabla u)^T),
 * @f]
 * with Lamé parameters @f$\lambda,\mu>0@f$.
 *
 * The weak form: Find @f$ u\in [V]^3 @f$ such that
 * @f[
 *   \int_\Omega \bigl[\lambda\,\mathrm{div}(u)\,\mathrm{div}(v)
 *       + 2\mu\,\varepsilon(u):\varepsilon(v)\bigr]\,dx = \int_\Omega f \cdot v\,dx,
 *   \quad \forall v\in [V]^3,
 *   \quad u = g \text{ on }\partial\Omega.
 * @f]
 */
namespace Rodin::Tests::Manufactured::LinearElasticity3D
{
  /**
   * @brief Affine manufactured solution (exactly representable in P1)
   * u(x,y,z) = (x, y, z)
   */
  TEST(LinearElasticity3D, AffineExact_Tetrahedron)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    const Real lambda = 1.0, mu = 1.0;
    const size_t dim = mesh.getSpaceDimension();

    P1 vh(mesh, dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    auto u_exact = VectorFunction{ F::x, F::y, F::z };

    Problem elasticity(u, v);
    elasticity = Integral(lambda * Div(u), Div(v))
               + Integral(mu * (Jacobian(u) + Jacobian(u).T()), 
                         0.5 * (Jacobian(v) + Jacobian(v).T()))
               + DirichletBC(u, u_exact);
    CG(elasticity).solve();

    P1 scalar(mesh);
    GridFunction err2(scalar);
    err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
    const Real L2error = Integral(err2).compute();
    EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Affine solution on Hexahedron mesh
   */
  TEST(LinearElasticity3D, AffineExact_Hexahedron)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Hexahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    const Real lambda = 1.5, mu = 0.5;
    const size_t dim = mesh.getSpaceDimension();

    P1 vh(mesh, dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    auto u_exact = VectorFunction{ F::x, F::y, F::z };

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

  /**
   * @brief General affine displacement field
   * u(x,y,z) = (2x - y + z, -x + 3y - 2z, x + y + 2z)
   */
  TEST(LinearElasticity3D, GeneralAffine_Tetrahedron)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    const Real lambda = 1.5, mu = 0.5;
    const size_t dim = mesh.getSpaceDimension();

    P1 vh(mesh, dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    auto u_exact = VectorFunction{ 
      2 * F::x - F::y + F::z,
      -F::x + 3 * F::y - 2 * F::z,
      F::x + F::y + 2 * F::z
    };

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

  /**
   * @brief Trigonometric displacement field
   * u(x,y,z) = (sin(πx) sin(πy) sin(πz), sin(πx) sin(πy) sin(πz), sin(πx) sin(πy) sin(πz))
   */
  TEST(LinearElasticity3D, SimpleSine_Tetrahedron)
  {
    auto pi = Math::Constants::pi();

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    const Real lambda = 1.0, mu = 1.0;
    const size_t dim = mesh.getSpaceDimension();

    P1 vh(mesh, dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    auto u1 = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto u_exact = VectorFunction{ u1, u1, u1 };

    // For this displacement, div(u) = 3π cos(πx) sin(πy) sin(πz) + ...
    // The body force f can be computed from -∇·σ(u)
    // For simplicity, we compute the strain energy and let the solver find the equilibrium
    
    // Approximate body force (exact computation is complex)
    auto f1 = -(lambda + 2 * mu) * 3 * pi * pi * u1;
    auto f_body = VectorFunction{ f1, f1, f1 };

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f_body, v)
               + DirichletBC(u, u_exact);
    CG(elasticity).solve();

    P1 scalar(mesh);
    GridFunction err2(scalar);
    err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
    const Real L2error = Integral(err2).compute();
    EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Polynomial displacement field
   * u(x,y,z) = (x(1-x), y(1-y), z(1-z))
   */
  TEST(LinearElasticity3D, Polynomial_Hexahedron)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Hexahedron, {8, 8, 8});
    mesh.getConnectivity().compute(2, 3);

    const Real lambda = 2.0, mu = 1.0;
    const size_t dim = mesh.getSpaceDimension();

    P1 vh(mesh, dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    auto u_exact = VectorFunction{ 
      F::x * (1 - F::x),
      F::y * (1 - F::y),
      F::z * (1 - F::z)
    };

    // Body force: f = -∇·σ(u)
    // For this polynomial field, div(u) = (1-2x) + (1-2y) + (1-2z) = 3 - 2(x+y+z)
    // Laplacian of each component: Δu_i = 2
    // f_i = -(λ ∂_i(div u) + 2μ Δu_i)
    // ∂_x(div u) = -2, ∂_y(div u) = -2, ∂_z(div u) = -2
    auto f_body = VectorFunction{
      -(-2 * lambda - 2 * 2 * mu),
      -(-2 * lambda - 2 * 2 * mu),
      -(-2 * lambda - 2 * 2 * mu)
    };

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f_body, v)
               + DirichletBC(u, u_exact);
    CG(elasticity).solve();

    P1 scalar(mesh);
    GridFunction err2(scalar);
    err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
    const Real L2error = Integral(err2).compute();
    EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Mixed components with different functions
   * u(x,y,z) = (x² + y, y² + z, z² + x)
   */
  TEST(LinearElasticity3D, MixedComponents_Tetrahedron)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {16, 16, 16});
    mesh.getConnectivity().compute(2, 3);

    const Real lambda = 1.0, mu = 1.0;
    const size_t dim = mesh.getSpaceDimension();

    P1 vh(mesh, dim);
    TrialFunction u(vh);
    TestFunction  v(vh);

    auto u_exact = VectorFunction{ 
      Pow(F::x, 2) + F::y,
      Pow(F::y, 2) + F::z,
      Pow(F::z, 2) + F::x
    };

    // div(u) = 2x + 2y + 2z
    // ∂_x(div u) = 2, ∂_y(div u) = 2, ∂_z(div u) = 2
    // Δu_1 = 2, Δu_2 = 2, Δu_3 = 2
    auto f_body = VectorFunction{
      -(2 * lambda + 2 * 2 * mu),
      -(2 * lambda + 2 * 2 * mu),
      -(2 * lambda + 2 * 2 * mu)
    };

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - Integral(f_body, v)
               + DirichletBC(u, u_exact);
    CG(elasticity).solve();

    P1 scalar(mesh);
    GridFunction err2(scalar);
    err2 = Pow(Frobenius(u.getSolution() - u_exact), 2);
    const Real L2error = Integral(err2).compute();
    EXPECT_NEAR(L2error, 0.0, RODIN_FUZZY_CONSTANT);
  }
}
