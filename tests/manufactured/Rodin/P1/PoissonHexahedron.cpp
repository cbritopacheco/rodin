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
 * @brief Manufactured solutions for the 3D Poisson problem on Hexahedron meshes.
 *
 * The system is given by:
 * @f[
 * \left\{
 * \begin{aligned}
 *   -\Delta u &= f \quad \text{in } \Omega,\\
 *  u &= g \quad \text{on } \partial\Omega.
 * \end{aligned}
 * \right.
 * @f]
 *
 * The weak formulation is: Find @f$ u \in V @f$ such that
 * @f[
 *   \int_\Omega \nabla u \cdot \nabla v \,dx = \int_\Omega f \, v \,dx,
 * @f]
 * for all @f$ v \in V @f$, with the essential boundary condition
 * @f[
 *   u = g \quad \text{on } \partial\Omega.
 * @f]
 */
namespace Rodin::Tests::Manufactured::PoissonHexahedron
{
  template <size_t M>
  class Manufactured_PoissonHexahedron_Test : public ::testing::Test
  {
    protected:
      Mesh<Context::Local> getMesh()
      {
        Mesh mesh;
        mesh = mesh.UniformGrid(Polytope::Type::Hexahedron, { M, M, M });
        mesh.scale(1.0 / (M - 1));
        mesh.getConnectivity().compute(2, 3);
        return mesh;
      }
  };

  using Manufactured_PoissonHexahedron_Test_8 =
    Rodin::Tests::Manufactured::PoissonHexahedron::Manufactured_PoissonHexahedron_Test<8>;
  using Manufactured_PoissonHexahedron_Test_16 =
    Rodin::Tests::Manufactured::PoissonHexahedron::Manufactured_PoissonHexahedron_Test<16>;
  using Manufactured_PoissonHexahedron_Test_32 =
    Rodin::Tests::Manufactured::PoissonHexahedron::Manufactured_PoissonHexahedron_Test<32>;

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  \Gamma = \partial \Omega
   * @f]
   *
   * @f[
   *  f(x, y, z) = 3 \pi^2 \sin(\pi x) \sin(\pi y) \sin(\pi z)
   * @f]
   *
   * @f[
   *  g(x, y, z) = 0
   * @f]
   *
   * @f[
   *  u(x, y, z) = \sin(\pi x) \sin(\pi y) \sin(\pi z)
   * @f]
   *
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_8, PoissonHexahedron_SimpleSine)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 3 * pi * pi * sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v);
    poisson.assemble();

    // Essential boundary conditions: u = 0 on boundary
    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u_exact = sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  f(x, y, z) = -6 x (1 - x) - 6 y (1 - y) - 6 z (1 - z)
   * @f]
   *
   * @f[
   *  g(x, y, z) = 0
   * @f]
   *
   * @f[
   *  u(x, y, z) = x(1 - x) y(1 - y) z(1 - z)
   * @f]
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_8, PoissonHexahedron_Polynomial)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = -6 * F::x * (1 - F::x) - 6 * F::y * (1 - F::y) - 6 * F::z * (1 - F::z);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v);
    poisson.assemble();

    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u_exact = F::x * (1 - F::x) * F::y * (1 - F::y) * F::z * (1 - F::z);

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  f(x, y, z) = x (1 - x) \pi^2 \sin(\pi y) \sin(\pi z) - 2 \sin(\pi y) \sin(\pi z)
   * @f]
   *
   * @f[
   *  g(x, y, z) = 0
   * @f]
   *
   * @f[
   *  u(x, y, z) = x(1 - x) \sin(\pi y) \sin(\pi z)
   * @f]
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_8, PoissonHexahedron_TrigonometricPolynomial)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = F::x * (1 - F::x) * pi * pi * sin(pi * F::y) * sin(pi * F::z) - 2 * sin(pi * F::y) * sin(pi * F::z);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v);
    poisson.assemble();

    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u_exact = F::x * (1 - F::x) * sin(pi * F::y) * sin(pi * F::z);

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  f(x, y, z) = 3 \pi^2 \cos(\pi x) \cos(\pi y) \cos(\pi z)
   * @f]
   *
   * @f[
   *  g(x, y, z) = \cos(\pi x) \cos(\pi y) \cos(\pi z)
   * @f]
   *
   * @f[
   *  u(x, y, z) = \cos(\pi x) \cos(\pi y) \cos(\pi z)
   * @f]
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_16, PoissonHexahedron_NonhomogeneousDirichlet)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 3 * pi * pi * cos(pi * F::x) * cos(pi * F::y) * cos(pi * F::z);
    auto g = cos(pi * F::x) * cos(pi * F::y) * cos(pi * F::z);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, g);
    poisson.assemble();

    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u_exact = cos(pi * F::x) * cos(pi * F::y) * cos(pi * F::z);

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  f(x, y, z) = (2 \pi^2 + 1) \sin(\pi x) \sin(\pi y) e^z
   * @f]
   *
   * @f[
   *  g(x, y, z) = \sin(\pi x) \sin(\pi y) e^z
   * @f]
   *
   * @f[
   *  u(x, y, z) = \sin(\pi x) \sin(\pi y) e^z
   * @f]
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_16, PoissonHexahedron_MixedBoundary)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = (2 * pi * pi + 1) * sin(pi * F::x) * sin(pi * F::y) * exp(F::z);
    auto g = sin(pi * F::x) * sin(pi * F::y) * exp(F::z);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, g);
    poisson.assemble();

    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u_exact = sin(pi * F::x) * sin(pi * F::y) * exp(F::z);

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  f(x, y, z) = 0
   * @f]
   *
   * @f[
   *  g(x, y, z) = x + y + z
   * @f]
   *
   * @f[
   *  u(x, y, z) = x + y + z
   * @f]
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_8, PoissonHexahedron_LinearNonhomogeneous)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = ScalarFunction(0.0);
    auto g = F::x + F::y + F::z;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, g);
    poisson.assemble();

    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u_exact = F::x + F::y + F::z;

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }

  // ===================== Vector Poisson Tests =====================

  /**
   * Vector Poisson with all components sin(pi x) sin(pi y) sin(pi z)
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_16, VectorPoissonHexahedron_SimpleSine)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh, 3);  // 3D vector space
    auto f = 3 * pi * pi * VectorFunction({
      sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z),
      sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z),
      sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z)
    });

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(Dot(f, v));
    poisson.assemble();

    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u_exact = VectorFunction({
      sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z),
      sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z),
      sin(pi * F::x) * sin(pi * F::y) * sin(pi * F::z)
    });

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }

  /**
   * Vector Poisson with polynomial components
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_8, VectorPoissonHexahedron_Polynomial)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh, 3);
    auto f_component = -6 * F::x * (1 - F::x) - 6 * F::y * (1 - F::y) - 6 * F::z * (1 - F::z);
    auto f = VectorFunction({f_component, f_component, f_component});

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(Dot(f, v));
    poisson.assemble();

    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u_component = F::x * (1 - F::x) * F::y * (1 - F::y) * F::z * (1 - F::z);
    auto u_exact = VectorFunction({u_component, u_component, u_component});

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }

  /**
   * Vector Poisson with mixed trigonometric-polynomial components
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_16, VectorPoissonHexahedron_TrigonometricPolynomial)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh, 3);
    
    auto f1 = F::x * (1 - F::x) * pi * pi * sin(pi * F::y) * sin(pi * F::z) - 2 * sin(pi * F::y) * sin(pi * F::z);
    auto f2 = F::y * (1 - F::y) * pi * pi * sin(pi * F::x) * sin(pi * F::z) - 2 * sin(pi * F::x) * sin(pi * F::z);
    auto f3 = F::z * (1 - F::z) * pi * pi * sin(pi * F::x) * sin(pi * F::y) - 2 * sin(pi * F::x) * sin(pi * F::y);
    auto f = VectorFunction({f1, f2, f3});

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(Dot(f, v));
    poisson.assemble();

    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u1 = F::x * (1 - F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto u2 = F::y * (1 - F::y) * sin(pi * F::x) * sin(pi * F::z);
    auto u3 = F::z * (1 - F::z) * sin(pi * F::x) * sin(pi * F::y);
    auto u_exact = VectorFunction({u1, u2, u3});

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }

  /**
   * Vector Poisson with non-homogeneous Dirichlet boundary conditions
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_16, VectorPoissonHexahedron_NonhomogeneousDirichlet)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh, 3);
    
    auto f1 = 3 * pi * pi * cos(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto f2 = 3 * pi * pi * sin(pi * F::x) * cos(pi * F::y) * sin(pi * F::z);
    auto f3 = 3 * pi * pi * sin(pi * F::x) * sin(pi * F::y) * cos(pi * F::z);
    auto f = VectorFunction({f1, f2, f3});
    
    auto g1 = cos(pi * F::x) * sin(pi * F::y) * sin(pi * F::z);
    auto g2 = sin(pi * F::x) * cos(pi * F::y) * sin(pi * F::z);
    auto g3 = sin(pi * F::x) * sin(pi * F::y) * cos(pi * F::z);
    auto g = VectorFunction({g1, g2, g3});

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(Dot(f, v)) + DirichletBC(u, g);
    poisson.assemble();

    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u_exact = VectorFunction({g1, g2, g3});

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }

  /**
   * Vector Poisson with quadratic components
   */
  TEST_F(Manufactured_PoissonHexahedron_Test_8, VectorPoissonHexahedron_Quadratic)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh, 3);
    
    auto f1 = ScalarFunction(-2.0);  // -d²/dx²(x²) = -2
    auto f2 = ScalarFunction(-2.0);  // -d²/dy²(y²) = -2
    auto f3 = ScalarFunction(-2.0);  // -d²/dz²(z²) = -2
    auto f = VectorFunction({f1, f2, f3});
    
    auto g1 = F::x * F::x;
    auto g2 = F::y * F::y;
    auto g3 = F::z * F::z;
    auto g = VectorFunction({g1, g2, g3});

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(Dot(f, v)) + DirichletBC(u, g);
    poisson.assemble();

    CG cg;
    poisson.solve(cg);

    GridFunction uh(vh);
    uh = poisson.getSolution();

    // Manufactured solution
    auto u_exact = VectorFunction({g1, g2, g3});

    // Compute error
    GridFunction error(vh);
    error.project(u_exact - uh);
    Real l2_error = error.getL2Norm();

    EXPECT_LT(l2_error, RODIN_FUZZY_CONSTANT);
  }
}
