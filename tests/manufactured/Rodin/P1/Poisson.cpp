/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Test/Random.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;
using namespace Rodin::Test::Random;


/**
 * @brief Manufactured solutions for the Poisson problem using P1 spaces.
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
namespace Rodin::Tests::Manufactured::Poisson
{
  template <size_t M>
  class ManufacturedPoissonTest : public ::testing::TestWithParam<Polytope::Type>
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

  using ManufacturedPoissonTest16x16 =
    Rodin::Tests::Manufactured::Poisson::ManufacturedPoissonTest<16>;
  using ManufacturedPoissonTest32x32 =
    Rodin::Tests::Manufactured::Poisson::ManufacturedPoissonTest<32>;
  using ManufacturedPoissonTest64x64 =
    Rodin::Tests::Manufactured::Poisson::ManufacturedPoissonTest<64>;

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  \Gamma = \partial \Omega
   * @f]
   *
   * @f[
   *  f(x, y) = 2 \pi^2 \sin(\pi x) \sin(\pi y)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   * @f[
   *  u(x, y) = \sin(\pi x) \sin(\pi y)
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_SimpleSine)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  \Gamma = \partial \Omega
   * @f]
   *
   * @f[
   *  f(x, y) = \left[ (\omega_1 \pi)^2 + (\omega_2 \pi)^2 \right] \sin(\omega_1 \pi x) \sin(\omega_2 \pi y)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   * @f[
   *  u(x,y) = \sin(\omega_1 \pi x) \sin(\omega_2 \pi y)
   * @f]
   *
   * @f[
   *  \omega_1, \omega_2 \in \mathbb{N}
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest64x64, Poisson_VariableFrequency)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    // For reproducibility, we choose fixed frequencies.
    const Real omega1 = 3;
    const Real omega2 = 2;

    P1 vh(mesh);
    auto f = (omega1 * pi * omega1 * pi + omega2 * pi * omega2 * pi) *
             sin(omega1 * pi * F::x) * sin(omega2 * pi * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto solution = sin(omega1 * pi * F::x) * sin(omega2 * pi * F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y) = x (1-x) \sin(2 \pi y)
   * @f]
   *
   * @f[
   *  f(x, y) = 2 \sin(2 \pi y) + 4 \pi^2 x (1-x) \sin(2 \pi y)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_VariableAmplitude)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 2 * sin(2 * pi * F::y)
           + 4 * pi * pi * F::x * (1 - F::x) * sin(2 * pi * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto solution = F::x * (1 - F::x) * sin(2 * pi * F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y) = x (1-x) y (1-y)
   * @f]
   *
   * @f[
   *  f(x, y) = 2 y (1-y) + 2 x (1-x)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_Polynomial)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 2 * F::y * (1 - F::y) + 2 * F::x * (1 - F::x);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto solution = F::x * (1 - F::x) * F::y * (1 - F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y) = x (1-x) \sin(\pi y)
   * @f]
   *
   * @f[
   *  f(x, y) = 2 \sin(\pi y) + \pi^2 x (1-x) \sin(\pi y)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_TrigonometricPolynomial)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 2 * sin(pi * F::y)
           + pi * pi * F::x * (1 - F::x) * sin(pi * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto solution = F::x * (1 - F::x) * sin(pi * F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y) = \cos(\pi x) \cos(\pi y)
   * @f]
   *
   * @f[
   *  f(x, y) = 2 \pi^2 \cos(\pi x) \cos(\pi y)
   * @f]
   *
   * @f[
   *  g(x, y) = u(x,y)
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_NonhomogeneousDirichlet)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = 2 * pi * pi * cos(pi * F::x) * cos(pi * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    // Use the manufactured solution as Dirichlet data.
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, cos(pi * F::x) * cos(pi * F::y));
    CG(poisson).solve();

    auto solution = cos(pi * F::x) * cos(pi * F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y) = \sin(\pi x) e^y
   * @f]
   *
   * @f[
   *  f(x, y) = (\pi^2-1) \sin(\pi x) e^y
   * @f]
   *
   * @f[
   *  g(x, y) = u(x,y)
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_MixedBoundary)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto f = (pi * pi - 1) * sin(pi * F::x) * exp(F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    // Apply Dirichlet conditions on the entire boundary.
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, sin(pi * F::x) * exp(F::y));
    CG(poisson).solve();

    auto solution = sin(pi * F::x) * exp(F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y) = e^{x(1-x)y(1-y)} - 1
   * @f]
   *
   * @f[
   *  f(x, y) = -e^{x(1-x)y(1-y)} \left(
   *    -2y(1-y)-2x(1-x) + (1-2x)^2 \left[y(1-y)\right]^2 + (1-2y)^2 \left[x(1-x)\right]^2
   *  \right)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_Exponential)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh);
    auto expr = F::x * (1 - F::x) * F::y * (1 - F::y);
    auto u_expr = exp(expr) - 1;
    auto f = - exp(expr) * (
               - 2 * F::y * (1 - F::y) - 2 * F::x * (1 - F::x)
               + pow((1 - 2 * F::x), 2) * pow(F::y * (1 - F::y), 2)
               + pow((1 - 2 * F::y), 2) * pow(F::x * (1 - F::x), 2)
             );

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto solution = u_expr;

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * Manufactured solution:
   * @f[
   *  u(x, y) = e^{\,x(1-x)y(1-y)}-1,
   * @f]
   *
   * Forcing function:
   * @f[
   *  f(x, y) = -e^{\,x(1-x)y(1-y)}\Bigl(
   *      -2\,y(1-y)-2\,x(1-x)
   *      +(1-2x)^2\,[y(1-y)]^2+(1-2y)^2\,[x(1-x)]^2
   *    \Bigr)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_Exponential_Revised)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh);

    // Define g(x,y) = x(1-x)y(1-y)
    auto g = F::x * (1 - F::x) * F::y * (1 - F::y);
    // Manufactured solution: u = e^{g} - 1  (so that u=0 on the boundary)
    auto u_expr = exp(g) - 1;

    // Forcing function f = -e^{g} (g_{xx}+g_{yy}+g_x^2+g_y^2)
    auto f_expr = - exp(g) * (
                     (-2 * F::y * (1 - F::y) - 2 * F::x * (1 - F::x))
                     + pow((1 - 2 * F::x), 2) * pow(F::y * (1 - F::y), 2)
                     + pow((1 - 2 * F::y), 2) * pow(F::x * (1 - F::x), 2)
                   );

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f_expr, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  A(x) = (1-x^2)\cos(3\pi x),\quad B(y) = (1-y^2)\sin(3\pi y)
   * @f]
   *
   * and define
   *
   * @f[
   *  u(x, y) = \begin{pmatrix}
   *    A(x)B(y)\\[1mm]
   *    A(x)B(y)
   *  \end{pmatrix}
   *  = \begin{pmatrix}
   *    (1-x^2)(1-y^2)\cos(3\pi x)\sin(3\pi y)\\[1mm]
   *    (1-x^2)(1-y^2)\cos(3\pi x)\sin(3\pi y)
   *  \end{pmatrix}
   * @f]
   *
   * @f[
   *  f(x, y) = \begin{pmatrix}
   *    -\Bigl( A''(x)B(y) + A(x)B''(y) \Bigr)\\[1mm]
   *    -\Bigl( A''(x)B(y) + A(x)B''(y) \Bigr)
   *  \end{pmatrix}
   * @f]
   *
   * @f[
   *  g(x, y) = u(x, y)
   * @f]
   */
  TEST_P(ManufacturedPoissonTest64x64, VectorPoisson_MixedTrigonometricExponential)
  {
    auto pi = Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());

    auto A = (1 - pow(F::x, 2)) * cos(3 * pi * F::x);
    auto B = (1 - pow(F::y, 2)) * sin(3 * pi * F::y);
    auto A_dd = cos(3 * pi * F::x) * (-2 - 9 * pi * pi * (1 - pow(F::x, 2)))
                + 12 * pi * F::x * sin(3 * pi * F::x);
    auto B_dd = sin(3 * pi * F::y) * (-2 - 9 * pi * pi * (1 - pow(F::y, 2)))
                - 12 * pi * F::y * cos(3 * pi * F::y);
    auto f_expr = - (A_dd * B + A * B_dd);

    VectorFunction f{ f_expr, f_expr };

    auto solution = A * B;
    VectorFunction sol{ solution, solution };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, sol);
    CG(poisson).solve();

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - sol), 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y) = \left(
   *    \begin{aligned}
   *      \sin(\pi x)\sin(\pi y)\\[1mm]
   *      \sin(2 \pi x)\sin(2 \pi y)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  f(x, y) = \left(
   *    \begin{aligned}
   *      2 \pi^2 \sin(\pi x)\sin(\pi y)\\[1mm]
   *      8 \pi^2 \sin(2 \pi x)\sin(2 \pi y)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest64x64, VectorPoisson_SimpleSine)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{
      2 * pi * pi * sin(pi * F::x) * sin(pi * F::y),
      8 * pi * pi * sin(2 * pi * F::x) * sin(2 * pi * F::y)
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();

    VectorFunction solution{
      sin(pi * F::x) * sin(pi * F::y),
      sin(2 * pi * F::x) * sin(2 * pi * F::y)
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * Let
   *
   * @f[
   *  u(x, y) = \left(
   *    \begin{aligned}
   *      \sin(\pi x^2)\sin(\pi y^2)\\[1mm]
   *      \sin(2 \pi x^2)\sin(2 \pi y^2)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  f(x, y) = \left(
   *    \begin{aligned}
   *      -\left(
   *        \left[2 \pi \cos(\pi x^2) - 4 \pi^2 x^2 \sin(\pi x^2)\right]\sin(\pi y^2)
   *        + \left[2 \pi \cos(\pi y^2) - 4 \pi^2 y^2 \sin(\pi y^2)\right]\sin(\pi x^2)
   *      \right)\\[1mm]
   *      -\left(
   *        \left[4 \pi \cos(2 \pi x^2) - 16 \pi^2 x^2 \sin(2 \pi x^2)\right]\sin(2 \pi y^2)
   *        + \left[4 \pi \cos(2 \pi y^2) - 16 \pi^2 y^2 \sin(2 \pi y^2)\right]\sin(2 \pi x^2)
   *      \right)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest64x64, VectorPoisson_VariableFrequency)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{
      -(
         (2 * pi * cos(pi * F::x * F::x) - 4 * pi * pi * F::x * F::x * sin(pi * F::x * F::x))
           * sin(pi * F::y * F::y)
         +
         (2 * pi * cos(pi * F::y * F::y) - 4 * pi * pi * F::y * F::y * sin(pi * F::y * F::y))
           * sin(pi * F::x * F::x)
       ),
      -(
         (4 * pi * cos(2 * pi * F::x * F::x) - 16 * pi * pi * F::x * F::x * sin(2 * pi * F::x * F::x))
           * sin(2 * pi * F::y * F::y)
         +
         (4 * pi * cos(2 * pi * F::y * F::y) - 16 * pi * pi * F::y * F::y * sin(2 * pi * F::y * F::y))
           * sin(2 * pi * F::x * F::x)
       )
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();

    VectorFunction solution{
      sin(pi * F::x * F::x) * sin(pi * F::y * F::y),
      sin(2 * pi * F::x * F::x) * sin(2 * pi * F::y * F::y)
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y) = \left(
   *    \begin{aligned}
   *      x (1-x) \sin(2 \pi y)\\[1mm]
   *      y (1-y) \sin(2 \pi x)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  f(x, y) = \left(
   *    \begin{aligned}
   *      2 \sin(2 \pi y) + 4 \pi^2 x (1-x) \sin(2 \pi y)\\[1mm]
   *      2 \sin(2 \pi x) + 4 \pi^2 y (1-y) \sin(2 \pi x)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, VectorPoisson_VariableAmplitude)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{
      2 * sin(2 * pi * F::y) + 4 * pi * pi * F::x * (1 - F::x) * sin(2 * pi * F::y),
      2 * sin(2 * pi * F::x) + 4 * pi * pi * F::y * (1 - F::y) * sin(2 * pi * F::x)
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();

    VectorFunction solution{
      F::x * (1 - F::x) * sin(2 * pi * F::y),
      F::y * (1 - F::y) * sin(2 * pi * F::x)
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y) = \left(
   *    \begin{aligned}
   *      x (1-x) y (1-y)\\[1mm]
   *      x (1-x) y (1-y)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  f(x, y) = \left(
   *    \begin{aligned}
   *      2 y (1-y) + 2 x (1-x)\\[1mm]
   *      2 y (1-y) + 2 x (1-x)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, VectorPoisson_Polynomial)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{
      2 * F::y * (1 - F::y) + 2 * F::x * (1 - F::x),
      2 * F::y * (1 - F::y) + 2 * F::x * (1 - F::x)
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();

    VectorFunction solution{
      F::x * (1 - F::x) * F::y * (1 - F::y),
      F::x * (1 - F::x) * F::y * (1 - F::y)
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0, 1] \times [0, 1]
   * @f]
   *
   * @f[
   *  u(x, y) = \left(
   *    \begin{aligned}
   *      x (1-x) \sin(\pi y)\\[1mm]
   *      y (1-y) \sin(\pi x)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  f(x, y) = \left(
   *    \begin{aligned}
   *      2 \sin(\pi y) + \pi^2 x (1-x) \sin(\pi y)\\[1mm]
   *      2 \sin(\pi x) + \pi^2 y (1-y) \sin(\pi x)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, VectorPoisson_TrigonometricPolynomial)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{
      2 * sin(pi * F::y) + pi * pi * F::x * (1 - F::x) * sin(pi * F::y),
      2 * sin(pi * F::x) + pi * pi * F::y * (1 - F::y) * sin(pi * F::x)
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();

    VectorFunction solution{
      F::x * (1 - F::x) * sin(pi * F::y),
      F::y * (1 - F::y) * sin(pi * F::x)
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x, y) = \left(
   *    \begin{aligned}
   *      \cos(\pi x) \cos(\pi y)\\[1mm]
   *      \sin(\pi x) \sin(\pi y)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  f(x, y) = \left(
   *    \begin{aligned}
   *      2 \pi^2 \cos(\pi x) \cos(\pi y)\\[1mm]
   *      2 \pi^2 \sin(\pi x) \sin(\pi y)
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  g(x, y) = u(x,y)
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest32x32, VectorPoisson_NonhomogeneousDirichlet)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{
      2 * pi * pi * cos(pi * F::x) * cos(pi * F::y),
      2 * pi * pi * sin(pi * F::x) * sin(pi * F::y)
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    // Use the manufactured solution as Dirichlet data.
    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, VectorFunction{
                cos(pi * F::x) * cos(pi * F::y),
                sin(pi * F::x) * sin(pi * F::y)
              });
    CG(poisson).solve();

    VectorFunction solution{
      cos(pi * F::x) * cos(pi * F::y),
      sin(pi * F::x) * sin(pi * F::y)
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x, y) = \left(
   *    \begin{aligned}
   *      \sin(\pi x) e^y\\[1mm]
   *      \sin(\pi y) e^x
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  f(x, y) = \left(
   *    \begin{aligned}
   *      (\pi^2-1)\sin(\pi x) e^y\\[1mm]
   *      (\pi^2-1)\sin(\pi y) e^x
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  g(x, y) = u(x,y)
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, VectorPoisson_MixedBoundary)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{
      (pi * pi - 1) * sin(pi * F::x) * exp(F::y),
      (pi * pi - 1) * sin(pi * F::y) * exp(F::x)
    };

    TrialFunction u(vh);
    TestFunction  v(vh);

    // Apply Dirichlet conditions on the entire boundary.
    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, VectorFunction{
                sin(pi * F::x) * exp(F::y),
                sin(pi * F::y) * exp(F::x)
              });
    CG(poisson).solve();

    VectorFunction solution{
      sin(pi * F::x) * exp(F::y),
      sin(pi * F::y) * exp(F::x)
    };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x, y) = \left(
   *    \begin{aligned}
   *      e^{x(1-x)y(1-y)} - 1\\[1mm]
   *      e^{x(1-x)y(1-y)} - 1
   *    \end{aligned}
   *  \right)
   * @f]
   *
   * @f[
   *  f(x, y) = -e^{x(1-x)y(1-y)} \left(
   *    -2y(1-y)-2x(1-x) + (1-2x)^2 \left[y(1-y)\right]^2 + (1-2y)^2 \left[x(1-x)\right]^2
   *  \right)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, VectorPoisson_Exponential)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    auto expr = F::x * (1 - F::x) * F::y * (1 - F::y);
    auto u_expr = exp(expr) - 1;
    auto f_expr = - exp(expr) * (
                    - 2 * F::y * (1 - F::y) - 2 * F::x * (1 - F::x)
                    + pow((1 - 2 * F::x), 2) * pow(F::y * (1 - F::y), 2)
                    + pow((1 - 2 * F::y), 2) * pow(F::x * (1 - F::x), 2)
                  );

    VectorFunction f{ f_expr, f_expr };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();

    VectorFunction solution{ u_expr, u_expr };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x, y) = x^2(1-x)^2 \, y^2(1-y)^2
   * @f]
   *
   * @f[
   *  f(x, y) = -\left(
   *    (2-12x+12x^2) \, y^2 (1-y)^2 + x^2 (1-x)^2 (2-12y+12y^2)
   *  \right)
   * @f]
   *
   * @f[
   *  g(x, y) = 0
   * @f]
   *
   */
  TEST_P(ManufacturedPoissonTest16x16, VectorPoisson_HigherOrder)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    auto A = pow(F::x, 2) * pow((1 - F::x), 2);
    auto B = pow(F::y, 2) * pow((1 - F::y), 2);
    auto f_expr = - ( (2 - 12 * F::x + 12 * pow(F::x, 2)) * B
                      + A * (2 - 12 * F::y + 12 * pow(F::y, 2)) );

    VectorFunction f{ f_expr, f_expr };

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();

    VectorFunction solution{ A * B, A * B };

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x,y) = x+y
   * @f]
   *
   * @f[
   *  f(x,y) = 0
   * @f]
   *
   * @f[
   *  g(x,y) = u(x,y)
   * @f]
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_LinearNonhomogeneous)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);
    auto f = Zero();
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, F::x + F::y);
    CG(poisson).solve();
    auto solution = F::x + F::y;
    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x,y) = x(1-x)+y(1-y)
   * @f]
   *
   * @f[
   *  f(x,y) = 4
   * @f]
   *
   * @f[
   *  g(x,y) = u(x,y)
   * @f]
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_QuadraticNonhomogeneous)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);
    auto f = RealFunction(4);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, F::x*(1-F::x) + F::y*(1-F::y));
    CG(poisson).solve();
    auto solution = F::x*(1-F::x) + F::y*(1-F::y);
    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x,y) = x^2(1-x)^2\, y^2(1-y)^2
   * @f]
   *
   * @f[
   *  f(x,y) = -\left( \left[2(1-2x)^2-4x(1-x)\right]\, y(1-y)^2 + x^2(1-x)^2\left[2(1-2y)^2-4y(1-y)\right] \right)
   * @f]
   *
   * @f[
   *  g(x,y) = u(x,y)
   * @f]
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_CubicPolynomial)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh);
    auto A = pow(F::x,2) * pow(1 - F::x,2);
    auto B = pow(F::y,2) * pow(1 - F::y,2);
    auto A_dd = 2 * pow(1 - 2*F::x, 2) - 4 * F::x * (1 - F::x);
    auto B_dd = 2 * pow(1 - 2*F::y, 2) - 4 * F::y * (1 - F::y);
    auto f = - (A_dd * B + A * B_dd);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, A*B);
    CG(poisson).solve();
    auto solution = A*B;
    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x,y) = \sin(\pi x)\sinh(\pi y)
   * @f]
   *
   * @f[
   *  f(x,y) = 0
   * @f]
   *
   * @f[
   *  g(x,y) = u(x,y)
   * @f]
   */
  TEST_P(ManufacturedPoissonTest64x64, Poisson_Harmonic)
  {
    auto pi = Math::Constants::pi();
    Mesh mesh = this->getMesh();
    P1 vh(mesh);
    auto f = Zero();
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, sin(pi*F::x)*sinh(pi*F::y));
    CG(poisson).solve();
    auto solution = sin(pi*F::x)*sinh(pi*F::y);
    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x,y) = \sin(2\pi x)\sin(2\pi y)
   * @f]
   *
   * @f[
   *  f(x,y) = 8 \pi^2 \sin(2\pi x)\sin(2\pi y)
   * @f]
   *
   * @f[
   *  g(x,y) = 0
   * @f]
   */
  TEST_P(ManufacturedPoissonTest64x64, Poisson_SineDouble)
  {
    auto pi = Math::Constants::pi();
    Mesh mesh = this->getMesh();
    P1 vh(mesh);
    auto f = 8 * pi * pi * sin(2*pi*F::x) * sin(2*pi*F::y);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();
    auto solution = sin(2*pi*F::x)*sin(2*pi*F::y);
    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x,y) = \sin(\pi x)\cos(\pi y)
   * @f]
   *
   * @f[
   *  f(x,y) = 2 \pi^2 \sin(\pi x)\cos(\pi y)
   * @f]
   *
   * @f[
   *  g(x,y) = u(x,y)
   * @f]
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_SineCosine_Nonhomogeneous)
  {
    auto pi = Math::Constants::pi();
    Mesh mesh = this->getMesh();
    P1 vh(mesh);
    auto f = 2 * pi * pi * sin(pi*F::x)*cos(pi*F::y);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, sin(pi*F::x)*cos(pi*F::y));
    CG(poisson).solve();
    auto solution = sin(pi*F::x)*cos(pi*F::y);
    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x,y) = \cos(\pi x)-\cos(\pi y)
   * @f]
   *
   * @f[
   *  f(x,y) = \pi^2\cos(\pi x)-\pi^2\cos(\pi y)
   * @f]
   *
   * @f[
   *  g(x,y) = u(x,y)
   * @f]
   */
  TEST_P(ManufacturedPoissonTest16x16, Poisson_Cosine_Nonhomogeneous)
  {
    auto pi = Math::Constants::pi();
    Mesh mesh = this->getMesh();
    P1 vh(mesh);
    auto f = pi*pi*cos(pi*F::x) - pi*pi*cos(pi*F::y);
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, cos(pi*F::x)-cos(pi*F::y));
    CG(poisson).solve();
    auto solution = cos(pi*F::x)-cos(pi*F::y);
    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x,y) = x^2(1-x)\, y(1-y)^2
   * @f]
   *
   * @f[
   *  f(x,y) = -\left( (2-6x)\, y(1-y)^2 + x^2(1-x)\, (-4+6y) \right)
   * @f]
   *
   * @f[
   *  g(x,y) = u(x,y)
   * @f]
   */
  TEST_P(ManufacturedPoissonTest16x16, VectorPoisson_RandomPolynomial)
  {
    Mesh mesh = this->getMesh();
    P1 vh(mesh, mesh.getSpaceDimension());
    auto A = pow(F::x,2) * (1 - F::x);
    auto B = F::y * pow(1 - F::y,2);
    auto f_expr = - ((2 - 6*F::x)*B + A * (-4 + 6*F::y));
    VectorFunction f{ f_expr, f_expr };
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, VectorFunction{A*B, A*B});
    CG(poisson).solve();
    VectorFunction solution{A*B, A*B};
    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - solution), 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x,y) = \begin{pmatrix} \sin(\pi x)\, y(1-y) \\ \sin(\pi x)\, y(1-y) \end{pmatrix}
   * @f]
   *
   * @f[
   *  f(x,y) = \begin{pmatrix} \pi^2\sin(\pi x)\, y(1-y) + 2\sin(\pi x) \\ \pi^2\sin(\pi x)\, y(1-y) + 2\sin(\pi x) \end{pmatrix}
   * @f]
   *
   * @f[
   *  g(x,y) = 0
   * @f]
   */
  TEST_P(ManufacturedPoissonTest16x16, VectorPoisson_TrigonometricMixed)
  {
    auto pi = Math::Constants::pi();
    Mesh mesh = this->getMesh();
    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{ pi*pi*sin(pi*F::x)*F::y*(1-F::y)+2*sin(pi*F::x),
                      pi*pi*sin(pi*F::x)*F::y*(1-F::y)+2*sin(pi*F::x) };
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();
    VectorFunction solution{ sin(pi*F::x)*F::y*(1-F::y),
                             sin(pi*F::x)*F::y*(1-F::y) };
    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution()-solution), 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @f[
   *  \Omega = [0,1] \times [0,1]
   * @f]
   *
   * @f[
   *  u(x,y) = \begin{pmatrix} \sin(\pi x)\sin(2\pi y) \\ \sin(2\pi x)\sin(\pi y) \end{pmatrix}
   * @f]
   *
   * @f[
   *  f(x,y) = \begin{pmatrix} 5\pi^2\sin(\pi x)\sin(2\pi y) \\ 5\pi^2\sin(2\pi x)\sin(\pi y) \end{pmatrix}
   * @f]
   *
   * @f[
   *  g(x,y) = 0
   * @f]
   */
  TEST_P(ManufacturedPoissonTest32x32, VectorPoisson_Sinusoidal)
  {
    auto pi = Math::Constants::pi();
    Mesh mesh = this->getMesh();
    P1 vh(mesh, mesh.getSpaceDimension());
    VectorFunction f{ 5 * pi * pi * sin(pi*F::x)*sin(2*pi*F::y),
                      5 * pi * pi * sin(2*pi*F::x)*sin(pi*F::y) };
    TrialFunction u(vh);
    TestFunction v(vh);
    Problem poisson(u, v);
    poisson = Integral(Jacobian(u), Jacobian(v))
            - Integral(f, v)
            + DirichletBC(u, Zero(mesh.getSpaceDimension()));
    CG(poisson).solve();
    VectorFunction solution{ sin(pi*F::x)*sin(2*pi*F::y),
                             sin(2*pi*F::x)*sin(pi*F::y) };
    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution()-solution), 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    ManufacturedPoissonTest16x16,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

  INSTANTIATE_TEST_SUITE_P(
    MeshParams32x32,
    ManufacturedPoissonTest32x32,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

  INSTANTIATE_TEST_SUITE_P(
    MeshParams64x64,
    ManufacturedPoissonTest64x64,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}

