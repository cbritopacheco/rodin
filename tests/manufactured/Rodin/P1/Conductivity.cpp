/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
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
 * @brief Manufactured solutions for the conductivity problem using P1 spaces.
 *
 * The system is given by:
 * @f[
 * \left\{
 * \begin{aligned}
 *   -\operatorname{div}( \gamma \nabla u) &= f \quad \text{in } \Omega,\\
 *  u &= g \quad \text{on } \partial\Omega.
 * \end{aligned}
 * \right.
 * @f]
 *
 *
 * The weak formulation is: Find @f$ u \in V @f$ such that
 * @f[
 *   \int_\Omega \gamma \nabla u \cdot \nabla v \,dx = \int_\Omega f \, v \,dx,
 * @f]
 * for all @f$ v \in V @f$, with the essential boundary condition
 * @f[
 *   u = g \quad \text{on } \partial\Omega.
 * @f]
 */
namespace Rodin::Tests::Manufactured::Conductivity
{
  /**
   * Test 1.
   *
   * @f[
   *  \Omega = [0,1]\times[0,1]
   * @f]
   *
   * @f[
   *  \Gamma = \partial \Omega
   * @f]
   *
   * @f[
   *  \gamma(x,y)=1+x+y
   * @f]
   *
   * @f[
   *  f(x,y)=2\pi^2(1+x+y)\sin(\pi x)\sin(\pi y)
   *    -\pi\Bigl[\sin(\pi y)\cos(\pi x)+\sin(\pi x)\cos(\pi y)\Bigr]
   * @f]
   *
   * @f[
   *  g(x,y)=0
   * @f]
   *
   * @f[
   *  u(x,y)=\sin(\pi x)\sin(\pi y)
   * @f]
   *
   */
  TEST(Rodin_Manufactured_P1, Conductivity_1)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, {16,16});
    mesh.scale(1.0/15);
    mesh.getConnectivity().compute(1,2);

    P1 vh(mesh);
    auto gamma = 1 + F::x + F::y;
    auto f =
      2 * pi * pi * gamma * sin(pi * F::x) * sin(pi * F::y)
      - pi * ( sin(pi * F::y)*cos(pi * F::x) + sin(pi * F::x)*cos(pi * F::y) )
      ;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem conductivity(u, v);
    conductivity = Integral(gamma * Grad(u), Grad(v))
                 - Integral(f, v)
                 + DirichletBC(u, Zero());
    CG(conductivity).solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * Conductivity with constant conductivity.
   *
   * @f[
   *  \gamma(x,y)=2
   * @f]
   *
   * Manufactured solution:
   * @f[
   *  u(x,y)=\sin(\pi x)\sin(\pi y)
   * @f]
   *
   * Since
   *
   * @f[
   *  \Delta u = -2\pi^2\sin(\pi x)\sin(\pi y),
   * @f]
   *
   * we have
   *
   * @f[
   *  f(x,y)=-\nabla\cdot(2\nabla u)=4\pi^2\sin(\pi x)\sin(\pi y).
   * @f]
   */
  TEST(Rodin_Manufactured_P1, Conductivity_Constant)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, {16,16});
    mesh.scale(1.0/15);
    mesh.getConnectivity().compute(1,2);

    P1 vh(mesh);
    const Real gamma = 2; // constant conductivity

    auto f = 4 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    TrialFunction u(vh);
    TestFunction v(vh);

    Problem conductivity(u, v);
    conductivity = Integral(gamma * Grad(u), Grad(v))
                 - Integral(f, v)
                 + DirichletBC(u, Zero());
    CG(conductivity).solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * Conductivity with an exponential manufactured solution.
   *
   * Let
   *
   * @f[
   *  \phi(x,y)= x(1-x)y(1-y),
   * @f]
   *
   * and choose
   *
   * @f[
   *  u(x,y)= e^{\phi(x,y)}-1,
   * @f]
   *
   * so that @f$u=0@f$ on @f$\partial\Omega@f$. With
   *
   * @f[
   *  \gamma(x,y)= 1+x^2,
   * @f]
   *
   * we compute
   *
   * @f[
   *  f(x,y)=-\nabla\cdot\Bigl(\gamma\nabla u\Bigr)
   *  = -e^{\phi}\Bigl[\gamma_x\phi_x+\gamma_y\phi_y+\gamma\Bigl(\phi_x^2+\phi_{xx}+\phi_y^2+\phi_{yy}\Bigr)\Bigr].
   * @f]
   */
  TEST(Rodin_Manufactured_P1, Conductivity_Exponential)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, {16,16});
    mesh.scale(1.0/15);
    mesh.getConnectivity().compute(1,2);

    P1 vh(mesh);

    auto phi   = F::x * (1 - F::x) * F::y * (1 - F::y);
    auto u_expr = exp(phi) - 1;
    auto gamma  = 1 + pow(F::x, 2);
    auto gamma_x = 2 * F::x;
    auto gamma_y = Zero();

    auto phi_x  = (1 - 2 * F::x) * F::y * (1 - F::y);
    auto phi_y  = (1 - 2 * F::y) * F::x * (1 - F::x);
    auto phi_xx = -2 * F::y * (1 - F::y);
    auto phi_yy = -2 * F::x * (1 - F::x);

    auto f = - exp(phi) * ( gamma_x * phi_x + gamma_y * phi_y
             + gamma * ( pow(phi_x, 2) + phi_xx + pow(phi_y, 2) + phi_yy) );

    TrialFunction u(vh);
    TestFunction v(vh);

    Problem conductivity(u, v);
    conductivity = Integral(gamma * Grad(u), Grad(v))
                 - Integral(f, v)
                 + DirichletBC(u, Zero());
    CG(conductivity).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * Conductivity with a polynomial manufactured solution.
   *
   * Let
   *
   * @f[
   *  u(x,y)=x(1-x)y(1-y),
   * @f]
   *
   * and choose
   *
   * @f[
   *  \gamma(x,y)=1+x+y.
   * @f]
   *
   * Then, since
   *
   * @f[
   *  u_x=(1-2x)y(1-y),\quad u_{xx}=-2y(1-y),
   * @f]
   *
   * @f[
   *  u_y=(1-2y)x(1-x),\quad u_{yy}=-2x(1-x),
   * @f]
   *
   * and with @f$\gamma_x=1,\;\gamma_y=1@f$, we have
   *
   * @f[
   *  f(x,y)=-\Bigl[\gamma_xu_x+\gamma\,u_{xx}+\gamma_yu_y+\gamma\,u_{yy}\Bigr].
   * @f]
   */
  TEST(Rodin_Manufactured_P1, Conductivity_Polynomial)
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, {16,16});
    mesh.scale(1.0/15);
    mesh.getConnectivity().compute(1,2);

    P1 vh(mesh);

    auto u_expr = F::x * (1 - F::x) * F::y * (1 - F::y);
    auto u_x  = (1 - 2 * F::x) * F::y * (1 - F::y);
    auto u_y  = (1 - 2 * F::y) * F::x * (1 - F::x);
    auto u_xx = -2 * F::y * (1 - F::y);
    auto u_yy = -2 * F::x * (1 - F::x);

    auto gamma = 1 + F::x + F::y;
    auto gamma_x = 1;
    auto gamma_y = 1;

    auto f = - ( gamma_x * u_x + gamma * u_xx + gamma_y * u_y + gamma * u_yy );

    TrialFunction u(vh);
    TestFunction v(vh);

    Problem conductivity(u, v);
    conductivity = Integral(gamma * Grad(u), Grad(v))
                 - Integral(f, v)
                 + DirichletBC(u, Zero());
    CG(conductivity).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * Conductivity with mixed boundary data.
   *
   * Manufactured solution:
   *
   * @f[
   *  u(x,y)=\sin(\pi x)e^y,
   * @f]
   *
   * with
   *
   * @f[
   *  \gamma(x,y)=1+x.
   * @f]
   *
   * Since
   *
   * @f[
   *  u_x=\pi\cos(\pi x)e^y,\quad u_y=\sin(\pi x)e^y,
   * @f]
   *
   * and
   *
   * @f[
   *  u_{xx}=-\pi^2\sin(\pi x)e^y,\quad u_{yy}=\sin(\pi x)e^y,
   * @f]
   *
   * with @f$\gamma_x=1,\;\gamma_y=0@f$, the forcing function is given by
   *
   * @f[
   *  f(x,y)=-\Bigl[u_x+(1+x)(u_{xx}+u_{yy})\Bigr].
   * @f]
   */
  TEST(Rodin_Manufactured_P1, Conductivity_MixedBoundary)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, {16,16});
    mesh.scale(1.0/15);
    mesh.getConnectivity().compute(1,2);

    P1 vh(mesh);

    auto u_expr = sin(pi * F::x) * exp(F::y);
    auto u_x  = pi * cos(pi * F::x) * exp(F::y);
    auto u_y  = sin(pi * F::x) * exp(F::y);
    auto u_xx = -pi * pi * sin(pi * F::x) * exp(F::y);
    auto u_yy = sin(pi * F::x) * exp(F::y);

    auto gamma = 1 + F::x; // gamma(x,y)=1+x
    // gamma_x = 1, gamma_y = 0.
    auto f = - ( u_x + (1 + F::x) * (u_xx + u_yy) );

    TrialFunction u(vh);
    TestFunction v(vh);

    Problem conductivity(u, v);
    conductivity = Integral(gamma * Grad(u), Grad(v))
                 - Integral(f, v)
                 + DirichletBC(u, sin(pi * F::x) * exp(F::y));
    CG(conductivity).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * Conductivity with nonhomogeneous Dirichlet data.
   *
   * Manufactured solution:
   *
   * @f[
   *  u(x,y)=\cos(\pi x)\cos(\pi y),
   * @f]
   *
   * with
   *
   * @f[
   *  \gamma(x,y)=1+x.
   * @f]
   *
   * Then,
   *
   * @f[
   *  f(x,y)=-\Bigl[(\gamma u_x)_x+(\gamma u_y)_y\Bigr],
   * @f]
   *
   * where
   *
   * @f[
   *  u_x=-\pi\sin(\pi x)\cos(\pi y),\quad u_{xx}=-\pi^2\cos(\pi x)\cos(\pi y),
   * @f]
   *
   * @f[
   *  u_y=-\pi\cos(\pi x)\sin(\pi y),\quad u_{yy}=-\pi^2\cos(\pi x)\cos(\pi y).
   * @f]
   *
   * @note @f$ \gamma_x=1,\; \gamma_y=0 @f$
   */
  TEST(Rodin_Manufactured_P1, Conductivity_NonhomogeneousDirichlet)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, {16,16});
    mesh.scale(1.0/15);
    mesh.getConnectivity().compute(1,2);

    P1 vh(mesh);

    auto u_expr = cos(pi * F::x) * cos(pi * F::y);
    auto u_x  = -pi * sin(pi * F::x) * cos(pi * F::y);
    auto u_y  = -pi * cos(pi * F::x) * sin(pi * F::y);
    auto u_xx = -pi * pi * cos(pi * F::x) * cos(pi * F::y);
    auto u_yy = -pi * pi * cos(pi * F::x) * cos(pi * F::y);
    auto gamma = 1 + F::x;
    // gamma_x = 1, gamma_y = 0.
    auto f = - ( 1 * u_x + gamma * (u_xx + u_yy) );

    TrialFunction u(vh);
    TestFunction v(vh);

    Problem conductivity(u, v);
    conductivity = Integral(gamma * Grad(u), Grad(v))
                 - Integral(f, v)
                 + DirichletBC(u, cos(pi * F::x) * cos(pi * F::y));
    CG(conductivity).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    diff.setWeights();
    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }
}
