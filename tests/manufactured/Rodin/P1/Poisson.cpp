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
 * Contains manufactured solution tests cases for the Poisson equation.
 *
 * @f[
 * \left\{
 *   \begin{aligned}
 *     - \Delta u &= f && \mathrm{in} \ \Omega, \\
 *     u &= g && \mathrm{on} \ \Gamma.
 *   \end{aligned}
 *  \right.
 * @f]
 */
namespace Rodin::Tests::Manufactured::Poisson
{
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
  TEST(Rodin_Manufactured_P1, Poisson_1)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 16, 16 });
    mesh.scale(1.0 / 15);
    mesh.getConnectivity().compute(1, 2);

    P1 vh(mesh);
    auto f = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    // Define problem
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    // Solve
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
   * u(x,y) = \sin(\omega_1 \pi x) \sin(\omega_2 \pi y)
   * @f]
   *
   * @f[
   * \omega_1, \omega_2 \in \mathbb{N}
   * @f]
   *
   */
  TEST(Rodin_Manufactured_P1, Poisson_2)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 32, 32 });
    mesh.scale(1.0 / 31);
    mesh.getConnectivity().compute(1, 2);

    const Real omega1 = RandomInteger(1, 8)();
    const Real omega2 = RandomInteger(1, 8)();

    P1 vh(mesh);
    auto f = (omega1 * pi * omega1 * pi + omega2 * pi * omega2 * pi) * sin(omega1 * pi * F::x) * sin(omega2 * pi * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    // Define problem
    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    // Solve
    CG(poisson).solve();

    auto solution = sin(omega1 * pi * F::x) * sin(omega2 * pi * F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();

    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }
}
