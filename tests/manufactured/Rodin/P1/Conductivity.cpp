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
 * Contains manufactured solution tests cases for the conductivity equation.
 *
 * @f[
 * \left\{
 *   \begin{aligned}
 *     - \nabla \cdot (\gamma \nabla u) &= f && \mathrm{in} \ \Omega, \\
 *     u &= g && \mathrm{on} \ \Gamma.
 *   \end{aligned}
 *  \right.
 * @f]
 */
namespace Rodin::Tests::Manufactured::Conductivity
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
   *  \gamma(x,y) = 1 + x + y
   * @f]
   *
   * @f[
   *  f(x, y) = 2\pi^2(1+x+y)\sin(\pi x)\sin(\pi y) -\pi\Bigl[\sin(\pi
   *  y)\cos(\pi x)+\sin(\pi x)\cos(\pi y)\Bigr]
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
  TEST(Rodin_Manufactured_P1, Conductivity_1)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 16, 16 });
    mesh.scale(1.0 / 15);
    mesh.getConnectivity().compute(1, 2);

    P1 vh(mesh);
    auto gamma = 1 + F::x + F::y;
    auto f =
      2 * pi * pi * gamma * sin(pi * F::x) * sin(pi * F::y)
      - pi * (sin(pi * F::y) * cos(pi * F::x) + sin(pi * F::x) * cos(pi * F::y))
      ;

    TrialFunction u(vh);
    TestFunction  v(vh);

    // Define problem
    Problem conductivity(u, v);
    conductivity = Integral(gamma * Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

    // Solve
    CG(conductivity).solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);
    diff.setWeights();

    Real error = Integral(diff).compute();

    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }
}
