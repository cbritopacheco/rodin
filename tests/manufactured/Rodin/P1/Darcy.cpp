/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

/**
 *
 * @brief Manufactured solution for the Darcy problem using P1 spaces.
 *
 * Let the pressure and flux satisfy
 * @f[
 *   \mathbf{u} = -\nabla p \quad \text{in } \Omega,
 * @f]
 * @f[
 *   \operatorname{div}\,\mathbf{u} = f \quad \text{in } \Omega,
 * @f]
 * with the Dirichlet boundary condition
 * @f[
 *   p = p_{\mathrm{exact}} \quad \text{on } \partial\Omega.
 * @f]
 *
 * The weak formulation is: Find @f$(\mathbf{u}, p) \in \mathbb{P}_1^2 \times \mathbb{P}_1@f$ such that
 * @f[
 *   (\mathbf{u}, \mathbf{v}) - (p, \operatorname{div}\,\mathbf{v}) -
 *   (\operatorname{div}\,\mathbf{u}, q) + (f, q) = 0,
 * @f]
 * for all test functions @f$(\mathbf{v},q)@f$ with the essential boundary condition
 * @f[
 *   p = \sin(\pi x)\sin(\pi y) \quad \text{on } \partial\Omega.
 * @f]
 *
 * @note Although mixed formulations typically use an H(div)–conforming space
 * for the flux, here a vector–valued P1 space.
 */
namespace Rodin::Tests::Manufactured::Darcy
{
  template <size_t M>
  class ManufacturedDarcyTest : public ::testing::TestWithParam<Polytope::Type>
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

  using ManufacturedDarcyTest16x16 = ManufacturedDarcyTest<16>;
  using ManufacturedDarcyTest32x32 = ManufacturedDarcyTest<32>;
  using ManufacturedDarcyTest64x64 = ManufacturedDarcyTest<64>;

  /**
   * @f[
   *   \Omega = [0,1]\times[0,1]
   * @f]
   *
   * @f[
   *   p(x,y)=\sin(\pi x)\sin(\pi y)
   * @f]
   *
   * @f[
   *   \mathbf{u}(x,y)=-\nabla p(x,y)
   *   = -\pi \begin{pmatrix}
   *       \cos(\pi x)\sin(\pi y) \\
   *       \sin(\pi x)\cos(\pi y)
   *     \end{pmatrix}
   * @f]
   *
   *
   * @f[
   *   f(x,y)=2\pi^2\sin(\pi x)\sin(\pi y)
   * @f]
   *
   */
  TEST_P(ManufacturedDarcyTest16x16, Darcy_SimpleSine_P1)
  {
    auto pi = Math::Constants::pi();
    Mesh mesh = this->getMesh();

    // Flux space (vector-valued)
    P1 uh(mesh, mesh.getSpaceDimension());

    // Pressure space (scalar)
    P1 ph(mesh);

    // Manufactured solution:
    auto p_exact = sin(pi * F::x) * sin(pi * F::y);
    auto u_exact = -pi * VectorFunction{
                      cos(pi * F::x) * sin(pi * F::y),
                      sin(pi * F::x) * cos(pi * F::y)
                    };
    auto f = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    // Define trial and test functions.
    TrialFunction u(uh); // Flux trial function.
    TrialFunction p(ph); // Pressure trial function.
    TestFunction  v(uh); // Flux test function.
    TestFunction  q(ph); // Pressure test function.

    // Assemble the weak form:
    Problem darcy(u, p, v, q);
    darcy = Integral(u, v)
          - Integral(p, Div(v))
          - Integral(Div(u), q)
          + Integral(f, q)
          + DirichletBC(p, p_exact);

    // Solve the system.
    CG cg(darcy);
    cg.setTolerance(1e-12);
    cg.setMaxIterations(1000);
    cg.solve();
    std::cout << cg.success() << std::endl;

    // Compute the L^2 error for pressure.
    GridFunction diff_p(ph);
    p.getSolution().save("p.gf");
    mesh.save("mesh.mesh");
    std::exit(1);
    diff_p = Pow(p.getSolution() - p_exact, 2);
    diff_p.setWeights();
    Real error_p = Integral(diff_p).compute();

    // Compute the L^2 error for flux.
    GridFunction diff_u(ph);
    diff_u = Pow(Frobenius(u.getSolution() - u_exact), 2);
    diff_u.setWeights();
    Real error_u = Integral(diff_u).compute();

    EXPECT_NEAR(error_p, 0, RODIN_FUZZY_CONSTANT);
    EXPECT_NEAR(error_u, 0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    ManufacturedDarcyTest16x16,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}

