#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

/**
 * @brief Manufactured solutions for a coupled system of two reaction–diffusion equations.
 *
 * The system is given by
 * @f[
 * \begin{aligned}
 *   -\Delta u + \alpha\,u + \beta\,v &= f_1(x,y) \quad \text{in } \Omega,\\[1mm]
 *   -\Delta v + \gamma\,u + \delta\,v &= f_2(x,y) \quad \text{in } \Omega,
 * \end{aligned}
 * @f]
 * with Dirichlet boundary conditions
 * @f[
 *   u(x,y)=g_1(x,y),\quad v(x,y)=g_2(x,y) \quad \text{on } \partial\Omega.
 * @f]
 *
 * The weak formulation is: Find @f$\mathbf{U}=(u,v)\in [P1]^2@f$ such that
 * @f[
 *   \int_\Omega \nabla u\cdot\nabla \phi + \nabla v\cdot\nabla\psi \,dx
 *   + \int_\Omega \Bigl[\alpha\,u\,\phi + \beta\,v\,\phi + \gamma\,u\,\psi + \delta\,v\,\psi\Bigr]dx
 *   = \int_\Omega f_1\,\phi+f_2\,\psi\,dx,
 * @f]
 * for all test functions @f$\mathbf{\Phi}=(\phi,\psi)\in [P1]^2@f$ and with
 * the essential condition @f$\mathbf{U}=(g_1,g_2)@f$ on @f$\partial\Omega@f$.
 */
namespace Rodin::Tests::Manufactured::ReactionDiffusion
{
  template <size_t M>
  class ManufacturedReactionDiffusionTest : public ::testing::TestWithParam<Polytope::Type>
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

  using ManufacturedReactionDiffusionTest16x16 = ManufacturedReactionDiffusionTest<16>;
  using ManufacturedReactionDiffusionTest32x32 = ManufacturedReactionDiffusionTest<32>;
  using ManufacturedReactionDiffusionTest64x64 = ManufacturedReactionDiffusionTest<64>;
  using ManufacturedReactionDiffusionTest128x128 = ManufacturedReactionDiffusionTest<128>;

  /**
   * @brief Sine/Cosine manufactured solution.
   *
   * @f[
   *   u(x,y)=\sin(\pi x)\sin(\pi y)
   * @f]
   *
   * Forcing term:
   * @f[
   *   f(x,y)=(2\pi^2+1)\sin(\pi x)\sin(\pi y)
   * @f]
   *
   * Weak formulation: Find \(u\in H^1_0(\Omega)\) such that
   * @f[
   *   \int_\Omega \nabla u\cdot\nabla v\,dx + \int_\Omega u\,v\,dx = \int_\Omega f\,v\,dx,
   * @f]
   * for all \(v\in H^1_0(\Omega)\).
   */
  TEST_P(ManufacturedReactionDiffusionTest16x16, ReactionDiffusion_Sine)
  {
      auto pi = Math::Constants::pi();
      Mesh mesh = this->getMesh();

      // Define a scalar P1 finite element space on the mesh.
      P1 vh(mesh);

      // Manufactured solution: u(x,y)= sin(pi*x)*sin(pi*y).
      // Forcing term: f(x,y) = (2*pi^2+1)* sin(pi*x)*sin(pi*y).
      auto f = (2 * pi * pi + 1) * sin(pi * F::x) * sin(pi * F::y);

      // Define trial and test functions.
      TrialFunction u(vh);
      TestFunction  v(vh);

      // Assemble the variational problem:
      // Find u such that ∫ (∇u·∇v + u*v) dx = ∫ f*v dx, with u=0 on ∂Ω.
      Problem rd(u, v);
      rd = Integral(Grad(u), Grad(v))
         + Integral(u, v)
         - Integral(f, v)
         + DirichletBC(u, Zero());

      // Solve the system using the conjugate–gradient (CG) method.
      CG(rd).solve();

      // The manufactured solution is u_exact = sin(pi*x)*sin(pi*y).
      auto u_exact = sin(pi * F::x) * sin(pi * F::y);

      // Compute the L² error.
      GridFunction diff(vh);
      diff = Pow(u.getSolution() - u_exact, 2);
      diff.setWeights();
      Real error = Integral(diff).compute();

      EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  TEST_P(ManufacturedReactionDiffusionTest64x64, ReactionDiffusion_VariableFrequency)
  {
    auto pi = Math::Constants::pi();
    Mesh mesh = this->getMesh();

    P1 vh(mesh);

    // Reaction coefficients
    const Real alpha = 0.5;
    const Real beta  = 1.0;
    const Real gamma = 1.5;
    const Real delta = 2.0;

    // Manufactured solution with variable frequencies:
    auto u0_exact = sin(3 * pi * F::x) * sin(2 * pi * F::y);
    auto u1_exact = sin(2 * pi * F::x) * sin(3 * pi * F::y);

    // Forcing terms computed to enforce the manufactured solution:
    auto f0 = (13 * pi * pi + alpha) * sin(3 * pi * F::x) * sin(2 * pi * F::y)
              + beta * sin(2 * pi * F::x) * sin(3 * pi * F::y);
    auto f1 = gamma * sin(3 * pi * F::x) * sin(2 * pi * F::y)
              + (13 * pi * pi + delta) * sin(2 * pi * F::x) * sin(3 * pi * F::y);

    // Define trial and test functions in the scalar P1 space for each variable.
    TrialFunction u0(vh), u1(vh);
    TestFunction  v0(vh), v1(vh);

    Problem rd(u0, u1, v0, v1);
    rd = Integral(Grad(u0), Grad(v0))
       + Integral(Grad(u1), Grad(v1))
       + Integral(alpha * u0, v0)
       + Integral(beta  * u1, v0)
       + Integral(gamma * u0, v1)
       + Integral(delta * u1, v1)
       - Integral(f0, v0)
       - Integral(f1, v1)
       + DirichletBC(u0, u0_exact)
       + DirichletBC(u1, u1_exact);

    CG(rd).solve();

    GridFunction diff0(vh);
    diff0 = Pow(u0.getSolution() - u0_exact, 2);
    diff0.setWeights();
    Real error0 = Integral(diff0).compute();
    EXPECT_NEAR(error0, 0, RODIN_FUZZY_CONSTANT);

    GridFunction diff1(vh);
    diff1 = Pow(u1.getSolution() - u1_exact, 2);
    diff1.setWeights();
    Real error1 = Integral(diff1).compute();
    EXPECT_NEAR(error1, 0, RODIN_FUZZY_CONSTANT);
  }

  TEST_P(ManufacturedReactionDiffusionTest16x16, WeaklyCoupledReactionDiffusion)
  {
      auto pi = Math::Constants::pi();
      Mesh mesh = this->getMesh();

      // Reaction coupling parameter (weak coupling).
      const Real alpha = 0.1;

      // Define a scalar P1 finite element space for both variables.
      P1 vh(mesh);

      // Manufactured solution:
      auto u_exact = sin(pi * F::x) * sin(pi * F::y);
      auto v_exact = cos(pi * F::x) * cos(pi * F::y);

      // Forcing terms:
      auto f1 = (2 * pi * pi + 1) * sin(pi * F::x) * sin(pi * F::y)
                + alpha * cos(pi * F::x) * cos(pi * F::y);
      auto f2 = (2 * pi * pi + 1) * cos(pi * F::x) * cos(pi * F::y)
                + alpha * sin(pi * F::x) * sin(pi * F::y);

      // Define trial and test functions.
      TrialFunction u(vh), v(vh);
      TestFunction  phi(vh), psi(vh);

      Problem rd(u, v, phi, psi);
      rd = Integral(Grad(u), Grad(phi))
         + Integral(u, phi)
         + Integral(alpha * v, phi)
         - Integral(f1, phi)
         + Integral(Grad(v), Grad(psi))
         + Integral(v, psi)
         + Integral(alpha * u, psi)
         - Integral(f2, psi)
         + DirichletBC(u, Zero())
         + DirichletBC(v, v_exact);

      // Solve the coupled system.
      CG(rd).solve();

      // Compute the L² error for u.
      GridFunction diff_u(vh);
      diff_u = Pow(u.getSolution() - u_exact, 2);
      diff_u.setWeights();
      Real error_u = Integral(diff_u).compute();

      // Compute the L² error for v.
      GridFunction diff_v(vh);
      diff_v = Pow(v.getSolution() - v_exact, 2);
      diff_v.setWeights();
      Real error_v = Integral(diff_v).compute();

      EXPECT_NEAR(error_u, 0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR(error_v, 0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
      MeshParams16x16,
      ManufacturedReactionDiffusionTest16x16,
      ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

  // INSTANTIATE_TEST_SUITE_P(
  //   MeshParams32x32,
  //   ManufacturedReactionDiffusionTest32x32,
  //   ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  // );

  INSTANTIATE_TEST_SUITE_P(
    MeshParams64x64,
    ManufacturedReactionDiffusionTest64x64,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}

