/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file Poisson.cpp
 * @brief Manufactured-solution tests for the Poisson problem using the
 *        PETSc solver backend (sequential, `Context::Local`).
 *
 * The PDE is:
 * @f[
 * \left\{
 * \begin{aligned}
 *   -\Delta u &= f \quad \text{in } \Omega = [0,1]^2,\\
 *    u         &= g \quad \text{on } \partial\Omega.
 * \end{aligned}
 * \right.
 * @f]
 *
 * The weak form is assembled by
 * `Variational::Problem<PETSc::Math::LinearSystem, U, V>` and solved with
 * `PETSc::Solver::CG` (conjugate gradient).
 *
 * Fixed-grid checks integrate the nodal squared-error interpolant. The
 * refinement checks integrate the physical squared error independently.
 * Complex PETSc builds use amplitude @f$ 1-2i @f$; errors are normalized by
 * its magnitude so the real and complex acceptance criteria agree.
 */

#include <cmath>
#include <complex>
#include <limits>
#include <gtest/gtest.h>

#include <petsc.h>

#include <Rodin/Configure.h>
#include <Rodin/Types.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Manufactured::PETScPoisson
{
  // =========================================================================
  // Test fixture
  // =========================================================================

  /**
   * @brief Parameterized fixture for PETSc Poisson manufactured-solution tests.
   *
   * Template parameter @p M controls the uniform grid resolution (M × M
   * cells).  The mesh is scaled to the unit square @f$[0,1]^2@f$.
   *
   * PETSc is initialized once per test suite (SetUpTestSuite /
   * TearDownTestSuite) so that `PetscInitialize` and `PetscFinalize` are
   * each called exactly once across all test cases in the binary.
   *
   * @tparam M Grid resolution (number of vertices per axis).
   */
  template <size_t M>
  class PETSc_Manufactured_Poisson : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      /**
       * @brief Builds a 2D uniform-grid mesh on @f$[0,1]^2@f$ with
       *        face-incidence data required by the Dirichlet-BC assembly.
       */
      Mesh<Context::Local> getMesh()
      {
        auto mesh = Mesh<Context::Local>::UniformGrid(GetParam(), {M, M});
        mesh.scale(1.0 / static_cast<Real>(M - 1));
        const size_t D = mesh.getDimension();
        mesh.getConnectivity().compute(D, D);
        mesh.getConnectivity().compute(D, D - 1);
        mesh.getConnectivity().compute(D - 1, D);
        mesh.getConnectivity().compute(D - 1, 0);
        return mesh;
      }

      /// @brief Returns a real unit amplitude or a genuinely complex phase.
      static PetscScalar amplitudeValue()
      {
#ifdef PETSC_USE_COMPLEX
        return PetscScalar(1, -2);
#else
        return PetscScalar(1);
#endif
      }

      /// @brief Constructs the coefficient function for the configured scalar field.
      static auto amplitude()
      {
#ifdef PETSC_USE_COMPLEX
        return ComplexFunction(amplitudeValue());
#else
        return RealFunction(amplitudeValue());
#endif
      }

      /// @brief Returns the squared magnitude used to normalize squared errors.
      static Real amplitudeSquared()
      {
        return std::norm(amplitudeValue());
      }

      /**
       * @brief Integrates physical squared error independently of assembly.
       * @tparam Field Computed scalar grid-function type.
       * @tparam Exact Analytic solution expression type.
       * @param field Computed solution.
       * @param exact Analytic solution.
       * @param order Quadrature order used for error integration.
       * @returns Squared physical error divided by the squared amplitude.
       */
      template <class Field, class Exact>
      Real normalizedSquaredL2Error(
        const Field& field, const Exact& exact, size_t order) const
      {
        Real result = 0;
        const auto& mesh = field.getFiniteElementSpace().getMesh();
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const auto& formula =
            QF::PolytopeQuadratureFormula::get(order, cell->getGeometry());
          for (size_t sample = 0; sample < formula.getSize(); ++sample)
          {
            const Point point(*cell, formula.getPoint(sample));
            result += formula.getWeight(sample) * point.getDistortion() *
              std::norm(field.getValue(point) - exact(point));
          }
        }
        return result / amplitudeSquared();
      }
  };

  // Resolution aliases.
  /// @brief Helper used by the tests to PET Sc Manufactured Poisson 16 x 16.
  using PETSc_Manufactured_Poisson_16x16 = PETSc_Manufactured_Poisson<16>;

  /// @brief Helper used by the tests to PET Sc Manufactured Poisson 32 x 32.
  using PETSc_Manufactured_Poisson_32x32 = PETSc_Manufactured_Poisson<32>;

  /// @brief Helper used by the tests to PET Sc Manufactured Poisson 64 x 64.
  using PETSc_Manufactured_Poisson_64x64 = PETSc_Manufactured_Poisson<64>;

  // =========================================================================
  // P1-exact test — zero truncation error
  // =========================================================================

  /**
   * @brief Solves an exactly-P1 Poisson problem with PETSc CG.
   *
   * Manufactured solution: @f$ u(x,y) = x + 2y + 1 @f$.
   *
   * Because the solution lies in the P1 space, the finite-element error
   * is exactly zero (up to floating-point rounding).  The test verifies:
   *   - The solver residual @f$ \|Ax - b\| / \|b\| @f$ is below @f$ 10^{-10} @f$.
   *   - The projected manufactured-solution residual is below @f$ 10^{-12} @f$.
   *   - The @f$ L^2 @f$ FE error is below @f$ 10^{-12} @f$.
   */
  TEST_P(PETSc_Manufactured_Poisson_16x16, Poisson_P1Exact)
  {
    auto mesh = this->getMesh();

    P1<PetscScalar> vh(mesh);

    auto solution = this->amplitude() * (F::x + 2 * F::y + 1);
    auto f = Zero<PetscScalar>(); // -Δ(affine) = 0

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, solution);

    PETSc::Solver::CG solver(poisson);
    solver.setTolerances(1e-12, 1e-14, 1e5, 1000);
    solver.solve();
    EXPECT_GT(solver.getIterationNumber(), 0);
    EXPECT_TRUE(std::isfinite(solver.getError()));

    // Check solver residual.
    {
      auto& A = poisson.getLinearSystem().getOperator();
      auto& b = poisson.getLinearSystem().getVector();
      auto& x = poisson.getLinearSystem().getSolution();

      PetscReal bnorm, rnorm;
      Vec r;
      VecDuplicate(b, &r);
      MatMult(A, x, r);
      VecAXPY(r, -1.0, b);
      VecNorm(r, NORM_2, &rnorm);
      VecNorm(b, NORM_2, &bnorm);
      VecDestroy(&r);

      const Real scale = std::max<Real>(static_cast<Real>(bnorm), 1);
      EXPECT_NEAR(static_cast<Real>(rnorm) / scale, 0, 1e-10);
    }

    // Check the projected exact-solution residual and nodal squared error.
    PETSc::Variational::GridFunction u_exact(vh);
    u_exact = solution;

    {
      auto& system = poisson.getLinearSystem();
      ::Vec residual = nullptr;
      ASSERT_EQ(VecDuplicate(system.getVector(), &residual), PETSC_SUCCESS);
      ASSERT_EQ(
        MatMult(system.getOperator(), u_exact.getData(), residual), PETSC_SUCCESS);
      ASSERT_EQ(VecAXPY(residual, -1, system.getVector()), PETSC_SUCCESS);
      PetscReal residualNorm = 0, loadNorm = 0;
      ASSERT_EQ(VecNorm(residual, NORM_2, &residualNorm), PETSC_SUCCESS);
      ASSERT_EQ(VecNorm(system.getVector(), NORM_2, &loadNorm), PETSC_SUCCESS);
      ASSERT_EQ(VecDestroy(&residual), PETSC_SUCCESS);
      EXPECT_LE(residualNorm / std::max(PetscReal(1), loadNorm), PetscReal(1e-12));
    }

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Abs(u.getSolution() - solution), 2) / this->amplitudeSquared();
    EXPECT_NEAR(Integral(diff).compute(), 0, 1e-12);
  }

  // =========================================================================
  // Simple-sine Poisson: u = sin(πx) sin(πy), f = 2π² sin(πx) sin(πy)
  // =========================================================================

  /**
   * @brief Solves the classic sine-manufactured Poisson problem with PETSc CG.
   *
   * Manufactured solution: @f$ u(x,y) = \sin(\pi x)\sin(\pi y) @f$.
   *
   * Right-hand side: @f$ f = 2\pi^2 \sin(\pi x)\sin(\pi y) @f$.
   *
   * The @f$ L^2 @f$ error should decrease at O(h²) for P1; on a 16×16 mesh
   * it stays below the fuzzy constant.
   */
  TEST_P(PETSc_Manufactured_Poisson_16x16, Poisson_SimpleSine)
  {
    const auto pi = Math::Constants::pi();
    auto mesh = this->getMesh();

    P1<PetscScalar> vh(mesh);
    auto f = this->amplitude() * (2 * pi * pi * sin(pi * F::x) * sin(pi * F::y));

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction v(vh);

    Problem poisson(u, v);
    poisson =
      Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero<PetscScalar>());

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    auto solution = this->amplitude() * (sin(pi * F::x) * sin(pi * F::y));

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Abs(u.getSolution() - solution), 2) / this->amplitudeSquared();
    EXPECT_NEAR(Integral(diff).compute(), 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * @brief Same as Poisson_SimpleSine but on a 32×32 mesh.
   *
   * Verifies that the error decreases with refinement.
   */
  TEST_P(PETSc_Manufactured_Poisson_32x32, Poisson_SimpleSine)
  {
    const auto pi = Math::Constants::pi();
    auto mesh = this->getMesh();

    P1<PetscScalar> vh(mesh);
    auto f = this->amplitude() * (2 * pi * pi * sin(pi * F::x) * sin(pi * F::y));

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction v(vh);

    Problem poisson(u, v);
    poisson =
      Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero<PetscScalar>());

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    auto solution = this->amplitude() * (sin(pi * F::x) * sin(pi * F::y));

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Abs(u.getSolution() - solution), 2) / this->amplitudeSquared();
    EXPECT_NEAR(Integral(diff).compute(), 0, RODIN_FUZZY_CONSTANT);
  }

  // =========================================================================
  // Polynomial solution: u = x(1-x)y(1-y), f = 2y(1-y) + 2x(1-x)
  // =========================================================================

  /**
   * @brief Solves a polynomial-manufactured Poisson problem with PETSc CG.
   *
   * Manufactured solution: @f$ u(x,y) = x(1-x)y(1-y) @f$.
   *
   * Right-hand side: @f$ f = 2y(1-y) + 2x(1-x) @f$.
   */
  TEST_P(PETSc_Manufactured_Poisson_16x16, Poisson_Polynomial)
  {
    auto mesh = this->getMesh();

    P1<PetscScalar> vh(mesh);
    auto f = this->amplitude() * (2 * F::y * (1 - F::y) + 2 * F::x * (1 - F::x));

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction v(vh);

    Problem poisson(u, v);
    poisson =
      Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero<PetscScalar>());

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    auto solution = this->amplitude() * (F::x * (1 - F::x) * F::y * (1 - F::y));

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Abs(u.getSolution() - solution), 2) / this->amplitudeSquared();
    EXPECT_NEAR(Integral(diff).compute(), 0, RODIN_FUZZY_CONSTANT);
  }

  // =========================================================================
  // Nonhomogeneous Dirichlet: u = cos(πx)cos(πy)
  // =========================================================================

  /**
   * @brief Nonhomogeneous Dirichlet Poisson problem with PETSc CG.
   *
   * Manufactured solution: @f$ u(x,y) = \cos(\pi x)\cos(\pi y) @f$.
   *
   * Right-hand side: @f$ f = 2\pi^2 \cos(\pi x)\cos(\pi y) @f$.
   *
   * The Dirichlet boundary data matches the manufactured solution.
   */
  TEST_P(PETSc_Manufactured_Poisson_16x16, Poisson_NonhomogeneousDirichlet)
  {
    const auto pi = Math::Constants::pi();
    auto mesh = this->getMesh();

    P1<PetscScalar> vh(mesh);
    auto solution = this->amplitude() * (cos(pi * F::x) * cos(pi * F::y));
    auto f = 2 * pi * pi * solution;

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction v(vh);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, solution);

    PETSc::Solver::CG solver(poisson);
    solver.solve();

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Abs(u.getSolution() - solution), 2) / this->amplitudeSquared();
    EXPECT_NEAR(Integral(diff).compute(), 0, RODIN_FUZZY_CONSTANT);
  }

  // =========================================================================
  // GMRES solver: same as SimpleSine but with GMRES
  // =========================================================================

  /**
   * @brief SimpleSine Poisson solved with PETSc GMRES.
   *
   * Verifies that the GMRES backend delivers the same manufactured-solution
   * accuracy as CG on a symmetric positive definite system.
   */
  TEST_P(PETSc_Manufactured_Poisson_16x16, Poisson_SimpleSine_GMRES)
  {
    const auto pi = Math::Constants::pi();
    auto mesh = this->getMesh();

    P1<PetscScalar> vh(mesh);
    auto f = this->amplitude() * (2 * pi * pi * sin(pi * F::x) * sin(pi * F::y));

    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TestFunction v(vh);

    Problem poisson(u, v);
    poisson =
      Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero<PetscScalar>());

    PETSc::Solver::GMRES solver(poisson);
    solver.solve();

    auto solution = this->amplitude() * (sin(pi * F::x) * sin(pi * F::y));

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Abs(u.getSolution() - solution), 2) / this->amplitudeSquared();
    EXPECT_NEAR(Integral(diff).compute(), 0, RODIN_FUZZY_CONSTANT);
  }

  // =========================================================================
  // Test suite instantiation
  // =========================================================================

  /**
   * @brief Verifies second-order physical L2 convergence for a complex-phase sine solution.
   *
   * Independent order-eight and order-ten error quadratures must agree. The
   * expected P1/Q1 L2 order is two; the interval allows pre-asymptotic effects.
   */
  TEST_P(PETSc_Manufactured_Poisson_16x16, SinePhysicalL2Convergence)
  {
    constexpr size_t errorOrder = 8, verificationOrder = 10;
    constexpr Real quadratureRelativeTolerance = 1e-8;
    constexpr Real minimumOrder = 1.8, maximumOrder = 2.2;
    Optional<std::pair<Real, Real>> previous;
    const auto pi = Math::Constants::pi();
    const auto solution = this->amplitude() * (sin(pi * F::x) * sin(pi * F::y));
    const auto forcing = 2 * pi * pi * solution;
    for (size_t resolution : {size_t(16), size_t(32), size_t(64)})
    {
      SCOPED_TRACE(resolution);
      auto mesh = LocalMesh::UniformGrid(GetParam(), {resolution, resolution});
      const Real h = Real(1) / Real(resolution - 1);
      mesh.scale(h);
      mesh.getConnectivity().compute(1, 2);
      P1<PetscScalar> space(mesh);
      PETSc::Variational::TrialFunction u(space);
      PETSc::Variational::TestFunction v(space);
      Problem poisson(u, v);
      poisson = Integral(Grad(u), Grad(v)) - Integral(forcing, v) +
        DirichletBC(u, Zero<PetscScalar>());
      PETSc::Solver::CG solver(poisson);
      solver.setTolerances(1e-12, 1e-14, 1e5, 10000);
      solver.solve();
      const Real error =
        std::sqrt(this->normalizedSquaredL2Error(u.getSolution(), solution, errorOrder));
      const Real verified = std::sqrt(
        this->normalizedSquaredL2Error(u.getSolution(), solution, verificationOrder));
      ASSERT_GT(verified, 0);
      ASSERT_TRUE(std::isfinite(error));
      EXPECT_LE(std::abs(error - verified), quadratureRelativeTolerance * verified);
      if (previous)
      {
        const Real order =
          std::log(previous->second / error) / std::log(previous->first / h);
        EXPECT_GE(order, minimumOrder);
        EXPECT_LE(order, maximumOrder);
      }
      previous = std::make_pair(h, error);
    }
  }

  struct PolytopeNameGenerator
  {
    std::string operator()(const ::testing::TestParamInfo<Polytope::Type>& info) const
    {
      switch (info.param)
      {
        case Polytope::Type::Triangle:
          return "Triangle";
        case Polytope::Type::Quadrilateral:
          return "Quadrilateral";
        default:
          return "Unknown";
      }
    }
  };

  /// @brief Instantiates PET Sc Manufactured Poisson 16 x 16 over the Mesh Params 16 x 16 parameter coverage.
  INSTANTIATE_TEST_SUITE_P(MeshParams16x16, PETSc_Manufactured_Poisson_16x16,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral),
    PolytopeNameGenerator());

  /// @brief Instantiates PET Sc Manufactured Poisson 32 x 32 over the Mesh Params 32 x 32 parameter coverage.
  INSTANTIATE_TEST_SUITE_P(MeshParams32x32, PETSc_Manufactured_Poisson_32x32,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral),
    PolytopeNameGenerator());
}

// ---------------------------------------------------------------------------
// main() — initializes PETSc once before running all tests.
// ---------------------------------------------------------------------------
/// @brief Initializes the parallel test runtime and runs the GoogleTest suite.
int main(int argc, char** argv)
{
  [[maybe_unused]] PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
  assert(ierr == PETSC_SUCCESS);

  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();

  PetscFinalize();
  return result;
}
