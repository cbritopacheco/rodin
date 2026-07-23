/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file KSPFactorizationReuseTest.cpp
 * @brief Factorization reuse contract of Rodin::Solver::KSP.
 *
 * A KSP that is reused across solves must not repeat work the operator did
 * not invalidate:
 *
 * - Same operator, same values: no symbolic and no numeric factorization.
 * - Same operator object, refilled values, same sparsity: symbolic reused,
 *   numeric redone.
 * - A different operator object: everything redone.
 *
 * The observable is PETSc's own event log: MatLUFactorSym and MatLUFactorNum
 * carry a call count, so the tests assert on how many times each ran rather
 * than on wall-clock time. Logging is enabled explicitly in main(); if PETSc
 * was configured without it, the tests skip rather than assert vacuously.
 *
 * Every reuse assertion is paired with a numerical assertion: skipping work is
 * only correct if the answer is still right. The systems here are mass
 * matrices, @f$ \alpha \int u v = \int c v @f$, whose exact discrete solution
 * is the constant @f$ c / \alpha @f$ at every P1 degree of freedom.
 */

#include <gtest/gtest.h>

#include <petsc.h>
#include <petscksp.h>
#include <petsclog.h>
#include <petscmat.h>

#include <Rodin/Geometry.h>
#include <Rodin/PETSc.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit::KSPFactorizationReuse
{
  namespace
  {
    /// @brief Number of symbolic and numeric LU factorizations logged so far.
    struct FactorCounts
    {
      int symbolic = 0;
      int numeric = 0;
    };

    /// @brief Reads the cumulative LU factorization event counts from PETSc.
    FactorCounts factorCounts()
    {
      const auto count = [](const char* name) -> int
      {
        PetscLogEvent event = -1;
        PetscErrorCode ierr = PetscLogEventGetId(name, &event);
        assert(ierr == PETSC_SUCCESS);
        if (event < 0)
          return 0;

        PetscEventPerfInfo info;
        ierr = PetscLogEventGetPerfInfo(PETSC_DETERMINE, event, &info);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
        return info.count;
      };

      FactorCounts counts;
      counts.symbolic = count("MatLUFactorSym");
      counts.numeric = count("MatLUFactorNum");
      return counts;
    }

    /// @brief Whether PETSc is accumulating the event counts the tests read.
    bool loggingActive()
    {
      PetscBool active = PETSC_FALSE;
      PetscErrorCode ierr = PetscLogIsActive(&active);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
      return active == PETSC_TRUE;
    }

    /// @brief Largest deviation of a grid function from the constant @p expected.
    template <class GridFunctionType>
    PetscReal maxDeviation(const GridFunctionType& gf, PetscReal expected)
    {
      ::Vec vec = gf.getData();
      PetscInt n = 0;
      PetscErrorCode ierr = VecGetLocalSize(vec, &n);
      assert(ierr == PETSC_SUCCESS);

      const PetscScalar* values = nullptr;
      ierr = VecGetArrayRead(vec, &values);
      assert(ierr == PETSC_SUCCESS);

      PetscReal deviation = 0;
      for (PetscInt i = 0; i < n; ++i)
        deviation = std::max(deviation, std::abs(PetscRealPart(values[i]) - expected));

      ierr = VecRestoreArrayRead(vec, &values);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
      return deviation;
    }

    /// @brief Forces every KSP in the test onto a direct LU solve.
    ///
    /// The reuse contract is only observable through a preconditioner that
    /// actually factorizes. PETSc's built-in LU is used so the test does not
    /// require MUMPS.
    void setDirectSolverOptions(const char* prefix)
    {
      const std::string kspType = std::string("-") + prefix + "ksp_type";
      const std::string pcType = std::string("-") + prefix + "pc_type";
      PetscErrorCode ierr =
        PetscOptionsSetValue(PETSC_NULLPTR, kspType.c_str(), "preonly");
      assert(ierr == PETSC_SUCCESS);
      ierr = PetscOptionsSetValue(PETSC_NULLPTR, pcType.c_str(), "lu");
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
    }

    constexpr PetscReal tolerance = 1e-10;
  }

  /**
   * A KSP re-solving an operator whose values never changed must not
   * factorize again, and must keep returning the right answer.
   */
  TEST(Rodin_PETSc_KSPFactorizationReuse, UnchangedOperatorSkipsFactorization)
  {
    if (!loggingActive())
      GTEST_SKIP() << "PETSc event logging is unavailable.";

    setDirectSolverOptions("reuse_unchanged_");

    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {8, 8});
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);

    constexpr PetscReal alpha = 3.0;
    constexpr PetscReal c = 6.0;

    Problem mass(u, v);
    mass = Integral(RealFunction(alpha) * u, v) - Integral(RealFunction(c), v);

    Rodin::Solver::KSP ksp(mass);
    ksp.setPrefix(std::string("reuse_unchanged_"));

    const auto before = factorCounts();
    mass.solve(ksp);
    const auto first = factorCounts();

    ASSERT_GT(first.symbolic, before.symbolic)
      << "the first solve did not factorize; the PC is not a direct solver.";
    EXPECT_NEAR(maxDeviation(u.getSolution(), c / alpha), 0.0, tolerance);

    // Re-solving the very same assembled system must be pure back-substitution.
    for (int i = 0; i < 3; ++i)
    {
      mass.solve(ksp);
      const auto again = factorCounts();
      EXPECT_EQ(again.symbolic, first.symbolic)
        << "symbolic factorization repeated on solve " << i;
      EXPECT_EQ(again.numeric, first.numeric)
        << "numeric factorization repeated on solve " << i;
      EXPECT_NEAR(maxDeviation(u.getSolution(), c / alpha), 0.0, tolerance);
    }
  }

  /**
   * Refilling the same operator object with new values must reuse the symbolic
   * factorization, redo the numeric one, and produce the new solution. A stale
   * factorization would return the old answer and is the failure this guards.
   */
  TEST(Rodin_PETSc_KSPFactorizationReuse, RefilledOperatorReusesSymbolicOnly)
  {
    if (!loggingActive())
      GTEST_SKIP() << "PETSc event logging is unavailable.";

    setDirectSolverOptions("reuse_refilled_");

    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {8, 8});
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);

    constexpr PetscReal c = 6.0;

    Problem mass(u, v);
    Rodin::Solver::KSP ksp(mass);
    ksp.setPrefix(std::string("reuse_refilled_"));

    mass = Integral(RealFunction(2.0) * u, v) - Integral(RealFunction(c), v);
    mass.solve(ksp);
    EXPECT_NEAR(maxDeviation(u.getSolution(), c / 2.0), 0.0, tolerance);

    const auto first = factorCounts();

    // Same sparsity pattern, different values: the ordering is still valid.
    mass = Integral(RealFunction(4.0) * u, v) - Integral(RealFunction(c), v);
    mass.solve(ksp);

    const auto second = factorCounts();
    EXPECT_EQ(second.symbolic, first.symbolic)
      << "symbolic factorization was not reused";
    EXPECT_EQ(second.numeric, first.numeric + 1)
      << "numeric factorization was skipped or repeated";
    EXPECT_NEAR(maxDeviation(u.getSolution(), c / 4.0), 0.0, tolerance)
      << "solve used a stale factorization";
  }

  /**
   * The wall-shear pattern: one operator, several right-hand sides. Assembling
   * only the RHS must leave the factorization untouched and still solve each
   * system correctly.
   */
  TEST(Rodin_PETSc_KSPFactorizationReuse, MultipleRightHandSidesFactorizeOnce)
  {
    if (!loggingActive())
      GTEST_SKIP() << "PETSc event logging is unavailable.";

    setDirectSolverOptions("reuse_multirhs_");

    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {8, 8});
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);

    constexpr PetscReal alpha = 5.0;

    Problem mass(u, v);
    Rodin::Solver::KSP ksp(mass);
    ksp.setPrefix(std::string("reuse_multirhs_"));

    mass = Integral(RealFunction(alpha) * u, v) - Integral(RealFunction(1.0), v);
    mass.solve(ksp);
    EXPECT_NEAR(maxDeviation(u.getSolution(), 1.0 / alpha), 0.0, tolerance);

    const auto first = factorCounts();

    for (const PetscReal c : {2.0, 3.0, 4.0})
    {
      mass = Integral(RealFunction(alpha) * u, v) - Integral(RealFunction(c), v);
      mass.assemble(AssemblyTarget::RHS);
      mass.solve(ksp);

      const auto again = factorCounts();
      EXPECT_EQ(again.symbolic, first.symbolic) << "symbolic repeated for rhs c = " << c;
      EXPECT_EQ(again.numeric, first.numeric) << "numeric repeated for rhs c = " << c;
      EXPECT_NEAR(maxDeviation(u.getSolution(), c / alpha), 0.0, tolerance);
    }
  }

  /**
   * Handing a KSP a genuinely different operator must re-run the whole setup
   * and solve the new system, not the old one.
   */
  TEST(Rodin_PETSc_KSPFactorizationReuse, DifferentOperatorRefactorizes)
  {
    if (!loggingActive())
      GTEST_SKIP() << "PETSc event logging is unavailable.";

    setDirectSolverOptions("reuse_switch_");

    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {8, 8});
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    PETSc::Variational::TrialFunction w(fes);
    PETSc::Variational::TestFunction z(fes);

    Problem first(u, v);
    first = Integral(RealFunction(2.0) * u, v) - Integral(RealFunction(1.0), v);

    Problem second(w, z);
    second = Integral(RealFunction(8.0) * w, z) - Integral(RealFunction(1.0), z);
    second.assemble();

    Rodin::Solver::KSP ksp(first);
    ksp.setPrefix(std::string("reuse_switch_"));

    first.solve(ksp);
    EXPECT_NEAR(maxDeviation(u.getSolution(), 1.0 / 2.0), 0.0, tolerance);
    const auto firstCounts = factorCounts();

    // Same KSP, different LinearSystem: the cached operator must be dropped.
    ksp.solve(second.getLinearSystem());
    w.getSolution().setData(second.getLinearSystem().getSolution());

    const auto secondCounts = factorCounts();
    EXPECT_GT(secondCounts.symbolic, firstCounts.symbolic)
      << "operator switch reused a factorization";
    EXPECT_NEAR(maxDeviation(w.getSolution(), 1.0 / 8.0), 0.0, tolerance)
      << "operator switch solved the wrong system";
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  // The reuse assertions read PETSc event counters, which only accumulate
  // while a log handler is running.
  PetscErrorCode ierr = PetscLogDefaultBegin();
  assert(ierr == PETSC_SUCCESS);
  (void)ierr;
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
