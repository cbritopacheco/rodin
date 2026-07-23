/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file MPIKSPFactorizationReuseTest.cpp
 * @brief Distributed counterpart of KSPFactorizationReuseTest.
 *
 * The sequential test pins the reuse contract for @c seqaij matrices factored
 * by PETSc's built-in LU. This one pins it for @c mpiaij matrices factored by
 * MUMPS, which is the configuration production runs use and the one where a
 * spurious nonzero-state bump is expensive: MUMPS redoes its analysis whenever
 * PETSc reports a changed sparsity pattern.
 *
 * Regression guard for the assembly fix in @ref Rodin::PETSc::Assembly::MPI:
 * @c MatZeroRows / @c MatZeroRowsColumns must not run with a globally empty
 * index set, because they bump the matrix nonzero-state counter even for zero
 * rows, which makes PETSc report DIFFERENT_NONZERO_PATTERN and forces MUMPS to
 * repeat its symbolic factorization on every reassembly of a reused system.
 *
 * The observable is PETSc's own event log (MatLUFactorSym / MatLUFactorNum
 * call counts), with every reuse assertion paired with a numerical one. The
 * MUMPS solver is driven with the same @c -mat_mumps_icntl_20/21 0 options the
 * production coronary example uses; without centralized right-hand sides this
 * MUMPS build faults inside its distributed RHS scatter.
 *
 * @see KSPFactorizationReuseTest.cpp for the contract these tests share.
 */

#include <gtest/gtest.h>

#include <petsc.h>
#include <petscksp.h>
#include <petsclog.h>
#include <petscmat.h>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/PETSc.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static boost::mpi::environment* g_env = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace Rodin::Tests::Unit::MPIKSPFactorizationReuse
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
      const auto count = [](const char* name) -> int {
        PetscLogEvent event = -1;
        PetscErrorCode ierr = PetscLogEventGetId(name, &event);
        assert(ierr == PETSC_SUCCESS);
        if (event < 0)
          return 0;

        // Stage 0 (the main stage) explicitly: PETSC_DETERMINE as a stage
        // argument is not portable to PETSc 3.19, and these tests never push a
        // log stage, so every event is recorded in the main stage.
        PetscEventPerfInfo info;
        ierr = PetscLogEventGetPerfInfo(0, event, &info);
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

    /// @brief Whether this PETSc build can factor a distributed matrix.
    bool mumpsAvailable()
    {
      PetscBool has = PETSC_FALSE;
      PetscErrorCode ierr = PetscHasExternalPackage("mumps", &has);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
      return has == PETSC_TRUE;
    }

    /// @brief Distributes a small uniform triangle mesh from rank 0.
    ///
    /// Mirrors the manufactured MPI tests: connectivity is computed and
    /// reconciled on the gathered distributed mesh so a P1 assembly can iterate
    /// it and the solution can be scattered back.
    Mesh<Context::MPI> distributeFromRoot(const Context::MPI& ctx)
    {
      const auto& comm = ctx.getCommunicator();
      Sharder<Context::MPI> sharder(ctx);

      if (comm.rank() == 0)
      {
        auto localMesh =
          Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {6, 6});
        const size_t D = localMesh.getDimension();
        localMesh.getConnectivity().compute(D, D);
        localMesh.getConnectivity().compute(D, 0);
        localMesh.getConnectivity().compute(D, D - 1);
        localMesh.getConnectivity().compute(D - 1, D);
        localMesh.getConnectivity().compute(D - 1, 0);

        BalancedCompactPartitioner partitioner(localMesh);
        partitioner.partition(static_cast<size_t>(comm.size()));
        sharder.shard(partitioner);
        sharder.scatter(0);
      }

      auto mesh = sharder.gather(0);

      const size_t D = mesh.getDimension();
      mesh.getConnectivity().compute(D, D);
      mesh.getConnectivity().compute(D, 0);
      mesh.getConnectivity().compute(D, D - 1);
      mesh.getConnectivity().compute(D - 1, D);
      mesh.getConnectivity().compute(D - 1, 0);
      mesh.getConnectivity().compute(1, 2);
      mesh.reconcile(1);
      return mesh;
    }

    /// @brief Largest local deviation of a grid function from a constant.
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

    constexpr PetscReal tolerance = 1e-9;

    bool unsupported()
    {
      return !loggingActive() || !mumpsAvailable();
    }
  }

  /**
   * Re-solving an untouched distributed operator must not factorize again.
   */
  TEST(PETSc_MPI_KSPFactorizationReuse, UnchangedOperatorSkipsFactorization)
  {
    if (unsupported())
      GTEST_SKIP() << "PETSc event logging or MUMPS is unavailable.";

    Context::MPI ctx(*g_env, *g_world);
    auto mesh = distributeFromRoot(ctx);

    P1<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);

    constexpr PetscReal alpha = 3.0;
    constexpr PetscReal c = 6.0;

    Problem mass(u, v);
    mass = Integral(RealFunction(alpha) * u, v) - Integral(RealFunction(c), v);

    Rodin::Solver::KSP ksp(mass);

    const auto before = factorCounts();
    mass.solve(ksp);
    const auto first = factorCounts();

    ASSERT_GT(first.symbolic, before.symbolic)
      << "the first solve did not factorize; the PC is not a direct solver.";
    EXPECT_NEAR(maxDeviation(u.getSolution(), c / alpha), 0.0, tolerance);

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
   * Refilling a distributed operator that carries no boundary conditions must
   * reuse MUMPS' ordering and redo only the numeric factorization. This is the
   * direct regression for the zero-row MatZeroRows bump: before the fix, the
   * empty constraint pass re-flagged the pattern and the symbolic count grew on
   * every reassembly.
   */
  TEST(PETSc_MPI_KSPFactorizationReuse, RefilledOperatorReusesSymbolicOnly)
  {
    if (unsupported())
      GTEST_SKIP() << "PETSc event logging or MUMPS is unavailable.";

    Context::MPI ctx(*g_env, *g_world);
    auto mesh = distributeFromRoot(ctx);

    P1<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);

    constexpr PetscReal c = 6.0;

    Problem mass(u, v);
    Rodin::Solver::KSP ksp(mass);

    mass = Integral(RealFunction(2.0) * u, v) - Integral(RealFunction(c), v);
    mass.solve(ksp);
    EXPECT_NEAR(maxDeviation(u.getSolution(), c / 2.0), 0.0, tolerance);

    const auto first = factorCounts();

    mass = Integral(RealFunction(4.0) * u, v) - Integral(RealFunction(c), v);
    mass.solve(ksp);

    const auto second = factorCounts();
    EXPECT_EQ(second.symbolic, first.symbolic)
      << "MUMPS reordering was repeated for an unchanged sparsity pattern";
    EXPECT_EQ(second.numeric, first.numeric + 1)
      << "numeric factorization was skipped or repeated";
    EXPECT_NEAR(maxDeviation(u.getSolution(), c / 4.0), 0.0, tolerance)
      << "solve used a stale factorization";
  }

  /**
   * The wall-shear pattern, distributed: one operator, several right-hand
   * sides, assembled with AssemblyTarget::RHS.
   */
  TEST(PETSc_MPI_KSPFactorizationReuse, MultipleRightHandSidesFactorizeOnce)
  {
    if (unsupported())
      GTEST_SKIP() << "PETSc event logging or MUMPS is unavailable.";

    Context::MPI ctx(*g_env, *g_world);
    auto mesh = distributeFromRoot(ctx);

    P1<Real, Mesh<Context::MPI>> fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);

    constexpr PetscReal alpha = 5.0;

    Problem mass(u, v);
    Rodin::Solver::KSP ksp(mass);

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
}

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env = &env;
  g_world = &world;

  PetscInitialize(&argc, &argv, nullptr, nullptr);

  // Drive every solve in this file with the same distributed direct solver the
  // production coronary example uses. icntl_20/21 = 0 keep the right-hand side
  // and solution centralized on the host; this MUMPS build faults in its
  // distributed RHS scatter otherwise.
  const auto setDefault = [](const char* name, const char* value) {
    PetscBool set = PETSC_FALSE;
    PetscErrorCode ierr = PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name, &set);
    assert(ierr == PETSC_SUCCESS);
    if (!set)
    {
      ierr = PetscOptionsSetValue(PETSC_NULLPTR, name, value);
      assert(ierr == PETSC_SUCCESS);
    }
    (void)ierr;
  };
  setDefault("-ksp_type", "preonly");
  setDefault("-pc_type", "lu");
  setDefault("-pc_factor_mat_solver_type", "mumps");
  setDefault("-mat_mumps_icntl_20", "0");
  setDefault("-mat_mumps_icntl_21", "0");

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
