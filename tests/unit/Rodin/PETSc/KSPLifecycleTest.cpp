/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <petsc.h>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include "Rodin/Geometry.h"
#include "Rodin/PETSc.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

template <class SolverType>
class KSPLifecycle : public ::testing::Test
{
  protected:
    using Scalar = typename SolverType::ScalarType;
    using Space = P1<Scalar>;
    using Trial =
      decltype(PETSc::Variational::TrialFunction(std::declval<const Space&>()));
    using Test = decltype(PETSc::Variational::TestFunction(std::declval<const Space&>()));
    using ProblemType = decltype(Problem(std::declval<Trial&>(), std::declval<Test&>()));

    KSPLifecycle()
      : m_mesh(LocalMesh::UniformGrid(Polytope::Type::Segment, {4})),
        m_space(m_mesh),
        m_trial(m_space),
        m_test(m_space),
        m_problem(m_trial, m_test)
    {}

    void SetUp() override
    {
      const auto source = [&] {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return ComplexFunction(expected());
        else
          return RealFunction(expected());
      }();
      m_problem = Integral(m_trial, m_test) - Integral(source, m_test);
      m_problem.assemble();
      ASSERT_EQ(MatDuplicate(m_problem.getLinearSystem().getOperator(), MAT_COPY_VALUES,
                  &m_preconditioner.get()),
        PETSC_SUCCESS);
    }

    /// @brief Returns the exact constant mass-system solution.
    static Scalar expected()
    {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -1);
      else
        return Scalar(2);
    }

    /// @brief Returns the expected algorithm for the concrete wrapper.
    static const char* expectedType()
    {
      if constexpr (std::is_same_v<SolverType, PETSc::Solver::CG>)
        return KSPCG;
      else if constexpr (std::is_same_v<SolverType, PETSc::Solver::GMRES>)
        return KSPGMRES;
      else
        return KSPBCGS;
    }

    /// @brief Configures observable programmatic settings and a distinct borrowed preconditioner matrix.
    void configure(SolverType& solver)
    {
      solver.setType(expectedType())
        .setPrefix(std::string("lifecycle_"))
        .setTolerances(
          relativeTolerance, absoluteTolerance, divergenceTolerance, maximumIterations)
        .setPreconditioner(m_preconditioner.get());
      ::PC pc = nullptr;
      ASSERT_EQ(KSPGetPC(solver.getHandle(), &pc), PETSC_SUCCESS);
      ASSERT_EQ(PCSetType(pc, PCJACOBI), PETSC_SUCCESS);
    }

    /// @brief Checks the exact solution, solver configuration and borrowed problem after a fresh solve.
    void check(SolverType& solver)
    {
      auto& system = m_problem.getLinearSystem();
      ASSERT_EQ(VecSet(system.getSolution(), 0), PETSC_SUCCESS);
      solver.solve(system);
      PetscInt size = 0;
      const PetscScalar* values = nullptr;
      ASSERT_EQ(VecGetLocalSize(system.getSolution(), &size), PETSC_SUCCESS);
      ASSERT_EQ(VecGetArrayRead(system.getSolution(), &values), PETSC_SUCCESS);
      for (PetscInt index = 0; index < size; ++index)
        EXPECT_LE(std::abs(values[index] - expected()), solutionTolerance);
      ASSERT_EQ(VecRestoreArrayRead(system.getSolution(), &values), PETSC_SUCCESS);
      ::KSPConvergedReason reason;
      ASSERT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
      EXPECT_GT(reason, 0);
      EXPECT_GT(solver.getIterationNumber(), 0u);
      EXPECT_LE(solver.getError(), solutionTolerance);
      PetscReal relative = 0, absolute = 0, divergence = 0;
      PetscInt iterations = 0;
      ASSERT_EQ(KSPGetTolerances(
                  solver.getHandle(), &relative, &absolute, &divergence, &iterations),
        PETSC_SUCCESS);
      EXPECT_EQ(relative, relativeTolerance);
      EXPECT_EQ(absolute, absoluteTolerance);
      EXPECT_EQ(divergence, divergenceTolerance);
      EXPECT_EQ(iterations, maximumIterations);
      const char* type = nullptr;
      const char* prefix = nullptr;
      ASSERT_EQ(KSPGetType(solver.getHandle(), &type), PETSC_SUCCESS);
      ASSERT_EQ(KSPGetOptionsPrefix(solver.getHandle(), &prefix), PETSC_SUCCESS);
      EXPECT_STREQ(type, expectedType());
      EXPECT_STREQ(prefix, "lifecycle_");
      ::Mat matrix = nullptr, preconditioner = nullptr;
      ASSERT_EQ(
        KSPGetOperators(solver.getHandle(), &matrix, &preconditioner), PETSC_SUCCESS);
      EXPECT_EQ(matrix, system.getOperator());
      EXPECT_EQ(preconditioner, m_preconditioner.get());
      EXPECT_EQ(&solver.getProblem(), &m_problem);
      EXPECT_EQ(&std::as_const(solver).getProblem(), &m_problem);
      EXPECT_EQ(std::as_const(solver).getHandle(), solver.getHandle());
      MPI_Comm communicator = MPI_COMM_NULL;
      std::as_const(solver).getComm(&communicator);
      int comparison = MPI_UNEQUAL;
      ASSERT_EQ(MPI_Comm_compare(communicator, system.getCommunicator(), &comparison),
        MPI_SUCCESS);
      EXPECT_TRUE(comparison == MPI_IDENT || comparison == MPI_CONGRUENT);
    }

    static constexpr PetscReal relativeTolerance = 1e-12;
    static constexpr PetscReal absoluteTolerance = 1e-14;
    static constexpr PetscReal divergenceTolerance = 1e4;
    static constexpr PetscInt maximumIterations = 50;
    static constexpr PetscReal solutionTolerance = 1e-11;
    LocalMesh m_mesh;
    Space m_space;
    Trial m_trial;
    Test m_test;
    ProblemType m_problem;
    NamedFormStorage<::Mat> m_preconditioner;
};

using KSPSolverTypes =
  ::testing::Types<Solver::KSP, PETSc::Solver::CG, PETSc::Solver::GMRES>;
TYPED_TEST_SUITE(KSPLifecycle, KSPSolverTypes);

/// @brief Copies and polymorphic clones retain the shared context and independent programmatic settings across source destruction.
TYPED_TEST(KSPLifecycle, CopiesAndClonesSurviveSourceDestruction)
{
  std::unique_ptr<TypeParam> survivor;
  {
    TypeParam original(this->m_problem);
    this->configure(original);
    this->check(original);
    auto copied = std::make_unique<TypeParam>(original);
    EXPECT_EQ(copied->getHandle(), original.getHandle());
    PetscInt references = 0;
    ASSERT_EQ(PetscObjectGetReference(
                reinterpret_cast<PetscObject>(original.getHandle()), &references),
      PETSC_SUCCESS);
    EXPECT_EQ(references, 2);
    survivor.reset(copied->copy());
    ASSERT_EQ(PetscObjectGetReference(
                reinterpret_cast<PetscObject>(original.getHandle()), &references),
      PETSC_SUCCESS);
    EXPECT_EQ(references, 3);
    copied.reset();
    original.setType(std::is_same_v<TypeParam, PETSc::Solver::CG> ? KSPGMRES : KSPCG)
      .setPrefix(std::string("changed_"))
      .setTolerances(1e-5, 1e-8, 1e7, 10)
      .setPreconditioner(this->m_problem.getLinearSystem().getOperator());
    original.solve(this->m_problem.getLinearSystem());
    this->check(*survivor);
  }
  this->check(*survivor);
}

/// @brief Moves transfer the raw handle, retain settings, and permit destruction or copying of the empty source.
TYPED_TEST(KSPLifecycle, MoveTransfersHandleAndPreservesConfiguration)
{
  TypeParam original(this->m_problem);
  this->configure(original);
  this->check(original);
  const auto handle = original.getHandle();
  TypeParam moved(std::move(original));
  EXPECT_EQ(original.getHandle(), nullptr);
  EXPECT_EQ(moved.getHandle(), handle);
  this->check(moved);
  TypeParam emptyCopy(original);
  EXPECT_EQ(emptyCopy.getHandle(), nullptr);
  TypeParam emptyMove(std::move(original));
  EXPECT_EQ(emptyMove.getHandle(), nullptr);
}

/// @brief Named and unnamed stored splits are installed when the runtime preconditioner selects field splitting.
TYPED_TEST(KSPLifecycle, StoredFieldSplitsAreInstalledAtSolve)
{
  auto& system = this->m_problem.getLinearSystem();
  ASSERT_EQ(this->m_space.getSize(), 4);
  ::IS first = nullptr, second = nullptr;
  ASSERT_EQ(ISCreateStride(system.getCommunicator(), 2, 0, 1, &first), PETSC_SUCCESS);
  ASSERT_EQ(ISCreateStride(system.getCommunicator(), 2, 2, 1, &second), PETSC_SUCCESS);
  system.setFieldSplits(
    PETSc::Math::LinearSystem::FieldSplits({{"left", first}, {"", second}}));
  TypeParam solver(this->m_problem);
  this->configure(solver);
  ::PC pc = nullptr;
  ASSERT_EQ(KSPGetPC(solver.getHandle(), &pc), PETSC_SUCCESS);
  ASSERT_EQ(PCSetType(pc, PCFIELDSPLIT), PETSC_SUCCESS);
  ASSERT_EQ(PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE), PETSC_SUCCESS);
  // Fixed positive diagonal block solves satisfy CG's preconditioner contract.
  for (const char* name : {"left", "1"})
  {
    const std::string prefix = std::string("-lifecycle_fieldsplit_") + name;
    ASSERT_EQ(PetscOptionsSetValue(nullptr, (prefix + "_ksp_type").c_str(), "preonly"),
      PETSC_SUCCESS);
    ASSERT_EQ(PetscOptionsSetValue(nullptr, (prefix + "_pc_type").c_str(), "jacobi"),
      PETSC_SUCCESS);
  }
  this->check(solver);
  for (const char* name : {"left", "1"})
  {
    const std::string prefix = std::string("-lifecycle_fieldsplit_") + name;
    ASSERT_EQ(
      PetscOptionsClearValue(nullptr, (prefix + "_ksp_type").c_str()), PETSC_SUCCESS);
    ASSERT_EQ(
      PetscOptionsClearValue(nullptr, (prefix + "_pc_type").c_str()), PETSC_SUCCESS);
  }
  PetscInt count = 0;
  ::KSP* subsolvers = nullptr;
  ASSERT_EQ(PCFieldSplitGetSubKSP(pc, &count, &subsolvers), PETSC_SUCCESS);
  EXPECT_EQ(count, 2);
  ASSERT_EQ(PetscFree(subsolvers), PETSC_SUCCESS);
}

/// @brief Reusing a solver clears obsolete prefixes and PREONLY-incompatible guesses.
TYPED_TEST(KSPLifecycle, ReconfigurationClearsStoredPETScState)
{
  TypeParam solver(this->m_problem);
  this->configure(solver);
  auto& system = this->m_problem.getLinearSystem();
  solver.solve(system);
  solver.setType(KSPPREONLY).setPrefix({});
  ASSERT_EQ(VecSet(system.getSolution(), 123), PETSC_SUCCESS);
  solver.solve(system);
  PetscBool nonzero = PETSC_TRUE;
  ASSERT_EQ(KSPGetInitialGuessNonzero(solver.getHandle(), &nonzero), PETSC_SUCCESS);
  EXPECT_EQ(nonzero, PETSC_FALSE);
  const char* prefix = nullptr;
  ASSERT_EQ(KSPGetOptionsPrefix(solver.getHandle(), &prefix), PETSC_SUCCESS);
  EXPECT_TRUE(prefix == nullptr || prefix[0] == '\0');
  PetscReal norm = 0;
  ASSERT_EQ(VecNorm(system.getSolution(), NORM_2, &norm), PETSC_SUCCESS);
  EXPECT_TRUE(std::isfinite(norm));
  EXPECT_LT(norm, 123);
  solver.setType(KSPGMRES);
  solver.solve(system);
  ASSERT_EQ(KSPGetInitialGuessNonzero(solver.getHandle(), &nonzero), PETSC_SUCCESS);
  EXPECT_EQ(nonzero, PETSC_TRUE);
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
