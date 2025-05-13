/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_PETSC_KSP_H
#define RODIN_SOLVER_PETSC_KSP_H

#include <petscksp.h>
#include "Rodin/Solver/Solver.h"
#include "Rodin/PETSc/Object.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Solver
{
  /**
   * @brief PETSc KSP (Krylov) linear solver wrapper.
   *
   * Inherits SolverBase<Mat,Vec,PetscScalar> for the generic interface,
   * and PETSc::Object for automatic cleanup of any forgotten handles.
   *
   * Combines programmatic configuration with command‐line overrides.
   */
  class KSP : public SolverBase< ::Mat, ::Vec, PetscScalar>, public PETSc::Object
  {
  public:
    using OperatorType = ::Mat;
    using VectorType   = ::Vec;
    using ScalarType   = PetscScalar;
    using Parent       = SolverBase<OperatorType, VectorType, ScalarType>;
    using ProblemType  = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;
    using Parent::solve;

    /**
     * @brief Construct and create the PETSc KSP object.
     *
     * Initializes to PETSC defaults.
     *
     * @param pb   Variational problem this solver will solve.
     * @param comm MPI communicator (default PETSC_COMM_WORLD).
     */
    explicit KSP(ProblemType& pb);

    ~KSP() override;

    ::PetscObject& getHandle() noexcept override;

    /**
     * @brief Solve @f$ Ax = b @f$, allocating @f$ x @f$ if null.
     *
     * If `x == PETSC_NULL`, automatically `VecDuplicate(b, &x)` and zero it.
     * Otherwise uses `x` contents as initial guess.
     *
     * Applies programmatic settings, then SetFromOptions, then KSPSolve.
     *
     * @param A Left‐hand side matrix.
     * @param x Solution vector (initial guess in; may be PETSC_NULL).
     * @param b Right‐hand side vector.
     */
    void solve(OperatorType& A, VectorType& x, VectorType& b) override;

    KSP& setType(::KSPType type) noexcept;

    KSP& setTolerances(PetscReal rtol,
                       PetscReal abstol,
                       PetscReal dtol,
                       PetscInt  maxIt) noexcept;

    KSP& setOperators(OperatorType A, OperatorType P = nullptr) noexcept;

  private:
    ::KSP       m_ksp{nullptr};
    ::KSPType   m_type{KSPGMRES};
    PetscReal   m_rtol   {PETSC_DEFAULT},
                m_abstol {PETSC_DEFAULT},
                m_dtol   {PETSC_DEFAULT};
    PetscInt    m_maxIt  {PETSC_DEFAULT};
    OperatorType m_A     {nullptr},
                 m_P     {nullptr};
  };

} // namespace Rodin::Solver

#endif // RODIN_SOLVER_PETSC_KSP_H
