#ifndef RODIN_SOLVER_PETSC_CG_H
#define RODIN_SOLVER_PETSC_CG_H

/**
 * @file
 * @brief PETSc specialization of the conjugate gradient (CG) solver.
 */

#include <petscksp.h>

#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "KSP.h"

namespace Rodin::Solver
{
  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for PETSc linear systems.
   *
   * This specialization targets @ref Rodin::PETSc::Math::LinearSystem and
   * configures the underlying PETSc KSP type to `KSPCG`.
   */
  template <>
  class CG<PETSc::Math::LinearSystem> final : public KSP
  {
    public:
      /// @brief PETSc matrix operator type.
      using OperatorType = ::Mat;
      /// @brief PETSc vector type.
      using VectorType = ::Vec;
      /// @brief Scalar type (PETSc scalar).
      using ScalarType = PetscScalar;
      /// @brief Linear system type for PETSc solvers.
      using LinearSystemType = PETSc::Math::LinearSystem;
      /// @brief Base problem type.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      /// @brief Parent class type.
      using Parent = KSP;
      using Parent::solve;

      /**
       * @brief Construct a CG solver and set the PETSc type to KSPCG.
       * @param pb The variational problem to solve.
       */
      explicit CG(ProblemBaseType& pb);

      /**
       * @brief Copy constructor.
       * @param other Another CG instance.
       */
      CG(const CG& other);

      /**
       * @brief Move constructor.
       * @param other Moved-from CG instance.
       */
      CG(CG&& other);

      /// @brief Creates a heap-allocated copy of this CG solver.
      /// @returns Pointer to the cloned CG instance.
      CG* copy() const noexcept override
      {
        return new CG(*this);
      }
  };
}

namespace Rodin::PETSc::Solver
{
  /**
   * @brief Convenient PETSc alias for
   * @ref Rodin::Solver::CG<Rodin::PETSc::Math::LinearSystem>.
   */
  using CG = Rodin::Solver::CG<Math::LinearSystem>;
}

#endif // RODIN_SOLVER_PETSC_CG_H
