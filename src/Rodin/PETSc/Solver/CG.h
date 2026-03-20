#ifndef RODIN_SOLVER_PETSC_CG_H
#define RODIN_SOLVER_PETSC_CG_H

/**
 * @file CG.h
 * @brief PETSc specialization of the conjugate gradient (CG) solver.
 *
 * Provides the `CG<PETSc::Math::LinearSystem>` specialization that sets
 * the PETSc KSP type to `KSPCG`.  The conjugate gradient method is
 * applicable to symmetric positive definite systems
 * @f$ A\mathbf{x} = \mathbf{b} @f$.
 *
 * @see Rodin::Solver::KSP, Rodin::PETSc::Solver::GMRES
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
   * Specialization of @ref Rodin::Solver::CG that targets
   * @ref Rodin::PETSc::Math::LinearSystem and sets the underlying KSP
   * algorithm to `KSPCG`.  Applicable only to symmetric positive definite
   * (SPD) systems.
   *
   * @see Rodin::Solver::KSP, Rodin::PETSc::Solver::CG
   */
  template <>
  class CG<PETSc::Math::LinearSystem> final : public KSP
  {
    public:
      /// @brief PETSc matrix type (`::Mat`) for the system operator.
      using OperatorType = ::Mat;
      /// @brief PETSc vector type (`::Vec`) for the RHS and solution.
      using VectorType = ::Vec;
      /// @brief Scalar type (`PetscScalar`).
      using ScalarType = PetscScalar;
      /// @brief Linear system type coupling @f$ A @f$, @f$ \mathbf{b} @f$, and @f$ \mathbf{x} @f$.
      using LinearSystemType = PETSc::Math::LinearSystem;
      /// @brief Base problem type that provides the linear system.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      /// @brief Parent class type (@ref Rodin::Solver::KSP).
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
