#ifndef RODIN_SOLVER_PETSC_GMRES_H
#define RODIN_SOLVER_PETSC_GMRES_H

/**
 * @file
 * @brief PETSc specialization of the GMRES solver.
 */

#include <petscksp.h>

#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "KSP.h"

namespace Rodin::Solver
{
  /**
   * @ingroup GMRESSpecializations
   * @brief Generalized minimal residual solver for PETSc linear systems.
   *
   * This specialization targets @ref Rodin::PETSc::Math::LinearSystem and
   * configures the underlying PETSc KSP type to `KSPGMRES`.
   */
  template <>
  class GMRES<PETSc::Math::LinearSystem> final : public KSP
  {
    public:
      using OperatorType = ::Mat;
      using VectorType = ::Vec;
      using ScalarType = PetscScalar;
      using LinearSystemType = PETSc::Math::LinearSystem;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using Parent = KSP;
      using Parent::solve;

      /**
       * @brief Construct a GMRES solver and set the PETSc type to KSPGMRES.
       * @param pb The variational problem to solve.
       */
      explicit GMRES(ProblemBaseType& pb);

      /**
       * @brief Copy constructor.
       * @param other Another GMRES instance.
       */
      GMRES(const GMRES& other);

      /**
       * @brief Move constructor.
       * @param other Moved-from GMRES instance.
       */
      GMRES(GMRES&& other);

      /// @brief Creates a heap-allocated copy of this GMRES solver.
      /// @returns Pointer to the cloned GMRES instance.
      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }
  };
}

namespace Rodin::PETSc::Solver
{
  /**
   * @brief Convenient PETSc alias for
   * @ref Rodin::Solver::GMRES<Rodin::PETSc::Math::LinearSystem>.
   */
  using GMRES = Rodin::Solver::GMRES<Math::LinearSystem>;
}

#endif // RODIN_SOLVER_PETSC_GMRES_H
