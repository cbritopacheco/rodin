#ifndef RODIN_SOLVER_PETSC_GMRES_H
#define RODIN_SOLVER_PETSC_GMRES_H

#include <petscksp.h>

#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "KSP.h"

namespace Rodin::Solver
{
  /**
   * @ingroup GMRESSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * PETSc::Matrix and PETSc::Vector.
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

      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }
  };
}

namespace Rodin::PETSc::Solver
{
  using GMRES = Rodin::Solver::GMRES<Math::LinearSystem>;
}

#endif // RODIN_SOLVER_PETSC_GMRES_H

