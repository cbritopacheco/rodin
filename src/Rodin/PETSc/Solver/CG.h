#ifndef RODIN_SOLVER_PETSC_CG_H
#define RODIN_SOLVER_PETSC_CG_H

#include <petscksp.h>

#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Solver/CG.h"
#include "KSP.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Solver
{
  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * PETSc::Matrix and PETSc::Vector.
   */
  template <>
  class CG<PETSc::LinearSystem> final : public KSP
  {
  public:
    using OperatorType = PETSc::Matrix;
    using VectorType = PETSc::Vector;
    using ScalarType = PetscScalar;
    using LinearSystemType = PETSc::LinearSystem;
    using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
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

    CG* copy() const noexcept override
    {
      return new CG(*this);
    }
  };
} // namespace Rodin::Solver

#endif // RODIN_SOLVER_PETSC_CG_H
