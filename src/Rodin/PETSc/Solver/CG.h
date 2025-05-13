// CG.h
#ifndef RODIN_SOLVER_PETSC_CG_H
#define RODIN_SOLVER_PETSC_CG_H

#include <petscksp.h>

#include "Rodin/Solver/CG.h"
#include "KSP.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Solver
{
  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * ::Mat and ::Vec.
   */
  template <>
  class CG<::Mat, ::Vec> final : public KSP
  {
  public:
    using OperatorType = ::Mat;
    using VectorType   = ::Vec;
    using ScalarType   = PetscScalar;
    using ProblemType  = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;
    using Parent       = KSP;
    using Parent::solve;

    /**
     * @brief Construct a CG solver and set the PETSc type to KSPCG.
     * @param pb The variational problem to solve.
     */
    explicit CG(ProblemType& pb);

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

    /**
     * @brief Destructor.
     */
    ~CG() override;

    CG* copy() const noexcept override
    {
      return new CG(*this);
    }
  };

  /**
   * @ingroup RodinCTAD
   * @brief Class template argument deduction for CG from ProblemBase.
   */
  CG(Variational::ProblemBase<::Mat, ::Vec, PetscScalar>&) -> CG<::Mat, ::Vec>;

} // namespace Rodin::Solver

#endif // RODIN_SOLVER_PETSC_CG_H
