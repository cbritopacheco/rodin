/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_PETSC_CG_H
#define RODIN_SOLVER_PETSC_CG_H

#include <petscksp.h>

#include "Rodin/Solver/CG.h"
#include "KSP.h"

namespace Rodin::Solver
{
  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * ::Mat and Math::Vector.
   */
  template <>
  class CG<::Mat, ::Vec> final : public KSP
  {
    public:
      using ScalarType = PetscScalar;
      using VectorType = ::Vec;
      using OperatorType = ::Mat;
      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;
      using Parent = KSP;
      using Parent::solve;

      /**
       * @brief Constructs the CG object with default parameters.
       */
      CG(ProblemType& pb)
        : Parent(pb)
      {
        this->setType(KSPCG);
      }

      CG(const CG& other)
        : Parent(other)
      {}

      CG(CG&& other)
        : Parent(std::move(other))
      {}

      ~CG() = default;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for CG
   */
  CG(Variational::ProblemBase<::Mat, ::Vec, PetscScalar>&) -> CG<::Mat, ::Vec>;
}

#endif
