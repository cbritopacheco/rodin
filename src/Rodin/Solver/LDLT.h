/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_LDLT_H
#define RODIN_SOLVER_LDLT_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup LDLTSpecializations LDLT Template Specializations
   * @brief Template specializations of the LDLT class.
   * @see LDLT
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for LDLT
   */
  template <class LinearSystem>
  LDLT(Variational::ProblemBase<LinearSystem>&) -> LDLT<LinearSystem>;

  /**
   * @ingroup LDLTSpecializations
   * @brief A direct sparse LDLT Cholesky factorizations without square root
   * for use with Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <class Scalar>
  class LDLT<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public SolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::Matrix<ScalarType>;

      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      using Parent = SolverBase<LinearSystemType>;

      using Parent::solve;

      LDLT(ProblemBaseType& pb)
        : Parent(pb)
      {}

      LDLT(const LDLT& other)
        : Parent(other)
      {}

      LDLT(LDLT&& other)
        : Parent(std::move(other))
      {}

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      LDLT* copy() const noexcept override
      {
        return new LDLT(*this);
      }

    private:
      Eigen::LDLT<OperatorType> m_solver;
  };
}

#endif



