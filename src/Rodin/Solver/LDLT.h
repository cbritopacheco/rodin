/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_LDLT_H
#define RODIN_SOLVER_LDLT_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

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
   * @ingroup LDLTSpecializations
   * @brief A direct sparse LDLT Cholesky factorizations without square root
   * for use with Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <class Scalar>
  class LDLT<Math::Matrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::Matrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::Matrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      LDLT(ProblemType& pb)
        : Parent(pb)
      {}

      LDLT(const LDLT& other)
        : Parent(other)
      {}

      LDLT(LDLT&& other)
        : Parent(std::move(other))
      {}

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        x = m_solver.compute(A).solve(b);
      }

      inline
      LDLT* copy() const noexcept override
      {
        return new LDLT(*this);
      }

    private:
      Eigen::LDLT<OperatorType> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for LDLT
   */
  template <class Scalar>
  LDLT(Variational::ProblemBase<Math::Matrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> LDLT<Math::Matrix<Scalar>, Math::Vector<Scalar>>;
}

#endif



