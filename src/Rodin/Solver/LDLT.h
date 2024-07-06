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
   * @ingroup RodinCTAD
   * @brief CTAD for LDLT
   */
  LDLT() -> LDLT<Math::Matrix<Real>, Math::Vector<Real>>;

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
  template <>
  class LDLT<Math::Matrix<Real>, Math::Vector<Real>> final
    : public SolverBase<Math::Matrix<Real>, Math::Vector<Real>>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::Matrix<Real>;

      /// Type of vector
      using VectorType = Math::Vector<Real>;

      /**
       * @brief Constructs the LDLT object with default parameters.
       */
      LDLT() = default;

      LDLT(const LDLT& other)
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
      Eigen::LDLT<Math::Matrix<Real>> m_solver;
  };

}

#endif



