/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SparseQR_H
#define RODIN_SOLVER_SparseQR_H

#include <Eigen/SparseQR>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SparseQR
   */
  SparseQR() -> SparseQR<Math::SparseMatrix, Math::Vector>;

  /**
   * @defgroup SparseQRSpecializations SparseQR Template Specializations
   * @brief Template specializations of the SparseQR class.
   * @see SparseQR
   */

  /**
   * @ingroup SparseQRSpecializations
   * @brief Sparse left-looking QR factorization with numerical column pivoting for use with Math::SparseMatrix and
   * Math::Vector.
   */
  template <>
  class SparseQR<Math::SparseMatrix, Math::Vector> final
    : public SolverBase<Math::SparseMatrix, Math::Vector>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::SparseMatrix;

      /// Type of vector
      using VectorType = Math::Vector;

      /**
       * @brief Constructs the SparseQR object with default parameters.
       */
      SparseQR() = default;

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

    private:
      Eigen::SparseQR<Math::SparseMatrix, Eigen::COLAMDOrdering<int>> m_solver;
  };

}

#endif



