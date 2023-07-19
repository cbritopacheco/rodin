/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SPARSELU_H
#define RODIN_SOLVER_SPARSELU_H

#include <Eigen/SparseLU>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SparseLU
   */
  SparseLU() -> SparseLU<Math::SparseMatrix, Math::Vector>;

  /**
   * @defgroup SparseLUSpecializations SparseLU Template Specializations
   * @brief Template specializations of the SparseLU class.
   * @see SparseLU
   */

  /**
   * @ingroup SparseLUSpecializations
   * @brief Sparse supernodal LU factorization for general matrices for use
   * with Math::SparseMatrix and Math::Vector.
   */
  template <>
  class SparseLU<Math::SparseMatrix, Math::Vector> final
    : public SolverBase<Math::SparseMatrix, Math::Vector>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::SparseMatrix;

      /// Type of vector
      using VectorType = Math::Vector;

      /**
       * @brief Constructs the SparseLU object with default parameters.
       */
      SparseLU() = default;

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

    private:
      Eigen::SparseLU<Math::SparseMatrix> m_solver;
  };

}

#endif



