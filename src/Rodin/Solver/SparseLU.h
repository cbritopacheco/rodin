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
  SparseLU() -> SparseLU<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  /**
   * @defgroup SparseLUSpecializations SparseLU Template Specializations
   * @brief Template specializations of the SparseLU class.
   * @see SparseLU
   */

  /**
   * @ingroup SparseLUSpecializations
   * @brief Sparse supernodal LU factorization for general matrices for use
   * with Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <>
  class SparseLU<Math::SparseMatrix<Real>, Math::Vector<Real>> final
    : public SolverBase<Math::SparseMatrix<Real>, Math::Vector<Real>>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::SparseMatrix<Real>;

      /// Type of vector
      using VectorType = Math::Vector<Real>;

      /**
       * @brief Constructs the SparseLU object with default parameters.
       */
      SparseLU() = default;

      SparseLU(const SparseLU& other)
      {}

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

      inline
      SparseLU* copy() const noexcept override
      {
        return new SparseLU(*this);
      }

    private:
      Eigen::SparseLU<Math::SparseMatrix<Real>> m_solver;
  };

}

#endif



