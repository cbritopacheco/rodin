/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SIMPLICIALLLT_H
#define RODIN_SOLVER_SIMPLICIALLLT_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SimplicialLLT
   */
  SimplicialLLT() -> SimplicialLLT<Math::SparseMatrix, Math::Vector>;

  /**
   * @defgroup SimplicialLLTSpecializations SimplicialLLT Template Specializations
   * @brief Template specializations of the SimplicialLLT class.
   * @see SimplicialLLT
   */

  /**
   * @ingroup SimplicialLLTSpecializations
   * @brief A direct sparse LLT Cholesky factorizations for use with
   * Math::SparseMatrix and Math::Vector.
   */
  template <>
  class SimplicialLLT<Math::SparseMatrix, Math::Vector> final
    : public SolverBase<Math::SparseMatrix, Math::Vector>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::SparseMatrix;

      /// Type of vector
      using VectorType = Math::Vector;

      /**
       * @brief Constructs the SimplicialLLT object with default parameters.
       */
      SimplicialLLT() = default;

      SimplicialLLT(const SimplicialLLT& other)
      {}

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        x = m_solver.compute(A).solve(b);
      }

      inline
      SimplicialLLT* copy() const noexcept override
      {
        return new SimplicialLLT(*this);
      }

    private:
      Eigen::SimplicialLLT<Math::SparseMatrix> m_solver;
  };

}

#endif



