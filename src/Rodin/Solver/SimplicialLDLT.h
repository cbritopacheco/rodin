/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SIMPLICIALLDLT_H
#define RODIN_SOLVER_SIMPLICIALLDLT_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SimplicialLDLT
   */
  SimplicialLDLT() -> SimplicialLDLT<Math::SparseMatrix, Math::Vector>;

  /**
   * @defgroup SimplicialLDLTSpecializations SimplicialLDLT Template Specializations
   * @brief Template specializations of the SimplicialLDLT class.
   * @see SimplicialLDLT
   */

  /**
   * @ingroup SimplicialLDLTSpecializations
   * @brief A direct sparse LDLT Cholesky factorizations without square root
   * for use with Math::SparseMatrix and Math::Vector.
   */
  template <>
  class SimplicialLDLT<Math::SparseMatrix, Math::Vector> final
    : public SolverBase<Math::SparseMatrix, Math::Vector>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::SparseMatrix;

      /// Type of vector
      using VectorType = Math::Vector;

      /**
       * @brief Constructs the SimplicialLDLT object with default parameters.
       */
      SimplicialLDLT() = default;

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        x = m_solver.compute(A).solve(b);
      }

    private:
      Eigen::ConjugateGradient<Math::SparseMatrix> m_solver;
  };

}

#endif


