/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_CG_H
#define RODIN_SOLVER_CG_H

#include <optional>
#include <functional>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @ingroup RodinCTAD
   * @brief CTAD for CG
   */
  CG() -> CG<Math::SparseMatrix, Math::Vector>;

  /**
   * @defgroup CGSpecializations CG Template Specializations
   * @brief Template specializations of the CG class.
   * @see CG
   */

  /**
   * @ingroup CGSpecializations
   * @brief Conjugate Gradient for use with Math::SparseMatrix and
   * Math::Vector.
   */
  template <>
  class CG<Math::SparseMatrix, Math::Vector>
    : public SolverBase<Math::SparseMatrix, Math::Vector>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::SparseMatrix;

      /// Type of vector
      using VectorType = Math::Vector;

      /**
       * @brief Constructs the CG object with default parameters.
       */
      CG()
      {}

      ~CG() = default;

      virtual
      void solve(OperatorType& A, VectorType& x, VectorType& b)
      const override
      {
        Eigen::ConjugateGradient<Math::SparseMatrix> solver;
        x = solver.compute(A).solve(b);
      }

    private:
  };

}

#endif

