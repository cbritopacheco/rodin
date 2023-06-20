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
#include<Eigen/SparseCholesky>

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
  class CG<Math::SparseMatrix, Math::Vector> final
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
      constexpr
      CG()
        : m_tolerance(std::numeric_limits<double>::epsilon())
      {}

      ~CG() = default;

      constexpr
      CG& setTolerance(double tol)
      {
        m_tolerance = tol;
        return *this;
      }

      constexpr
      CG& setMaxIterations(size_t maxIt)
      {
        m_maxIterations = maxIt;
        return *this;
      }

      void solve(OperatorType& A, VectorType& x, VectorType& b) const override
      {
        Eigen::ConjugateGradient<Math::SparseMatrix> solver;

        solver.setTolerance(m_tolerance);

        if (m_maxIterations.has_value())
          solver.setMaxIterations(m_maxIterations.value());

        x = solver.compute(A).solve(b);
      }

    private:
      double m_tolerance;
      std::optional<size_t> m_maxIterations;
  };

}

#endif

