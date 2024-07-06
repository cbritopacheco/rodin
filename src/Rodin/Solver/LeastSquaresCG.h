/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_LeastSquaresCG_H
#define RODIN_SOLVER_LeastSquaresCG_H

#include <optional>
#include <functional>
#include <Eigen/SparseCholesky>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @ingroup RodinCTAD
   * @brief CTAD for LeastSquaresCG
   */
  LeastSquaresCG() -> LeastSquaresCG<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  /**
   * @defgroup LeastSquaresCGSpecializations LeastSquaresCG Template Specializations
   * @brief Template specializations of the LeastSquaresCG class.
   * @see LeastSquaresCG
   */

  /**
   * @ingroup LeastSquaresCGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <>
  class LeastSquaresCG<Math::SparseMatrix<Real>, Math::Vector<Real>> final
    : public SolverBase<Math::SparseMatrix<Real>, Math::Vector<Real>>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::SparseMatrix<Real>;

      /// Type of vector
      using VectorType = Math::Vector<Real>;

      /**
       * @brief Constructs the LeastSquaresCG object with default parameters.
       */
      LeastSquaresCG() = default;

      LeastSquaresCG(const LeastSquaresCG& other)
      {}

      ~LeastSquaresCG() = default;

      LeastSquaresCG& setTolerance(double tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      LeastSquaresCG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        x = m_solver.compute(A).solve(b);
      }

      inline
      LeastSquaresCG* copy() const noexcept override
      {
        return new LeastSquaresCG(*this);
      }

    private:
      Eigen::LeastSquaresConjugateGradient<Math::SparseMatrix<Real>> m_solver;
  };

  /**
   * @defgroup LeastSquaresCGSpecializations LeastSquaresCG Template Specializations
   * @brief Template specializations of the LeastSquaresCG class.
   * @see LeastSquaresCG
   */

  /**
   * @ingroup LeastSquaresCGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <>
  class LeastSquaresCG<Math::Matrix<Real>, Math::Vector<Real>> final
    : public SolverBase<Math::Matrix<Real>, Math::Vector<Real>>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::Matrix<Real>;

      /// Type of vector
      using VectorType = Math::Vector<Real>;

      /**
       * @brief Constructs the LeastSquaresCG object with default parameters.
       */
      LeastSquaresCG() = default;

      LeastSquaresCG(const LeastSquaresCG& other)
      {}

      ~LeastSquaresCG() = default;

      LeastSquaresCG& setTolerance(double tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      LeastSquaresCG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

      inline
      LeastSquaresCG* copy() const noexcept override
      {
        return new LeastSquaresCG(*this);
      }

    private:
      Eigen::LeastSquaresConjugateGradient<Math::Matrix<Real>> m_solver;
  };
}

#endif


