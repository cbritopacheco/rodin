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
#include <Eigen/SparseCholesky>

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
  CG() -> CG<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;

  /**
   * @defgroup CGSpecializations CG Template Specializations
   * @brief Template specializations of the CG class.
   * @see CG
   */

  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::SparseMatrix<Scalar> and Math::Vector<Scalar>.
   */
  template <>
  class CG<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::SparseMatrix<Scalar>;

      /// Type of vector
      using VectorType = Math::Vector<Scalar>;

      /**
       * @brief Constructs the CG object with default parameters.
       */
      CG() = default;

      CG(const CG& other)
      {}

      ~CG() = default;

      CG& setTolerance(double tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      CG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        x = m_solver.compute(A).solve(b);
      }

      inline
      bool success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      inline
      CG* copy() const noexcept override
      {
        return new CG(*this);
      }

    private:
      Eigen::ConjugateGradient<
        Math::SparseMatrix<Scalar>, Eigen::Lower | Eigen::Upper> m_solver;
  };

  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::Matrix<Scalar> and Math::Vector<Scalar>.
   */
  template <>
  class CG<Math::Matrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::Matrix<Scalar>, Math::Vector<Scalar>>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::Matrix<Scalar>;

      /// Type of vector
      using VectorType = Math::Vector<Scalar>;

      /**
       * @brief Constructs the CG object with default parameters.
       */
      CG() = default;

      CG(const CG& other)
      {}

      ~CG() = default;

      CG& setTolerance(double tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      CG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        x = m_solver.compute(A).solve(b);
      }

      inline
      bool success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      inline
      CG* copy() const noexcept override
      {
        return new CG(*this);
      }

    private:
      Eigen::ConjugateGradient<
        Math::Matrix<Scalar>, Eigen::Lower | Eigen::Upper> m_solver;
  };

}

#endif

