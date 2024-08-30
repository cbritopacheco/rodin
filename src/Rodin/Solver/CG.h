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
   * @defgroup CGSpecializations CG Template Specializations
   * @brief Template specializations of the CG class.
   * @see CG
   */

  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class CG<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      /**
       * @brief Constructs the CG object with default parameters.
       */
      CG(ProblemType& pb)
        : Parent(pb)
      {}

      CG(const CG& other)
        : Parent(other)
      {}

      CG(CG&& other)
        : Parent(std::move(other))
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

      bool success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      CG* copy() const noexcept override
      {
        return new CG(*this);
      }

    private:
      Eigen::ConjugateGradient<OperatorType, Eigen::Lower | Eigen::Upper> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for CG
   */
  template <class Scalar>
  CG(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> CG<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;

  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::Matrix and Math::Vector.
   */
  template <class Scalar>
  class CG<Math::Matrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::Matrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::Matrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      CG(ProblemType& pb)
        : Parent(pb)
      {}

      CG(const CG& other)
        : Parent(other)
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

      bool success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      CG* copy() const noexcept override
      {
        return new CG(*this);
      }

    private:
      Eigen::ConjugateGradient<OperatorType, Eigen::Lower | Eigen::Upper> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for CG
   */
  template <class Scalar>
  CG(Variational::ProblemBase<Math::Matrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> CG<Math::Matrix<Scalar>, Math::Vector<Scalar>>;
}

#endif

