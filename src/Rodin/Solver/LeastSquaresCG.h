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
   * @defgroup LeastSquaresCGSpecializations LeastSquaresCG Template Specializations
   * @brief Template specializations of the LeastSquaresCG class.
   * @see LeastSquaresCG
   */

  /**
   * @ingroup LeastSquaresCGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class LeastSquaresCG<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      LeastSquaresCG(ProblemType& pb)
        : Parent(pb)
      {}

      LeastSquaresCG(const LeastSquaresCG& other)
        : Parent(other)
      {}

      LeastSquaresCG(LeastSquaresCG&& other)
        : Parent(std::move(other))
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
      Eigen::LeastSquaresConjugateGradient<OperatorType> m_solver;
  };

  template <class Scalar>
  LeastSquaresCG(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> LeastSquaresCG<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;

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
  template <class Scalar>
  class LeastSquaresCG<Math::Matrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::Matrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      /// Type of vector
      using VectorType = Math::Vector<Scalar>;

      /// Type of linear operator
      using OperatorType = Math::Matrix<Scalar>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      LeastSquaresCG(ProblemType& pb)
        : Parent(pb)
      {}

      LeastSquaresCG(const LeastSquaresCG& other)
        : Parent(other)
      {}

      LeastSquaresCG(LeastSquaresCG&& other)
        : Parent(std::move(other))
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
      Eigen::LeastSquaresConjugateGradient<OperatorType> m_solver;
  };

  template <class Scalar>
  LeastSquaresCG(Variational::ProblemBase<Math::Matrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> LeastSquaresCG<Math::Matrix<Scalar>, Math::Vector<Scalar>>;
}

#endif


