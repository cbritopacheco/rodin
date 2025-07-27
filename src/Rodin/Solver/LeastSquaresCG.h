/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_LEASTSQUARESCG_H
#define RODIN_SOLVER_LEASTSQUARESCG_H

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

  template <class LinearSystem>
  LeastSquaresCG(Variational::ProblemBase<LinearSystem>&) -> LeastSquaresCG<LinearSystem>;

  /**
   * @ingroup LeastSquaresCGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class LeastSquaresCG<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
    : public SolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      using Parent = SolverBase<LinearSystemType>;

      using Parent::solve;

      LeastSquaresCG(ProblemBaseType& pb)
        : Parent(pb)
      {}

      LeastSquaresCG(const LeastSquaresCG& other)
        : Parent(other)
      {}

      LeastSquaresCG(LeastSquaresCG&& other)
        : Parent(std::move(other))
      {}

      ~LeastSquaresCG() = default;

      LeastSquaresCG& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      LeastSquaresCG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      LeastSquaresCG* copy() const noexcept override
      {
        return new LeastSquaresCG(*this);
      }

    private:
      Eigen::LeastSquaresConjugateGradient<OperatorType> m_solver;
  };

  /**
   * @ingroup LeastSquaresCGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <class Scalar>
  class LeastSquaresCG<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public SolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;

      /// Type of vector
      using VectorType = Math::Vector<Scalar>;

      /// Type of linear operator
      using OperatorType = Math::Matrix<Scalar>;

      /// Type of linear system
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      using ProblemType = Variational::ProblemBase<LinearSystemType>;

      using Parent = SolverBase<LinearSystemType>;

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

      LeastSquaresCG& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      LeastSquaresCG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getOperator());
      }

      LeastSquaresCG* copy() const noexcept override
      {
        return new LeastSquaresCG(*this);
      }

    private:
      Eigen::LeastSquaresConjugateGradient<OperatorType> m_solver;
  };
}

#endif


