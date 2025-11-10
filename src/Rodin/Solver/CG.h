/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_CG_H
#define RODIN_SOLVER_CG_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Rodin/Types.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup CGSpecializations CG Template Specializations
   * @brief Template specializations of the CG class.
   * @see CG
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for CG
   */
  template <class LinearSystem>
  CG(Variational::ProblemBase<LinearSystem>&) -> CG<LinearSystem>;

  /**
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class CG<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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

      /**
       * @brief Constructs the CG object with default parameters.
       */
      CG(ProblemBaseType& pb)
        : Parent(pb)
      {}

      CG(const CG& other)
        : Parent(other)
      {}

      CG(CG&& other)
        : Parent(std::move(other))
      {}

      ~CG() = default;

      CG& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      CG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      Boolean success() const
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
   * @ingroup CGSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::Matrix and Math::Vector.
   */
  template <class Scalar>
  class CG<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public SolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::Matrix<ScalarType>;

      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      using ProblemType = Variational::ProblemBase<LinearSystemType>;

      using Parent = SolverBase<LinearSystemType>;

      using Parent::solve;

      CG(ProblemType& pb)
        : Parent(pb)
      {}

      CG(const CG& other)
        : Parent(other)
      {}

      ~CG() = default;

      CG& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      CG& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      Boolean success() const
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
}

#endif

