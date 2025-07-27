/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_IDRSTABL_H
#define RODIN_SOLVER_IDRSTABL_H

#include <unsupported/Eigen/IterativeSolvers>

#include "Rodin/Configure.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup IDRSTABLSpecializations IDRSTABL Template Specializations
   * @brief Template specializations of the IDRSTABL class.
   * @see IDRSTABL
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for IDRSTABL
   */
  template <class LinearSystem>
  IDRSTABL(Variational::ProblemBase<LinearSystem>&) -> IDRSTABL<LinearSystem>;

  /**
   * @ingroup IDRSTABLSpecializations
   * @brief IDRSTABL for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class IDRSTABL<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
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

      IDRSTABL(ProblemBaseType& pb)
        : Parent(pb)
      {}

      IDRSTABL(const IDRSTABL& other)
        : Parent(other)
      {}

      IDRSTABL(IDRSTABL&& other)
        : Parent(std::move(other))
      {}

      ~IDRSTABL() = default;

      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      IDRSTABL& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      IDRSTABL& setTolerance(Real tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      IDRSTABL* copy() const noexcept override
      {
        return new IDRSTABL(*this);
      }

    private:
      Eigen::IDRSTABL<OperatorType> m_solver;
  };
}

#endif




