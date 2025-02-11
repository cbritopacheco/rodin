/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_IDRSTABL_H
#define RODIN_SOLVER_IDRSTABL_H

#include <optional>
#include <functional>

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
   * @ingroup IDRSTABLSpecializations
   * @brief IDRSTABL for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class IDRSTABL<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      IDRSTABL(ProblemType& pb)
        : Parent(pb)
      {}

      IDRSTABL(const IDRSTABL& other)
        : Parent(other)
      {}

      IDRSTABL(IDRSTABL&& other)
        : Parent(std::move(other))
      {}

      ~IDRSTABL() = default;

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

      IDRSTABL& setMaxIterations(size_t it)
      {
        m_solver.setMaxIterations(it);
        return *this;
      }

      IDRSTABL& setTolerance(Real tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      bool success() const
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

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for IDRSTABL
   */
  template <class Scalar>
  IDRSTABL(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> IDRSTABL<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;
}

#endif




