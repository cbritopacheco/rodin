/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_DGMRES_H
#define RODIN_SOLVER_DGMRES_H

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
   * @defgroup DGMRESSpecializations DGMRES Template Specializations
   * @brief Template specializations of the DGMRES class.
   * @see DGMRES
   */

  /**
   * @ingroup DGMRESSpecializations
   * @brief DGMRES for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class DGMRES<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      DGMRES(ProblemType& pb)
        : Parent(pb)
      {}

      DGMRES(const DGMRES& other)
        : Parent(other)
      {}

      DGMRES(DGMRES&& other)
        : Parent(std::move(other))
      {}

      ~DGMRES() = default;

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

      DGMRES& setTolerance(Real tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      DGMRES& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      DGMRES* copy() const noexcept override
      {
        return new DGMRES(*this);
      }

    private:
      Eigen::DGMRES<OperatorType> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for DGMRES
   */
  template <class Scalar>
  DGMRES(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> DGMRES<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;
}

#endif




