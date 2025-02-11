/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_GMRES_H
#define RODIN_SOLVER_GMRES_H

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
   * @defgroup GMRESSpecializations GMRES Template Specializations
   * @brief Template specializations of the GMRES class.
   * @see GMRES
   */

  /**
   * @ingroup GMRESSpecializations
   * @brief GMRES for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class GMRES<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      GMRES(ProblemType& pb)
        : Parent(pb)
      {}

      GMRES(const GMRES& other)
        : Parent(other)
      {}

      GMRES(GMRES&& other)
        : Parent(std::move(other))
      {}

      ~GMRES() = default;

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

      GMRES& setTolerance(Real tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      GMRES& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      GMRES& setRestart(size_t restart)
      {
        m_solver.set_restart(restart);
        return *this;
      }

      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }

    private:
      Eigen::GMRES<OperatorType> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for GMRES
   */
  template <class Scalar>
  GMRES(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> GMRES<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;
}

#endif



