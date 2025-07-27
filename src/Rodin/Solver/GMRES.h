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
#include "Rodin/Math/ForwardDecls.h"
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
   * @ingroup RodinCTAD
   * @brief CTAD for GMRES
   */
  template <class LinearSystem>
  GMRES(Variational::ProblemBase<LinearSystem>&) -> GMRES<LinearSystem>;

  /**
   * @ingroup GMRESSpecializations
   * @brief GMRES for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class GMRES<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
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

      GMRES(ProblemBaseType& pb)
        : Parent(pb)
      {}

      GMRES(const GMRES& other)
        : Parent(other)
      {}

      GMRES(GMRES&& other)
        : Parent(std::move(other))
      {}

      ~GMRES() = default;

      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      GMRES& setTolerance(const Real& tol)
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
}

#endif



