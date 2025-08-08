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
#include "Rodin/Math/ForwardDecls.h"
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
   * @ingroup RodinCTAD
   * @brief CTAD for DGMRES
   */
  template <class LinearSystem>
  DGMRES(Variational::ProblemBase<LinearSystem>&) -> DGMRES<LinearSystem>;

  /**
   * @ingroup DGMRESSpecializations
   * @brief DGMRES for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class DGMRES<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
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

      DGMRES(ProblemBaseType& pb)
        : Parent(pb)
      {}

      DGMRES(const DGMRES& other)
        : Parent(other)
      {}

      DGMRES(DGMRES&& other)
        : Parent(std::move(other))
      {}

      ~DGMRES() = default;

      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      DGMRES& setTolerance(const Real& tol)
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
}

#endif




