/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SPQR_H
#define RODIN_SOLVER_SPQR_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_SPQR

#include <optional>
#include <functional>

#include <Eigen/SPQRSupport>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup SPQRSpecializations SPQR Template Specializations
   * @brief Template specializations of the SPQR class.
   * @see SPQR
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SPQR
   */
  template <class LinearSystem>
  SPQR(Variational::ProblemBase<LinearSystem>&) -> SPQR<LinearSystem>;

  /**
   * @ingroup SPQRSpecializations
   * @brief SPQR for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class SPQR<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
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

      SPQR(ProblemBaseType& pb)
        : Parent(pb)
      {}

      SPQR(const SPQR& other)
        : Parent(other)
      {}

      SPQR(SPQR&& other)
        : Parent(std::move(other))
      {}

      ~SPQR() = default;

      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      SPQR& setTolerance(Real tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      SPQR& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      SPQR* copy() const noexcept override
      {
        return new SPQR(*this);
      }

    private:
      Eigen::SPQR<OperatorType> m_solver;
  };
}

#endif
#endif




