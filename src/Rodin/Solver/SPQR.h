/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SPQR_H
#define RODIN_SOLVER_SPQR_H

#ifdef RODIN_USE_SPQR

#include <optional>
#include <functional>

#include <Eigen/SPQRSupport>

#include "Rodin/Configure.h"
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
   * @ingroup SPQRSpecializations
   * @brief SPQR for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class SPQR<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      SPQR(ProblemType& pb)
        : Parent(pb)
      {}

      SPQR(const SPQR& other)
        : Parent(other)
      {}

      SPQR(SPQR&& other)
        : Parent(std::move(other))
      {}

      ~SPQR() = default;

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
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

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SPQR
   */
  template <class Scalar>
  SPQR(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> SPQR<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;
}

#endif
#endif




