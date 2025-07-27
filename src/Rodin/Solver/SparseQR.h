/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SPARSEQR_H
#define RODIN_SOLVER_SPARSEQR_H

#include <Eigen/SparseQR>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup SparseQRSpecializations SparseQR Template Specializations
   * @brief Template specializations of the SparseQR class.
   * @see SparseQR
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SparseQR
   */
  template <class LinearSystem>
  SparseQR(Variational::ProblemBase<LinearSystem>&) -> SparseQR<LinearSystem>;

  /**
   * @ingroup SparseQRSpecializations
   * @brief Sparse left-looking QR factorization with numerical column pivoting
   * for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class SparseQR<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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

      SparseQR(ProblemBaseType& pb)
        : Parent(pb)
      {}

      SparseQR(const SparseQR& other)
        : Parent(other)
      {}

      SparseQR(SparseQR&& other)
        : Parent(std::move(other))
      {}

      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      SparseQR* copy() const noexcept override
      {
        return new SparseQR(*this);
      }

    private:
      Eigen::SparseQR<OperatorType, Eigen::COLAMDOrdering<int>> m_solver;
  };
}

#endif



