/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SPARSEQR_H
#define RODIN_SOLVER_SPARSEQR_H

#include <Eigen/SparseQR>

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
   * @ingroup SparseQRSpecializations
   * @brief Sparse left-looking QR factorization with numerical column pivoting
   * for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class SparseQR<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      SparseQR(ProblemType& pb)
        : Parent(pb)
      {}

      SparseQR(const SparseQR& other)
        : Parent(other)
      {}

      SparseQR(SparseQR&& other)
        : Parent(std::move(other))
      {}

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

      inline
      SparseQR* copy() const noexcept override
      {
        return new SparseQR(*this);
      }

    private:
      Eigen::SparseQR<OperatorType, Eigen::COLAMDOrdering<int>> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SparseQR
   */
  template <class Scalar>
  SparseQR(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> SparseQR<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;
}

#endif



