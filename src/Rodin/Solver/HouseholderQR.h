/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_HOUSEHOLDERQR_H
#define RODIN_SOLVER_HOUSEHOLDERQR_H

#include <Eigen/Dense>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup HouseholderQRSpecializations HouseholderQR Template Specializations
   * @brief Template specializations of the HouseholderQR class.
   * @see HouseholderQR
   */

  /**
   * @ingroup HouseholderQRSpecializations
   * @brief A direct sparse HouseholderQR Cholesky factorizations without square root
   * for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class HouseholderQR<Math::Matrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::Matrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::Matrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      HouseholderQR(ProblemType& pb)
        : Parent(pb)
      {}

      HouseholderQR(const HouseholderQR& other)
        : Parent(other)
      {}

      HouseholderQR(HouseholderQR&& other)
        : Parent(std::move(other))
      {}

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        x = m_solver.compute(A).solve(b);
      }

      inline
      HouseholderQR* copy() const noexcept override
      {
        return new HouseholderQR(*this);
      }

    private:
      Eigen::HouseholderQR<OperatorType> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for HouseholderQR
   */
  template <class Scalar>
  HouseholderQR(Variational::ProblemBase<Math::Matrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> HouseholderQR<Math::Matrix<Scalar>, Math::Vector<Scalar>>;
}

#endif




