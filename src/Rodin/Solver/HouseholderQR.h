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
   * @ingroup RodinCTAD
   * @brief CTAD for HouseholderQR
   */
  HouseholderQR() -> HouseholderQR<Math::Matrix<Scalar>, Math::Vector<Scalar>>;

  /**
   * @defgroup HouseholderQRSpecializations HouseholderQR Template Specializations
   * @brief Template specializations of the HouseholderQR class.
   * @see HouseholderQR
   */

  /**
   * @ingroup HouseholderQRSpecializations
   * @brief A direct sparse HouseholderQR Cholesky factorizations without square root
   * for use with Math::SparseMatrix<Scalar> and Math::Vector<Scalar>.
   */
  template <>
  class HouseholderQR<Math::Matrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::Matrix<Scalar>, Math::Vector<Scalar>>
  {
    public:
      /// Type of linear operator
      using OperatorType = Math::Matrix<Scalar>;

      /// Type of vector
      using VectorType = Math::Vector<Scalar>;

      /**
       * @brief Constructs the HouseholderQR object with default parameters.
       */
      HouseholderQR() = default;

      HouseholderQR(const HouseholderQR& other)
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
      Eigen::HouseholderQR<Math::Matrix<Scalar>> m_solver;
  };

}

#endif




