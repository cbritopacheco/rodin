/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_HOUSEHOLDERQR_H
#define RODIN_SOLVER_HOUSEHOLDERQR_H

#include <Eigen/Dense>

#include "Rodin/Math/ForwardDecls.h"
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
   * @ingroup RodinCTAD
   * @brief CTAD for HouseholderQR
   */
  template <class LinearSystem>
  HouseholderQR(Variational::ProblemBase<LinearSystem>&) -> HouseholderQR<LinearSystem>;

  /**
   * @ingroup HouseholderQRSpecializations
   * @brief A direct sparse HouseholderQR Cholesky factorizations without square root
   * for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class HouseholderQR<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public SolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::Matrix<ScalarType>;

      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      using Parent = SolverBase<LinearSystemType>;

      using Parent::solve;

      HouseholderQR(ProblemBaseType& pb)
        : Parent(pb)
      {}

      HouseholderQR(const HouseholderQR& other)
        : Parent(other)
      {}

      HouseholderQR(HouseholderQR&& other)
        : Parent(std::move(other))
      {}

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      HouseholderQR* copy() const noexcept override
      {
        return new HouseholderQR(*this);
      }

    private:
      Eigen::HouseholderQR<OperatorType> m_solver;
  };
}

#endif




