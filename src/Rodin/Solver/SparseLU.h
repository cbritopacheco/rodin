/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SPARSELU_H
#define RODIN_SOLVER_SPARSELU_H

#include <Eigen/SparseLU>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup SparseLUSpecializations SparseLU Template Specializations
   * @brief Template specializations of the SparseLU class.
   * @see SparseLU
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SparseLU
   */
  template <class LinearSystem>
  SparseLU(Variational::ProblemBase<LinearSystem>&) -> SparseLU<LinearSystem>;

  /**
   * @ingroup SparseLUSpecializations
   * @brief Sparse supernodal LU factorization for general matrices for use
   * with Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <class Scalar>
  class SparseLU<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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

      SparseLU(ProblemBaseType& pb)
        : Parent(pb)
      {}

      SparseLU(const SparseLU& other)
        : Parent(other)
      {}

      SparseLU(SparseLU&& other)
        : Parent(std::move(other))
      {}

      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      SparseLU* copy() const noexcept override
      {
        return new SparseLU(*this);
      }

    private:
      Eigen::SparseLU<OperatorType> m_solver;
  };
}

#endif



