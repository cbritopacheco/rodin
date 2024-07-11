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
   * @ingroup SparseLUSpecializations
   * @brief Sparse supernodal LU factorization for general matrices for use
   * with Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <class Scalar>
  class SparseLU<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      SparseLU(ProblemType& pb)
        : Parent(pb)
      {}

      SparseLU(const SparseLU& other)
        : Parent(other)
      {}

      SparseLU(SparseLU&& other)
        : Parent(std::move(other))
      {}

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

      SparseLU* copy() const noexcept override
      {
        return new SparseLU(*this);
      }

    private:
      Eigen::SparseLU<OperatorType> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SparseLU
   */
  template <class Scalar>
  SparseLU(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> SparseLU<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;
}

#endif



