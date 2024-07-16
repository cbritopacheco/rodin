/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SIMPLICIALLLT_H
#define RODIN_SOLVER_SIMPLICIALLLT_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup SimplicialLLTSpecializations SimplicialLLT Template Specializations
   * @brief Template specializations of the SimplicialLLT class.
   * @see SimplicialLLT
   */

  /**
   * @ingroup SimplicialLLTSpecializations
   * @brief A direct sparse LLT Cholesky factorizations for use with
   * Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <class Scalar>
  class SimplicialLLT<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      SimplicialLLT(ProblemType& pb)
        : Parent(pb)
      {}

      SimplicialLLT(const SimplicialLLT& other)
        : Parent(other)
      {}

      SimplicialLLT(SimplicialLLT&& other)
        : Parent(std::move(other))
      {}

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        x = m_solver.compute(A).solve(b);
      }

      inline
      SimplicialLLT* copy() const noexcept override
      {
        return new SimplicialLLT(*this);
      }

    private:
      Eigen::SimplicialLLT<OperatorType> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SimplicialLLT
   */
  template <class Scalar>
  SimplicialLLT(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> SimplicialLLT<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;
}

#endif



