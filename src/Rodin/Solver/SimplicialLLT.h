/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SIMPLICIALLLT_H
#define RODIN_SOLVER_SIMPLICIALLLT_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/ForwardDecls.h"
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
   * @ingroup RodinCTAD
   * @brief CTAD for SimplicialLLT
   */
  template <class LinearSystem>
  SimplicialLLT(Variational::ProblemBase<LinearSystem>&) -> SimplicialLLT<LinearSystem>;

  /**
   * @ingroup SimplicialLLTSpecializations
   * @brief A direct sparse LLT Cholesky factorizations for use with
   * Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <class Scalar>
  class SimplicialLLT<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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

      SimplicialLLT(ProblemBaseType& pb)
        : Parent(pb)
      {}

      SimplicialLLT(const SimplicialLLT& other)
        : Parent(other)
      {}

      SimplicialLLT(SimplicialLLT&& other)
        : Parent(std::move(other))
      {}

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      SimplicialLLT* copy() const noexcept override
      {
        return new SimplicialLLT(*this);
      }

    private:
      Eigen::SimplicialLLT<OperatorType> m_solver;
  };
}

#endif



