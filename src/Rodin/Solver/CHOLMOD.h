/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_CHOLMOD_H
#define RODIN_SOLVER_CHOLMOD_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_CHOLMOD

#include <optional>
#include <functional>

#include <Eigen/CholmodSupport>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver::CHOLMOD
{
  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SupernodalLLT
   */
  template <class LinearSystemType>
  SupernodalLLT(Variational::ProblemBase<LinearSystemType>&) -> SupernodalLLT<LinearSystemType>;

  template <class Scalar>
  class SupernodalLLT<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
    : public SolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      using Parent = SolverBase<LinearSystemType>;

      using Parent::solve;

      SupernodalLLT(ProblemBaseType& pb)
        : Parent(pb)
      {}

      SupernodalLLT(const SupernodalLLT& other)
        : Parent(other)
      {}

      SupernodalLLT(SupernodalLLT&& other)
        : Parent(std::move(other))
      {}

      ~SupernodalLLT() = default;

      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      Real success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      SupernodalLLT* copy() const noexcept override
      {
        return new SupernodalLLT(*this);
      }

    private:
      Eigen::CholmodSupernodalLLT<OperatorType> m_solver;
  };
}

#endif // #ifdef RODIN_USE_CHOLMOD
#endif



