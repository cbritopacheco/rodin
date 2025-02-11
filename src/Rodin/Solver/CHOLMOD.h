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
  template <class Scalar>
  class SupernodalLLT<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      SupernodalLLT(ProblemType& pb)
        : Parent(pb)
      {}

      SupernodalLLT(const SupernodalLLT& other)
        : Parent(other)
      {}

      SupernodalLLT(SupernodalLLT&& other)
        : Parent(std::move(other))
      {}

      ~SupernodalLLT() = default;

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

      bool success() const
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

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SupernodalLLT
   */
  template <class Scalar>
  SupernodalLLT(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> SupernodalLLT<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;
}

#endif // #ifdef RODIN_USE_CHOLMOD
#endif



