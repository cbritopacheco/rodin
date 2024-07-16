/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SOLVER_H
#define RODIN_SOLVER_SOLVER_H

#include "Rodin/Configure.h"
#include "Rodin/Copyable.h"
#include "Rodin/Math.h"

#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

namespace Rodin::Solver
{
  template <class Operator, class Vector, class Scalar>
  class SolverBase : public Copyable
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Vector;

      using OperatorType = Operator;

      using Parent = Copyable;

      /**
       * @brief Default virtual destructor.
       */
      virtual ~SolverBase() = default;

      SolverBase(Variational::ProblemBase<OperatorType, VectorType, ScalarType>& pb)
        : m_pb(pb)
      {}

      SolverBase(const SolverBase& other)
        : Parent(other),
          m_pb(other.m_pb)
      {}

      SolverBase(SolverBase&& other)
        : Parent(other),
          m_pb(std::move(other.m_pb))
      {}

      void solve()
      {
        m_pb.get().solve(*this);
      }

      /**
       * @brief Solves the linear algebra system.
       *
       * Solves the following system:
       * @f[
       *  Ax = b ,
       * @f]
       * where @f$ A @f$ has type @f$ \text{OperatorType} @f$, the solution
       * @f$ x @f$ has type @f$ \text{VectorType} @f$, and the right hand side
       * @f$ b @f$ has type @f$ \text{VectorType} @f$.
       *
       * @param[in,out] A Left hand side operator
       * @param[in,out] x Input for the initial guess, and output for the
       * solution
       * @param[in,out] b Right hand side vector
       */
      virtual void solve(OperatorType& A, VectorType& x, VectorType& b) = 0;

    private:
      std::reference_wrapper<Variational::ProblemBase<OperatorType, VectorType, ScalarType>> m_pb;
  };
}

#endif
