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
  /**
   * @defgroup RodinSolvers Linear System Solvers
   * @brief Solver classes for linear systems of equations
   *
   * This module provides various solver implementations for solving linear
   * systems @f$ Ax = b @f$ that arise from finite element discretizations.
   */

  /**
   * @ingroup RodinSolvers
   * @brief Base class for linear system solvers.
   *
   * This abstract base class defines the interface for solving linear systems
   * of equations. Derived classes implement specific solver algorithms such as
   * direct factorization methods or iterative methods.
   *
   * @tparam LinearSystem Type of linear system to solve
   */
  template <class LinearSystem>
  class SolverBase : public Copyable
  {
    public:
      /// Scalar type used in the linear system
      using ScalarType =
        typename FormLanguage::Traits<LinearSystem>::ScalarType;

      /// Vector type for solution and right-hand side
      using VectorType =
        typename FormLanguage::Traits<LinearSystem>::VectorType;

      /// Operator (matrix) type for the linear system
      using OperatorType =
        typename FormLanguage::Traits<LinearSystem>::OperatorType;

      /// Base problem type 
      using ProblemBaseType =
        Variational::ProblemBase<LinearSystem>;

      /// Parent class type
      using Parent = Copyable;

      /**
       * @brief Default virtual destructor.
       */
      virtual ~SolverBase() = default;

      /**
       * @brief Constructs a solver for the given problem.
       * @param pb Reference to the problem to solve
       */
      SolverBase(ProblemBaseType& pb)
        : m_pb(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      SolverBase(const SolverBase& other)
        : Parent(other),
          m_pb(other.m_pb)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from  
       */
      SolverBase(SolverBase&& other)
        : Parent(other),
          m_pb(std::move(other.m_pb))
      {}

      /**
       * @brief Solves the associated problem.
       *
       * This method delegates solving to the problem's solve method,
       * passing this solver as the argument.
       */
      void solve()
      {
        m_pb.get().solve(*this);
      }

      /**
       * @brief Solves the linear system.
       *
       * Solves the linear system @f$ Ax = b @f$ where @f$ A @f$ is the matrix
       * operator, @f$ x @f$ is the solution vector, and @f$ b @f$ is the 
       * right-hand side vector.
       *
       * @param[in,out] system The linear system to solve. The solution is 
       *                       stored in the system's solution vector.
       */
      virtual void solve(LinearSystem& system) = 0;

    private:
      /// Reference to the problem being solved
      std::reference_wrapper<ProblemBaseType> m_pb;
  };
}

#endif
