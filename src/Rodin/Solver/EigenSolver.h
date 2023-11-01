/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EIGENSOLVER_H
#define RODIN_EIGENSOLVER_H

#include <functional>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @ingroup RodinCTAD
   * @brief CTAD for EigenSolver
   */
  template <class EigenSolverType>
  EigenSolver(std::reference_wrapper<EigenSolverType>)
    -> EigenSolver<EigenSolverType, Math::SparseMatrix, Math::Vector>;

  /**
   * @defgroup EigenSolverSpecializations EigenSolver Template Specializations
   * @brief Template specializations of the EigenSolver class.
   * @see EigenSolver
   */

  template <class EigenSolverType, class Operator, class Vector>
  class EigenSolver final
    : public SolverBase<Operator, Vector>
  {
    public:
      /// Type of linear operator
      using OperatorType = Operator;

      /// Type of vector
      using VectorType = Vector;

      explicit
      EigenSolver(std::reference_wrapper<EigenSolverType> solver)
        : m_solver(solver)
      {}

      ~EigenSolver() = default;

      /**
       * @brief Gets the reference to the underlying solver instance.
       */
      inline
      EigenSolverType& getSolver()
      {
        return m_solver.get();
      }

      /**
       * @brief Gets the constant reference to underlying solver instance.
       */
      inline
      const EigenSolverType& getSolver() const
      {
        return m_solver.get();
      }

      /**
       * @brief Solves the linear algebra system.
       *
       * This method has the following definition:
       * @code{cpp}
       *   inline
       *   void solve(OperatorType& A, VectorType& x, VectorType& b) override
       *   {
       *     x = m_solver.get().compute(A).solve(b);
       *   }
       * @endcode
       */
      inline
      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        x = getSolver().compute(A).solve(b);
      }

    private:
      std::reference_wrapper<EigenSolverType> m_solver;
  };
}

#endif
