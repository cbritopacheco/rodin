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
   * @brief CTAD (Class Template Argument Deduction) guide for SparseLU
   */
  template <class LinearSystem>
  SparseLU(Variational::ProblemBase<LinearSystem>&) -> SparseLU<LinearSystem>;

  /**
   * @ingroup SparseLUSpecializations
   * @brief Sparse supernodal LU factorization solver for general matrices.
   *
   * This class provides a sparse direct solver using LU decomposition with
   * supernodal techniques. It performs the factorization @f$ PA = LU @f$
   * where @f$ P @f$ is a permutation matrix, @f$ L @f$ is lower triangular,
   * and @f$ U @f$ is upper triangular.
   *
   * The solver is suitable for general (non-symmetric) sparse matrices and
   * provides robust solutions for systems arising from finite element
   * discretizations.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   *
   * **Example usage:**
   * @code{.cpp}
   * Problem problem(u, v);
   * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
   * 
   * Solver::SparseLU solver(problem);
   * solver.solve();
   * @endcode
   */
  template <class Scalar>
  class SparseLU<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
    : public SolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      /// Type of scalar values
      using ScalarType = Scalar;

      /// Type of solution and right-hand side vectors
      using VectorType = Math::Vector<ScalarType>;

      /// Type of system matrix (operator)
      using OperatorType = Math::SparseMatrix<ScalarType>;

      /// Type of linear system being solved
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      /// Type of problem being solved
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      /// Parent class type
      using Parent = SolverBase<LinearSystemType>;

      using Parent::solve;

      /**
       * @brief Constructs a SparseLU solver for the given problem.
       * @param pb Reference to the problem to solve
       */
      SparseLU(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copy constructor.
       * @param other Solver to copy from
       */
      SparseLU(const SparseLU& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Solver to move from
       */
      SparseLU(SparseLU&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Solves the linear system using sparse LU factorization.
       * @param[in,out] axb The linear system to solve. The solution is stored
       *                    in the system's solution vector.
       *
       * This method computes the LU factorization of the system matrix and
       * uses it to solve for the solution vector.
       */
      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new SparseLU instance
       */
      SparseLU* copy() const noexcept override
      {
        return new SparseLU(*this);
      }

    private:
      /// Underlying Eigen solver implementation
      Eigen::SparseLU<OperatorType> m_solver;
  };
}

#endif



