/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PartialPivLU.h
 * @brief Dense LU factorization solver with partial pivoting.
 *
 * This header provides the PartialPivLU solver class, which implements
 * dense LU factorization for general square matrices using Eigen's
 * PartialPivLU solver.
 *
 * ## Algorithm
 * The solver computes
 * @f[
 *   A = P^{-1} L U
 * @f]
 * or equivalently
 * @f[
 *   P A = L U
 * @f]
 * where @f$ P @f$ is a permutation matrix, @f$ L @f$ is unit lower triangular,
 * and @f$ U @f$ is upper triangular.
 *
 * ## Applicability
 * - General dense square matrices
 * - Non-symmetric systems
 * - Newton systems with dense Jacobians
 *
 * ## Notes
 * - This solver does not assume symmetry.
 * - Eigen's PartialPivLU assumes the matrix is square and invertible.
 * - If rank deficiency is possible and must be handled robustly, FullPivLU is safer.
 *
 * @see Eigen::PartialPivLU
 */
#ifndef RODIN_SOLVER_PARTIALPIVLU_H
#define RODIN_SOLVER_PARTIALPIVLU_H

#include <Eigen/LU>

#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"

#include "ForwardDecls.h"
#include "LinearSolver.h"

namespace Rodin::FormLanguage
{
  template <class LinearSystem>
  struct Traits<Solver::PartialPivLU<LinearSystem>>
  {
    using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @ingroup RodinCTAD
   * @brief CTAD guide for PartialPivLU.
   */
  template <class LinearSystem>
  PartialPivLU(Variational::ProblemBase<LinearSystem>&) -> PartialPivLU<LinearSystem>;

  /**
   * @brief Dense LU factorization solver with partial pivoting for general matrices.
   *
   * @tparam Scalar Scalar type.
   */
  template <class Scalar>
  class PartialPivLU<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public LinearSolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;
      using VectorType = Math::Vector<ScalarType>;
      using OperatorType = Math::Matrix<ScalarType>;
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      PartialPivLU(ProblemBaseType& pb)
        : Parent(pb)
      {}

      PartialPivLU(const PartialPivLU& other)
        : Parent(other)
      {}

      PartialPivLU(PartialPivLU&& other)
        : Parent(std::move(other))
      {}

      void solve(LinearSystemType& axb) override
      {
        const auto& A = axb.getOperator();
        const auto& b = axb.getVector();

        assert(A.rows() == A.cols());
        assert(A.rows() == b.size());

        m_solver.compute(A);

        // Eigen documents that info() is always Success for PartialPivLU,
        // so this does not certify invertibility. It only reflects API consistency.
        axb.getSolution() = m_solver.solve(b);
      }

      PartialPivLU* copy() const noexcept override
      {
        return new PartialPivLU(*this);
      }

      /**
       * @brief Returns the underlying Eigen solver.
       *
       * Useful for diagnostics such as permutationP(), matrixLU(), rcond().
       */
      const Eigen::PartialPivLU<OperatorType>& getSolver() const noexcept
      {
        return m_solver;
      }

    private:
      Eigen::PartialPivLU<OperatorType> m_solver;
  };
}

#endif
