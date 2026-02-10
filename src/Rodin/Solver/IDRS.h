/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file IDRS.h
 * @brief IDRS(s) solver for (generally) non-symmetric linear systems.
 *
 * This header provides the IDRS(s) (Induced Dimension Reduction) solver class,
 * an iterative method for solving linear systems:
 * @f[
 *   Ax = b
 * @f]
 * for general (possibly non-symmetric) matrices.
 *
 * IDR(s) is typically competitive with BiCGSTAB/GMRES for certain nonsymmetric
 * problems, with short recurrences and modest memory footprint.
 *
 * ## Notes
 * - This implementation uses Eigen's unsupported IDRS.
 *
 * ## Usage Example
 * ```cpp
 * Solver::IDRS solver(problem);
 * solver.setTolerance(1e-10).setMaxIterations(1000).setS(8).solve();
 * ```
 */
#ifndef RODIN_SOLVER_IDRS_H
#define RODIN_SOLVER_IDRS_H

#include <Eigen/Core>
#include <unsupported/Eigen/IterativeSolvers>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Matrix.h"

#include "ForwardDecls.h"
#include "Rodin/Types.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup IDRSSpecializations IDRS Template Specializations
   * @brief Template specializations of the IDRS class.
   * @see IDRS
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for IDRS
   */
  template <class LinearSystem>
  IDRS(Variational::ProblemBase<LinearSystem>&) -> IDRS<LinearSystem>;

  /**
   * @ingroup IDRSSpecializations
   * @brief IDRS(s) solver for sparse systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class IDRS<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>> final
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

      IDRS(ProblemBaseType& pb)
        : Parent(pb)
      {}

      IDRS(const IDRS& other)
        : Parent(other)
      {}

      IDRS(IDRS&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      ~IDRS() = default;

      IDRS& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      IDRS& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      /**
       * @brief Sets the IDR(s) parameter s (dimension of the shadow space).
       *
       * Typical values: 4, 8, 16. Larger s can improve robustness but increases
       * work per iteration.
       */
      IDRS& setS(size_t s)
      {
        m_solver.setS(s);
        return *this;
      }

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      IDRS* copy() const noexcept override
      {
        return new IDRS(*this);
      }

    private:
      Eigen::IDRS<OperatorType> m_solver;
  };

  /**
   * @ingroup IDRSSpecializations
   * @brief IDRS(s) solver for dense systems.
   *
   * @tparam Scalar The scalar type (e.g., Real, Complex)
   */
  template <class Scalar>
  class IDRS<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>> final
    : public SolverBase<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;
      using VectorType = Math::Vector<ScalarType>;
      using OperatorType = Math::Matrix<ScalarType>;
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      using ProblemType = Variational::ProblemBase<LinearSystemType>;
      using Parent = SolverBase<LinearSystemType>;

      using Parent::solve;

      IDRS(ProblemType& pb)
        : Parent(pb)
      {}

      IDRS(const IDRS& other)
        : Parent(other),
          m_solver(other.m_solver)
      {}

      IDRS(IDRS&& other)
        : Parent(std::move(other)),
          m_solver(std::move(other.m_solver))
      {}

      ~IDRS() = default;

      IDRS& setTolerance(const Real& tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      IDRS& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      IDRS& setS(size_t s)
      {
        m_solver.setS(s);
        return *this;
      }

      void solve(LinearSystemType& axb) override
      {
        axb.getSolution() = m_solver.compute(axb.getOperator()).solve(axb.getVector());
      }

      Boolean success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      IDRS* copy() const noexcept override
      {
        return new IDRS(*this);
      }

    private:
      Eigen::IDRS<OperatorType> m_solver;
  };
}

#endif
