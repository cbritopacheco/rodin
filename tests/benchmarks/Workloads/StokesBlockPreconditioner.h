/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_WORKLOADS_STOKESBLOCKPRECONDITIONER_H
#define RODIN_WORKLOADS_STOKESBLOCKPRECONDITIONER_H

#include <cassert>
#include <cmath>
#include <memory>
#include <utility>
#include <Eigen/SparseCholesky>
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Tests::Benchmarks
{
  /**
   * @brief Fixed linear block preconditioner for the unit-viscosity Stokes workload.
   *
   * For velocity, pressure and mean blocks, incomplete Cholesky
   * approximates the velocity inverse and the pressure mass matrix approximates
   * @f$-BK^{-1}C@f$. The bordered mean constraint is solved explicitly using
   * pressure-mass solves; no dense Stokes Schur matrix or pressure pin is used.
   * Each application performs a full approximate block factorization.
   *
   * The immutable pressure mass matrix is owned and shared between configuration
   * copies. Factors and work arrays belong to each instance. Copying preserves
   * configuration only; compute must succeed before solve. This is a benchmark
   * fixture for the scalar unit-viscosity forms, not a general Stokes solver.
   */
  class StokesBlockPreconditioner
  {
    public:
      /// @brief Scalar type required by Eigen's preconditioner protocol.
      using Scalar = Real;

      /// @brief Sparse index type of the native operator.
      using StorageIndex = Math::SparseMatrix<Real>::StorageIndex;

      /// @brief Native sparse operator type.
      using Matrix = Math::SparseMatrix<Real>;

      /// @brief Native coefficient vector type.
      using Vector = Math::Vector<Real>;
      enum
      {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic
      };

      /// @brief Constructs an unconfigured preconditioner.
      StokesBlockPreconditioner() = default;

      /// @brief Copies owned configuration into an uncomputed instance.
      StokesBlockPreconditioner(const StokesBlockPreconditioner& other)
        : m_mass(other.m_mass),
          m_velocitySize(other.m_velocitySize)
      {}

      /// @brief Copies configuration and invalidates any previous factorization.
      StokesBlockPreconditioner& operator=(const StokesBlockPreconditioner& other)
      {
        m_mass = other.m_mass;
        m_velocitySize = other.m_velocitySize;
        m_status = Eigen::InvalidInput;
        return *this;
      }

      /**
       * @brief Owns the pressure mass matrix and records the explicit block split.
       * @param mass Pressure mass matrix, moved into immutable owned storage.
       * @param velocitySize First pressure offset in the assembled system.
       */
      void configure(Matrix mass, size_t velocitySize)
      {
        m_mass = std::make_shared<const Matrix>(std::move(mass));
        m_velocitySize = velocitySize;
        m_status = Eigen::InvalidInput;
      }

      /// @brief Returns the expected total number of unknowns, including the mean.
      Eigen::Index rows() const
      {
        return m_mass ? Eigen::Index(m_velocitySize) + m_mass->rows() + 1 : 0;
      }

      /// @brief Returns the expected total number of unknowns.
      Eigen::Index cols() const
      {
        return rows();
      }

      /**
       * @brief Builds independent velocity and pressure factors from the full system.
       * @param matrix Assembled velocity, pressure and mean operator in that order.
       * @returns This preconditioner; inspect info before applying it.
       */
      template <class Operator>
      StokesBlockPreconditioner& compute(const Operator& matrix)
      {
        m_status = Eigen::InvalidInput;
        if (!m_mass || m_velocitySize == 0 || m_mass->rows() == 0 ||
          m_mass->rows() != m_mass->cols() || matrix.rows() != rows() ||
          matrix.cols() != rows())
          return *this;
        const Eigen::Index nv = m_velocitySize, np = m_mass->rows();
        Matrix velocity = matrix.topLeftCorner(nv, nv);
        velocity.prune([](StorageIndex, StorageIndex, Real value) { return value != 0; });
        m_velocity.compute(velocity);
        m_pressure.compute(*m_mass);
        m_status = Eigen::NumericalIssue;
        if (m_velocity.info() != Eigen::Success || m_pressure.info() != Eigen::Success ||
          (m_pressure.vectorD().array() <= 0).any())
          return *this;
        m_divergence = matrix.block(nv, 0, np, nv);
        m_momentum = matrix.block(0, nv, nv, np);
        m_mean = matrix.block(nv, nv + np, np, 1);
        m_inverseMean = m_pressure.solve(m_mean);
        m_meanDenominator = m_mean.dot(m_inverseMean);
        if (std::isfinite(m_meanDenominator) && m_meanDenominator > 0)
          m_status = Eigen::Success;
        return *this;
      }

      /**
       * @brief Applies the fixed approximate inverse, including the mean equation.
       * @param rhs Right-hand side of all three blocks.
       * @param result Destination coefficient vector.
       */
      template <class Rhs, class Destination>
      void _solve_impl(const Rhs& rhs, Destination& result) const
      {
        assert(m_status == Eigen::Success);
        const Eigen::Index nv = m_velocitySize, np = m_mass->rows();
        const Vector velocity = m_velocity.solve(rhs.head(nv));
        const Vector pressure =
          m_pressure.solve(rhs.segment(nv, np) - m_divergence * velocity);
        const Real mean = (m_mean.dot(pressure) - rhs(nv + np)) / m_meanDenominator;
        const Vector correctedPressure = pressure - m_inverseMean * mean;
        result.head(nv) = velocity - m_velocity.solve(m_momentum * correctedPressure);
        result.segment(nv, np) = correctedPressure;
        result(nv + np) = mean;
      }

      /// @brief Returns an Eigen solve expression after successful setup.
      template <class Rhs>
      const Eigen::Solve<StokesBlockPreconditioner, Rhs> solve(
        const Eigen::MatrixBase<Rhs>& rhs) const
      {
        assert(m_status == Eigen::Success && rhs.rows() == rows());
        return Eigen::Solve<StokesBlockPreconditioner, Rhs>(*this, rhs.derived());
      }

      /// @brief Reports setup success, missing configuration or a failed factor.
      Eigen::ComputationInfo info() const
      {
        return m_status;
      }

    private:
      std::shared_ptr<const Matrix> m_mass;
      size_t m_velocitySize = 0;
      Eigen::IncompleteCholesky<Real> m_velocity;
      Eigen::SimplicialLDLT<Matrix> m_pressure;
      Matrix m_divergence, m_momentum;
      Vector m_mean, m_inverseMean;
      Real m_meanDenominator = 0;
      Eigen::ComputationInfo m_status = Eigen::InvalidInput;
  };
}
#endif
