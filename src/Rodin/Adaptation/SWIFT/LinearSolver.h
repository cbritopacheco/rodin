/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_LINEARSOLVER_H
#define RODIN_ADAPTATION_SWIFT_LINEARSOLVER_H
#include "Rodin/Alert/MemberFunctionException.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <limits>
#include <Eigen/IterativeLinearSolvers>
#include "Rodin/Solver/LinearSolver.h"
#include "Rodin/Solver/MUMPS.h"
#include "Rodin/Solver/SparseLU.h"
#include "Parameters.h"
#include "Report.h"

namespace Rodin::Adaptation::SWIFT
{
  template <class Scalar>
  class LowRankUpdate;
}

namespace Eigen::internal
{
  template <class Scalar>
  struct traits<Rodin::Adaptation::SWIFT::LowRankUpdate<Scalar>>
    : traits<Rodin::Math::SparseMatrix<Scalar>>
  {};
}

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Non-owning operator @f$M=A-UU^*@f$ with a sparse core.
   *
   * ## Architecture
   *
   * Matrix-vector products retain the low-rank subtraction without forming a
   * dense matrix. Direct solvers can use @ref getAugmented, whose auxiliary
   * equation eliminates to the same operator. Neither representation adds a
   * gauge or regularization. The core and factor must outlive this view and
   * must not change during a solve.
   *
   * @tparam ScalarType Scalar of the core and factor.
   */
  template <class ScalarType>
  class LowRankUpdate final : public Eigen::EigenBase<LowRankUpdate<ScalarType>>
  {
    public:
      /// @brief Scalar required by Eigen's matrix-free solver interface.
      using Scalar = ScalarType;

      /// @brief Real scalar associated with the operator.
      using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

      /// @brief Sparse index type.
      using StorageIndex = typename Math::SparseMatrix<Scalar>::StorageIndex;

      /// @brief Dynamic matrix dimensions for Eigen expressions.
      enum
      {
        ColsAtCompileTime = Eigen::Dynamic, ///< Runtime column count.
        MaxColsAtCompileTime = Eigen::Dynamic, ///< Unbounded column count.
        IsRowMajor = false ///< Column-major operator interface.
      };

      /**
       * @brief Binds the sparse core and low-rank factor.
       * @param core Square sparse matrix.
       * @param factor Factor with one row per core row, or an empty matrix.
       */
      LowRankUpdate(
        const Math::SparseMatrix<Scalar>& core, const Math::Matrix<Scalar>& factor)
        : m_core(core),
          m_factor(factor)
      {
        if (core.rows() != core.cols() ||
          (factor.size() != 0 && factor.rows() != core.rows()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "A low-rank update requires a square core and matching factor rows."
            << Alert::Raise;
        }
      }

      /**
       * @brief Returns the number of rows.
       * @returns Core row count.
       */
      Eigen::Index rows() const
      {
        return m_core.rows();
      }

      /**
       * @brief Returns the number of columns.
       * @returns Core column count.
       */
      Eigen::Index cols() const
      {
        return m_core.cols();
      }

      /**
       * @brief Applies @f$Ax-U(U^*x)@f$ without dense assembly.
       * @tparam Derived Input vector expression.
       * @param x Vector to multiply.
       * @returns Operator image.
       */
      template <class Derived>
      Math::Vector<Scalar> operator*(const Eigen::MatrixBase<Derived>& x) const
      {
        Math::Vector<Scalar> image = m_core * x;
        if (m_factor.size() != 0)
          image -= m_factor * (m_factor.adjoint() * x);
        return image;
      }

      /**
       * @brief Returns the diagonal of the complete operator.
       * @returns Core diagonal minus squared factor-row norms.
       */
      Math::Vector<Scalar> diagonal() const
      {
        Math::Vector<Scalar> result = m_core.diagonal();
        if (m_factor.size() != 0)
          result -= m_factor.cwiseAbs2().rowwise().sum().template cast<Scalar>();
        return result;
      }

      /**
       * @brief Forms @f$\left(\begin{smallmatrix}A&U\\U^*&I\end{smallmatrix}\right)@f$.
       * @returns Sparse augmented matrix, or the core for an empty factor.
       */
      Math::SparseMatrix<Scalar> getAugmented() const
      {
        if (m_factor.size() == 0)
          return m_core;
        const auto n = rows(), rank = m_factor.cols();
        Math::SparseMatrix<Scalar> result(n + rank, n + rank);
        result.reserve(m_core.nonZeros() + 2 * n * rank + rank);
        // Columns are already ordered. Write compressed storage directly,
        // preserving structural zeros without a second triplet-sized buffer.
        for (Eigen::Index column = 0; column < n; ++column)
        {
          result.startVec(column);
          for (typename Math::SparseMatrix<Scalar>::InnerIterator entry(m_core, column);
            entry; ++entry)
            result.insertBack(entry.row(), column) = entry.value();
          for (Eigen::Index k = 0; k < rank; ++k)
            result.insertBack(n + k, column) = Eigen::numext::conj(m_factor(column, k));
        }
        for (Eigen::Index k = 0; k < rank; ++k)
        {
          result.startVec(n + k);
          for (Eigen::Index row = 0; row < n; ++row)
            result.insertBack(row, n + k) = m_factor(row, k);
          result.insertBack(n + k, n + k) = Scalar(1);
        }
        result.finalize();
        return result;
      }

    private:
      /// @brief Borrowed sparse core.
      const Math::SparseMatrix<Scalar>& m_core;

      /// @brief Borrowed low-rank factor.
      const Math::Matrix<Scalar>& m_factor;
  };

  /**
   * @brief Linear solve of the local SWIFT metric or hinge tangent.
   *
   * ## Architecture
   *
   * The centered metric is represented as A - U U^T. Direct solvers use the
   * equivalent augmented system with a positive identity auxiliary block;
   * Matrix-free CG applies the low-rank subtraction without forming a dense matrix.
   * Native direct solvers retain factorizations where supported. Acceptance
   * delegates to the selected solver's success() contract; this adapter adds
   * no independent residual threshold. No gauge is applied.
   */
  template <class LinearSystem>
  class LinearSolver final : public Solver::LinearSolverBase<LinearSystem>
  {
    public:
      /// @brief Native linear-solver interface.
      using Parent = Solver::LinearSolverBase<LinearSystem>;

      /// @brief Linear-system type handled by this solver.
      using LinearSystemType = LinearSystem;
      using Parent::solve;

      /**
       * @brief Binds a solver to its variational problem.
       * @param problem Problem whose lifetime includes all solves.
       */
      explicit LinearSolver(typename Parent::ProblemBaseType& problem)
        : Parent(problem)
      {}

      /**
       * @brief Copies configuration and problem binding, not backend factorization resources.
       * @param other Solver whose configuration is copied.
       */
      LinearSolver(const LinearSolver& other)
        : Parent(other),
          m_centering(other.m_centering),
          m_parameters(other.m_parameters),
          m_report(other.m_report)
      {}

      /**
       * @brief Copies solver and convergence controls.
       * @param parameters Backend and residual-tolerance configuration.
       * @returns This solver.
       */
      LinearSolver& setParameters(const Parameters& parameters)
      {
        m_parameters = parameters;
        return *this;
      }

      /**
       * @brief Copies the centered metric's low-rank factor.
       * @param couplings Factor of the mean-strain subtraction.
       * @returns This solver.
       */
      LinearSolver& setCentering(const Math::Matrix<Real>& couplings)
      {
        m_centering = couplings;
        return *this;
      }

      /**
       * @brief Binds an optional destination for factorization counters.
       * @param report Diagnostics destination, or null to disable counters.
       * @returns This solver.
       */
      LinearSolver& setReport(Report* report)
      {
        m_report = report;
        return *this;
      }
      void solve(LinearSystemType& system) override
      {
        const auto start = std::chrono::steady_clock::now();
        m_iterations = 0;
        m_error = std::numeric_limits<Real>::infinity();
        m_success = solveSystem(system);
        m_seconds =
          std::chrono::duration<Real>(std::chrono::steady_clock::now() - start).count();
      }

      /**
       * @brief Reports the selected backend's solve status.
       * @returns Whether the selected solver reported success.
       */
      bool success() const noexcept
      {
        return m_success;
      }

      /**
       * @brief Returns the last iterative solve count.
       * @returns CG iterations, or zero for a direct solve.
       */
      size_t getIterations() const noexcept
      {
        return m_iterations;
      }

      /**
       * @brief Returns the iterative backend's relative residual estimate.
       * @returns CG's estimate, or NaN for a direct solve without that diagnostic.
       */
      Real getError() const noexcept
      {
        return m_error;
      }

      /**
       * @brief Returns the last solve duration.
       * @returns Elapsed wall-clock seconds including factorization.
       */
      Real getSolveSeconds() const noexcept
      {
        return m_seconds;
      }
      LinearSolver* copy() const noexcept override
      {
        return new LinearSolver(*this);
      }

    private:
      /// @brief Jacobi scaling from the full centered operator diagonal.
      class Preconditioner final : public Eigen::DiagonalPreconditioner<Real>
      {
        public:
          /**
           * @brief Computes inverse diagonal entries, using one for zero entries.
           * @param op Complete operator to precondition.
           * @returns This preconditioner.
           */
          Preconditioner& compute(const LowRankUpdate<Real>& op)
          {
            m_invdiag = op.diagonal();
            for (Eigen::Index i = 0; i < m_invdiag.size(); ++i)
            {
              const Real value = m_invdiag(i);
              m_invdiag(i) = value != Real(0) ? Real(1) / value : Real(1);
            }
            m_isInitialized = true;
            return *this;
          }
      };

      /**
       * @brief Adapts the centered operator and delegates its solution.
       * @param axb System whose solution is overwritten by the selected backend.
       * @returns The selected backend's success status.
       */
      bool solveSystem(LinearSystemType& axb)
      {
        auto& iterations = m_iterations;
        auto& error = m_error;
        auto* report = m_report;
        using OperatorType =
          typename FormLanguage::Traits<LinearSystemType>::OperatorType;
        using VectorType = typename FormLanguage::Traits<LinearSystemType>::VectorType;
        if constexpr (std::is_same_v<OperatorType, Math::SparseMatrix<Real>> &&
          std::is_same_v<VectorType, Math::Vector<Real>>)
        {
          const auto& rhs = axb.getVector();
          const LowRankUpdate<Real> op(axb.getOperator(), m_centering);
          if (m_parameters.linear.solver != Parameters::LinearSolver::CG)
          {
            const auto solveDirect = [&](auto& direct, StringView backend) {
              const auto solveSystem = [&](LinearSystemType& system,
                                         bool transferOperator) {
                if constexpr (requires { direct.factorize(system); })
                {
                  auto& next = system.getOperator();
                  next.makeCompressed();
                  const auto& previous = m_directSystem.getOperator();
                  const bool samePattern = next.rows() == previous.rows() &&
                    next.cols() == previous.cols() &&
                    next.nonZeros() == previous.nonZeros() &&
                    std::equal(next.outerIndexPtr(),
                      next.outerIndexPtr() + next.outerSize() + 1,
                      previous.outerIndexPtr()) &&
                    std::equal(next.innerIndexPtr(),
                      next.innerIndexPtr() + next.nonZeros(), previous.innerIndexPtr());
                  const bool sameValues = samePattern &&
                    std::equal(next.valuePtr(), next.valuePtr() + next.nonZeros(),
                      previous.valuePtr());
                  if (!samePattern)
                    direct.clear(Solver::Factorization::Symbolic);
                  const bool numeric = direct.success() &&
                    direct.getInfo().factorization == Solver::Factorization::Numeric;
                  // The augmented system is temporary; transfer its storage.
                  // A caller-owned core must remain unchanged after the solve.
                  if (transferOperator)
                    m_directSystem.getOperator() = std::move(next);
                  else
                    m_directSystem.getOperator() = next;
                  m_directSystem.getVector() = system.getVector();
                  if (!sameValues || !numeric)
                  {
                    if (report)
                    {
                      ++report->directFactorizations;
                      if (!direct.getInfo().factorization)
                        ++report->directAnalyses;
                    }
                    direct.factorize(m_directSystem);
                  }
                  if (direct.success())
                    direct.solve(m_directSystem);
                  if (direct.success())
                    system.getSolution() = m_directSystem.getSolution();
                }
                else
                  direct.solve(system);
              };
              const auto& matrix = axb.getOperator();
              if (m_centering.size() == 0)
                solveSystem(axb, false);
              else
              {
                const auto n = matrix.rows();
                LinearSystemType augmented;
                augmented.getOperator() = op.getAugmented();
                augmented.getVector() =
                  Math::Vector<Real>::Zero(augmented.getOperator().rows());
                augmented.getVector().head(n) = rhs;
                solveSystem(augmented, true);
                if (direct.success())
                  axb.getSolution() = augmented.getSolution().head(n);
              }
              iterations = 0;
              error = std::numeric_limits<Real>::quiet_NaN();
              const bool ok = direct.success();
              if (m_parameters.trace)
              {
                std::cout << "        metric " << backend << ": ok=" << ok
                          << "  status=" << direct.getInfo().status;
                if constexpr (requires { direct.getResources().instance.infog[1]; })
                  std::cout << "  detail=" << direct.getResources().instance.infog[1];
                std::cout << '\n';
              }
              return ok;
            };
#ifdef RODIN_USE_MUMPS
            if (m_parameters.linear.solver == Parameters::LinearSolver::MUMPS)
            {
              if (!m_mumps)
              {
                m_mumps =
                  std::make_unique<Solver::MUMPS<LinearSystemType>>(this->getProblem());
              }
              m_mumps->setSymmetric(Solver::MUMPS<LinearSystemType>::Symmetry::General)
                .setMaxThreads(m_parameters.linear.threads);
              return solveDirect(*m_mumps, "MUMPS");
            }
#endif
            Solver::SparseLU<LinearSystemType> direct(this->getProblem());
            return solveDirect(direct, "LU");
          }
          const std::size_t maxIterations = m_parameters.convergence.iterations.linear > 0
            ? m_parameters.convergence.iterations.linear
            : std::min<std::size_t>(AutomaticCGMaxIterations,
                std::max<std::size_t>(AutomaticCGMinIterations,
                  AutomaticCGIterationsPerDOF * axb.getOperator().rows()));
          Eigen::ConjugateGradient<LowRankUpdate<Real>, Eigen::Lower | Eigen::Upper,
            Preconditioner>
            cg;
          cg.setTolerance(m_parameters.convergence.tolerance.linearRelative)
            .setMaxIterations(maxIterations);
          cg.compute(op);
          if (axb.getSolution().size() == rhs.size())
            axb.getSolution() = cg.solveWithGuess(rhs, axb.getSolution());
          else
            axb.getSolution() = cg.solve(rhs);
          iterations = cg.iterations();
          error = cg.error();
          const bool ok = cg.info() == Eigen::Success;
          if (m_parameters.trace && !ok)
          {
            std::cout << "        cg failure: it=" << iterations
                      << "  max_it=" << maxIterations << "  err=" << error << '\n';
          }
          return ok;
        }
        else
        {
          Alert::MemberFunctionException(*this, __func__)
            << "SWIFT requires the local Eigen backend." << Alert::Raise;
          return false;
        }
      }

      /// Bounds on the automatic CG budget when no explicit cap is provided.
      static constexpr size_t AutomaticCGMinIterations = 100;
      static constexpr size_t AutomaticCGMaxIterations = 2000;

      /// Iterations per algebraic DOF in the automatic-budget heuristic.
      static constexpr size_t AutomaticCGIterationsPerDOF = 2;
      Math::Matrix<Real> m_centering;
      LinearSystemType m_directSystem;
#ifdef RODIN_USE_MUMPS
      std::unique_ptr<Solver::MUMPS<LinearSystemType>> m_mumps;
#endif
      Parameters m_parameters;
      Report* m_report = nullptr;
      bool m_success = false;
      size_t m_iterations = 0;
      Real m_error = std::numeric_limits<Real>::infinity();
      Real m_seconds = 0;
  };
}
namespace Rodin::FormLanguage
{
  /// @brief Native solver traits for the SWIFT linear adapter.
  template <class LinearSystem>
  struct Traits<Adaptation::SWIFT::LinearSolver<LinearSystem>>
  {
    /// @brief Linear system assembled by the variational problem.
    using LinearSystemType = LinearSystem;
  };
}
#endif
