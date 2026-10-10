/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_LINEARSOLVER_H
#define RODIN_ADAPTATION_SWIFT_LINEARSOLVER_H
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <limits>
#include <vector>
#include "Rodin/Solver/LinearSolver.h"
#include "Rodin/Solver/MUMPS.h"
#include "Rodin/Solver/SparseLU.h"
#include "Parameters.h"
#include "Report.h"

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Linear solve of the local SWIFT metric or hinge tangent.
   * @par Architecture
   * The centered metric is represented as A - U U^T. Direct solvers use the
   * equivalent augmented system with a positive identity auxiliary block;
   * CG applies the low-rank subtraction without forming a dense matrix.
   * Native direct solvers retain factorizations where supported. Residuals
   * are checked against the centered physical operator. No gauge is applied.
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
       * @brief Reports whether the last solve met the physical residual tolerance.
       * @returns Whether the last solution is finite and residual-validated.
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
       * @brief Returns the last centered-system residual.
       * @returns Relative Euclidean residual norm.
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
      /**
       * @brief Applies the centered physical operator.
       * @param matrix Sparse uncentered operator.
       * @param x Input vector.
       * @param y Output vector, overwritten by the centered operator action.
       */
      void applyMetric(const Math::SparseMatrix<Real>& matrix,
        const Math::Vector<Real>& x, Math::Vector<Real>& y) const
      {
        y = matrix * x;
        if (m_centering.size() != 0)
          y -= m_centering * (m_centering.transpose() * x);
      }

      /**
       * @brief Solves the centered system with diagonally preconditioned CG.
       * @param A Sparse uncentered operator.
       * @param b Right-hand side.
       * @param x Initial guess and returned solution.
       * @param maxIterations Maximum CG corrections.
       * @param relativeTolerance Residual threshold relative to the load norm.
       * @param iterations Returned correction count.
       * @param error Returned relative residual norm.
       * @returns Whether the recurrence residual satisfies the requested tolerance.
       */
      bool metricConjugateGradient(const Math::SparseMatrix<Real>& A,
        const Math::Vector<Real>& b, Math::Vector<Real>& x, std::size_t maxIterations,
        Real relativeTolerance, std::size_t& iterations, Real& error) const
      {
        Math::Vector<Real> jacobi(A.rows());
        for (Eigen::Index i = 0; i < A.rows(); ++i)
        {
          const Real d = A.coeff(i, i) -
            (m_centering.size() == 0 ? Real(0) : m_centering.row(i).squaredNorm());
          jacobi(i) = (std::abs(d) > Real(0)) ? Real(1) / d : Real(1);
        }
        auto applyPreconditioner = [&](const Math::Vector<Real>& residual) {
          return Math::Vector<Real>(jacobi.cwiseProduct(residual));
        };
        const Real rhsNorm = b.norm();
        iterations = 0;
        if (!(rhsNorm > Real(0)))
        {
          x.setZero();
          error = Real(0);
          return true;
        }
        Math::Vector<Real> Ax;
        applyMetric(A, x, Ax);
        Math::Vector<Real> r = b - Ax;
        Math::Vector<Real> z = applyPreconditioner(r);
        Math::Vector<Real> p = z;
        Real rz = r.dot(z);
        const Real threshold = relativeTolerance * rhsNorm;
        Math::Vector<Real> Ap;
        for (std::size_t it = 0; it < maxIterations; ++it)
        {
          if (r.norm() <= threshold)
            break;
          applyMetric(A, p, Ap);
          const Real pAp = p.dot(Ap);
          if (!(pAp > Real(0)))
            break;
          const Real alpha = rz / pAp;
          x += alpha * p;
          r -= alpha * Ap;
          z = applyPreconditioner(r);
          const Real rzNext = r.dot(z);
          p = z + (rzNext / rz) * p;
          rz = rzNext;
          ++iterations;
        }
        error = r.norm() / rhsNorm;
        return std::isfinite(error) && error <= relativeTolerance;
      }

      /**
       * @brief Solves and validates a centered physical system.
       * @param axb System whose solution is overwritten by the selected backend.
       * @returns Whether the backend succeeded with a finite, residual-validated solution.
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
          if (m_parameters.linear.solver != Parameters::LinearSolver::CG)
          {
            const auto solveDirect = [&](auto& direct, StringView backend) {
              const auto solveSystem = [&](LinearSystemType& system) {
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
                solveSystem(axb);
              else
              {
                // Eliminating the auxiliary variables gives the centered operator A - U U^T.
                const auto n = matrix.rows();
                const auto rank = m_centering.cols();
                std::vector<Math::SparseTriplet<Real>> entries;
                entries.reserve(matrix.nonZeros() + 2 * n * rank + rank);
                for (Eigen::Index column = 0; column < matrix.outerSize(); ++column)
                {
                  for (Math::SparseMatrix<Real>::InnerIterator entry(matrix, column);
                       entry; ++entry)
                    entries.emplace_back(entry.row(), entry.col(), entry.value());
                }
                for (Eigen::Index k = 0; k < rank; ++k)
                {
                  for (Eigen::Index i = 0; i < n; ++i)
                  {
                    const Real value = m_centering(i, k);
                    entries.emplace_back(i, n + k, value);
                    entries.emplace_back(n + k, i, value);
                  }
                  entries.emplace_back(n + k, n + k, Real(1));
                }
                LinearSystemType augmented;
                augmented.getOperator().resize(n + rank, n + rank);
                augmented.getOperator().setFromTriplets(entries.begin(), entries.end());
                augmented.getVector() = Math::Vector<Real>::Zero(n + rank);
                augmented.getVector().head(n) = rhs;
                solveSystem(augmented);
                if (direct.success())
                  axb.getSolution() = augmented.getSolution().head(n);
              }
              iterations = 0;
              Math::Vector<Real> image;
              if (direct.success())
                applyMetric(matrix, axb.getSolution(), image);
              error = direct.success() ? (image - rhs).norm() /
                  std::max(rhs.norm(), std::numeric_limits<Real>::min())
                                       : std::numeric_limits<Real>::infinity();
              const bool ok = direct.success() && axb.getSolution().allFinite() &&
                error <= m_parameters.convergence.tolerance.linearRelative;
              if (m_parameters.trace)
              {
                std::cout << "        metric " << backend << ": ok=" << ok
                          << "  rel=" << error << "  status=" << direct.getInfo().status;
                if constexpr (requires { direct.getResources().instance.infog[1]; })
                  std::cout << "  detail=" << direct.getResources().instance.infog[1];
                std::cout << '\n';
              }
              if (!ok)
                return false;
              return true;
            };
#ifdef RODIN_USE_MUMPS
            if (m_parameters.linear.solver == Parameters::LinearSolver::MUMPS)
            {
              if (!m_mumps)
                m_mumps =
                  std::make_unique<Solver::MUMPS<LinearSystemType>>(this->getProblem());
              m_mumps->setSymmetric(Solver::MUMPS<LinearSystemType>::Symmetry::General)
                .setMaxThreads(m_parameters.linear.threads);
              return solveDirect(*m_mumps, "MUMPS");
            }
#endif
            Solver::SparseLU<LinearSystemType> direct(this->getProblem());
            return solveDirect(direct, "LU");
          }
          const auto& guess = axb.getSolution();
          const std::size_t maxIterations = m_parameters.convergence.iterations.linear > 0
            ? m_parameters.convergence.iterations.linear
            : std::min<std::size_t>(AutomaticCGMaxIterations,
                std::max<std::size_t>(AutomaticCGMinIterations,
                  AutomaticCGIterationsPerDOF * axb.getOperator().rows()));
          Math::Vector<Real> solution =
            (guess.size() == rhs.size()) ? guess : Math::Vector<Real>::Zero(rhs.size());
          const bool solved =
            metricConjugateGradient(axb.getOperator(), rhs, solution, maxIterations,
              m_parameters.convergence.tolerance.linearRelative, iterations, error);
          Math::Vector<Real> image;
          applyMetric(axb.getOperator(), solution, image);
          error =
            (image - rhs).norm() / std::max(rhs.norm(), std::numeric_limits<Real>::min());
          const bool ok = solved && std::isfinite(error) &&
            error <= m_parameters.convergence.tolerance.linearRelative;
          axb.getSolution() = solution;
          if (m_parameters.trace && (!ok || !solution.allFinite()))
            std::cout << "        cg failure: ok=" << ok << "  it=" << iterations
                      << "  max_it=" << maxIterations << "  err=" << error
                      << "  finite=" << solution.allFinite() << '\n';
          return ok && solution.allFinite();
        }
        else
        {
          Alert::Exception() << "SWIFT requires the local Eigen backend." << Alert::Raise;
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
