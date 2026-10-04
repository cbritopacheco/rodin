/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIR_LINEARSOLVER_H
#define RODIN_ADAPTATION_WNGIR_LINEARSOLVER_H
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <limits>
#include <vector>
#include <Eigen/Eigenvalues>
#include "Rodin/Solver/LinearSolver.h"
#include "Rodin/Solver/MUMPS.h"
#include "Rodin/Solver/SparseLU.h"
#include "Parameters.h"
#include "Report.h"

namespace Rodin::Adaptation
{
  /**
   * @brief Linear solve of the local WNGIR metric or hinge tangent.
   * @par Architecture
   * Compatible unresolved similarity modes are selected in a small restricted
   * eigensolve and gauged only for the solve. Native direct solvers retain
   * factorizations where supported. CG applies the sparse matrix and gauges
   * without a dense low-rank matrix. Acceptance checks the original residual.
   */
  template <class LinearSystem>
  class WNGIRLinearSolver final : public Solver::LinearSolverBase<LinearSystem>
  {
    public:
      using Parent = Solver::LinearSolverBase<LinearSystem>;
      using LinearSystemType = LinearSystem;
      using Parent::solve;
      explicit WNGIRLinearSolver(typename Parent::ProblemBaseType& problem)
        : Parent(problem)
      {}
      /// Copies configuration and problem binding, not backend factorization resources.
      WNGIRLinearSolver(const WNGIRLinearSolver& other)
        : Parent(other),
          m_similarityModes(other.m_similarityModes),
          m_parameters(other.m_parameters),
          m_report(other.m_report)
      {}
      WNGIRLinearSolver& setParameters(const WNGIRParameters& parameters)
      {
        m_parameters = parameters;
        return *this;
      }
      WNGIRLinearSolver& setSimilarityModes(const Math::Matrix<Real>& modes)
      {
        m_similarityModes = modes;
        return *this;
      }
      WNGIRLinearSolver& setReport(WNGIRReport* report)
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
      bool success() const noexcept
      {
        return m_success;
      }
      size_t getIterations() const noexcept
      {
        return m_iterations;
      }
      Real getError() const noexcept
      {
        return m_error;
      }
      Real getSolveSeconds() const noexcept
      {
        return m_seconds;
      }
      WNGIRLinearSolver* copy() const noexcept override
      {
        return new WNGIRLinearSolver(*this);
      }

    private:
      struct SimilarityGauge
      {
          std::vector<Math::Vector<Real>> modes;
          std::vector<Real> weights;
      };
      size_t gaugeSimilarityModes(const Math::SparseMatrix<Real>& matrix,
        const Math::Vector<Real>& force, SimilarityGauge& projection) const
      {
        const auto& modes = m_similarityModes;
        const size_t count = modes.cols();
        if (count == 0)
          return 0;
        Math::Matrix<Real> images(matrix.rows(), count);
        for (size_t column = 0; column < count; ++column)
        {
          Math::Vector<Real> image;
          applyMetric(matrix, projection, modes.col(column), image);
          images.col(column) = image;
        }
        Math::Matrix<Real> gram = modes.transpose() * images;
        gram = (Real(0.5) * (gram + gram.transpose())).eval();
        Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eigen(gram);
        if (eigen.info() != Eigen::Success)
          return 0;
        const Real scale = matrix.diagonal().cwiseAbs().maxCoeff();
        const Real tolerance =
          NullModeRoundoffFactor * std::numeric_limits<Real>::epsilon();
        size_t unresolved = 0;
        for (size_t column = 0; column < count; ++column)
        {
          const Math::Vector<Real> candidate = modes * eigen.eigenvectors().col(column);
          if ((images * eigen.eigenvectors().col(column)).norm() <= tolerance * scale &&
            std::abs(force.dot(candidate)) <= tolerance * force.norm())
          {
            projection.modes.push_back(candidate);
            projection.weights.push_back(scale);
            ++unresolved;
          }
        }
        return unresolved;
      }

      void applyMetric(const Math::SparseMatrix<Real>& A,
        const SimilarityGauge& projection, const Math::Vector<Real>& x,
        Math::Vector<Real>& y) const
      {
        y = A * x;
        for (std::size_t k = 0; k < projection.weights.size(); ++k)
          y += projection.weights[k] * projection.modes[k].dot(x) * projection.modes[k];
      }

      bool metricConjugateGradient(const Math::SparseMatrix<Real>& A,
        const SimilarityGauge& projection, const Math::Vector<Real>& b,
        Math::Vector<Real>& x, std::size_t maxIterations, Real relativeTolerance,
        std::size_t& iterations, Real& error) const
      {
        Math::Vector<Real> jacobi(A.rows());
        for (Eigen::Index i = 0; i < A.rows(); ++i)
        {
          const Real d = A.coeff(i, i);
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
        applyMetric(A, projection, x, Ax);
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
          applyMetric(A, projection, p, Ap);
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
          SimilarityGauge gauged;
          const size_t unresolved = gaugeSimilarityModes(axb.getOperator(), rhs, gauged);
          if (report)
            report->unresolvedSimilarityModes = unresolved;
          if (m_parameters.directSolver != WNGIRParameters::DirectSolver::CG)
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
              const auto& rigid = gauged;
              if (rigid.weights.empty())
                solveSystem(axb);
              else
              {
                // Positive solve-only gauges use an auxiliary negative identity block.
                const auto n = matrix.rows();
                const auto rank = static_cast<Eigen::Index>(rigid.weights.size());
                std::vector<Math::SparseTriplet<Real>> entries;
                entries.reserve(matrix.nonZeros() + 2 * n * rank + rank);
                for (Eigen::Index column = 0; column < matrix.outerSize(); ++column)
                  for (Math::SparseMatrix<Real>::InnerIterator entry(matrix, column);
                       entry; ++entry)
                    entries.emplace_back(entry.row(), entry.col(), entry.value());
                for (Eigen::Index k = 0; k < rank; ++k)
                {
                  const Real scale = std::sqrt(rigid.weights[k]);
                  for (Eigen::Index i = 0; i < n; ++i)
                  {
                    const Real value = scale * rigid.modes[k](i);
                    entries.emplace_back(i, n + k, value);
                    entries.emplace_back(n + k, i, value);
                  }
                  entries.emplace_back(n + k, n + k, Real(-1));
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
                image = matrix * axb.getSolution();
              error = direct.success() ? (image - rhs).norm() /
                  std::max(rhs.norm(), std::numeric_limits<Real>::min())
                                       : std::numeric_limits<Real>::infinity();
              const bool ok = direct.success() && axb.getSolution().allFinite() &&
                error <= m_parameters.cgRelativeTolerance;
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
            if (m_parameters.directSolver == WNGIRParameters::DirectSolver::MUMPS)
            {
              if (!m_mumps)
                m_mumps =
                  std::make_unique<Solver::MUMPS<LinearSystemType>>(this->getProblem());
              m_mumps->setSymmetric(Solver::MUMPS<LinearSystemType>::Symmetry::General)
                .setMaxThreads(m_parameters.directSolverThreads);
              return solveDirect(*m_mumps, "MUMPS");
            }
#endif
            Solver::SparseLU<LinearSystemType> direct(this->getProblem());
            return solveDirect(direct, "LU");
          }
          const auto& guess = axb.getSolution();
          const std::size_t maxIterations = m_parameters.cgMaxIterations > 0
            ? m_parameters.cgMaxIterations
            : std::min<std::size_t>(AutomaticCGMaxIterations,
                std::max<std::size_t>(AutomaticCGMinIterations,
                  AutomaticCGIterationsPerDOF * axb.getOperator().rows()));
          Math::Vector<Real> solution =
            (guess.size() == rhs.size()) ? guess : Math::Vector<Real>::Zero(rhs.size());
          const bool solved = metricConjugateGradient(axb.getOperator(), gauged, rhs,
            solution, maxIterations, m_parameters.cgRelativeTolerance, iterations, error);
          Math::Vector<Real> image;
          image = axb.getOperator() * solution;
          error =
            (image - rhs).norm() / std::max(rhs.norm(), std::numeric_limits<Real>::min());
          const bool ok =
            solved && std::isfinite(error) && error <= m_parameters.cgRelativeTolerance;
          axb.getSolution() = solution;
          if (m_parameters.trace && (!ok || !solution.allFinite()))
            std::cout << "        cg failure: ok=" << ok << "  it=" << iterations
                      << "  max_it=" << maxIterations << "  err=" << error
                      << "  finite=" << solution.allFinite() << '\n';
          return ok && solution.allFinite();
        }
        else
        {
          Alert::Exception() << "WNGIR requires the local Eigen backend." << Alert::Raise;
          return false;
        }
      }

      /// Relative roundoff threshold for compatible null modes; heuristic.
      static constexpr Real NullModeRoundoffFactor = Real(256);
      /// Bounds on the automatic CG budget when no explicit cap is provided.
      static constexpr size_t AutomaticCGMinIterations = 100;
      static constexpr size_t AutomaticCGMaxIterations = 2000;
      /// Iterations per algebraic DOF in the automatic-budget heuristic.
      static constexpr size_t AutomaticCGIterationsPerDOF = 2;
      Math::Matrix<Real> m_similarityModes;
      LinearSystemType m_directSystem;
#ifdef RODIN_USE_MUMPS
      std::unique_ptr<Solver::MUMPS<LinearSystemType>> m_mumps;
#endif
      WNGIRParameters m_parameters;
      WNGIRReport* m_report = nullptr;
      bool m_success = false;
      size_t m_iterations = 0;
      Real m_error = std::numeric_limits<Real>::infinity();
      Real m_seconds = 0;
  };
}
namespace Rodin::FormLanguage
{
  template <class LinearSystem>
  struct Traits<Adaptation::WNGIRLinearSolver<LinearSystem>>
  {
      using LinearSystemType = LinearSystem;
  };
}
#endif
