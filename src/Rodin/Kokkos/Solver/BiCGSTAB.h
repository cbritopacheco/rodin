/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_BICGSTAB_H
#define RODIN_KOKKOS_SOLVER_BICGSTAB_H

/**
 * @file
 * @brief Jacobi preconditioned BiCGSTAB for backend sparse systems.
 */

#include <cmath>
#include <limits>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Solver/BiCGSTAB.h"
#include "Rodin/Kokkos/Math/LinearSystem.h"
#include "Rodin/Kokkos/Utility.h"

#include "Krylov.h"

namespace Rodin::Solver
{
  /**
   * @ingroup BiCGSTABSpecializations
   * @brief Right Jacobi preconditioned BiCGSTAB in backend memory.
   *
   * ## Architecture
   * The operator is read through its native CSR view. Backend right-hand sides
   * and solutions are borrowed; explicitly selected host vectors use retained
   * transfer buffers. Vector updates and all work vectors remain on the backend.
   * Host iteration control receives scalar reductions. Copies keep settings and the last
   * diagnostic but start with independent work storage.
   *
   * The matching-size solution is the initial guess. Defaults match the native
   * solver: machine epsilon tolerance and twice the number of rows for the
   * iteration limit. A zero right-hand side returns the zero solution.
   * Convergence certifies @f$ \|b-Ax\|/\|b\| \leq \mathrm{tol} @f$ using a
   * freshly computed residual. If a recursive residual meets the tolerance but
   * the true residual does not, the recurrence restarts from the true residual
   * within the same iteration budget. Zero or non-finite recurrence divisors
   * report failure through @ref success rather than raising an exception.
   *
   * @tparam Backend Kokkos backend of the matrix.
   * @tparam VectorType_ Right-hand side and solution vector container.
   */
  template <class Backend, class VectorType_>
  class BiCGSTAB<::Rodin::Math::LinearSystem<
    ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>, VectorType_>>
    final : public LinearSolverBase<::Rodin::Math::LinearSystem<
              ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>, VectorType_>>
  {
    using Krylov = Internal::KokkosKrylov<Backend>;
    using WorkVector = typename Krylov::Vector;
    using ExecutionSpace = typename Krylov::ExecutionSpace;
    using Policy = ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<size_t>>;

    public:
      /// @brief Scalar type.
      using ScalarType = Real;

      /// @brief System coefficient vector type.
      using VectorType = VectorType_;

      /// @brief Backend sparse matrix type.
      using OperatorType = ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>;

      /// @brief Linear system type.
      using LinearSystemType = ::Rodin::Math::LinearSystem<OperatorType, VectorType>;

      /// @brief Problem type.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      /// @brief Parent class type.
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /// @param problem Problem whose system is solved.
      BiCGSTAB(ProblemBaseType& problem)
        : Parent(problem)
      {}

      /// @param other Solver whose settings and diagnostic are copied.
      BiCGSTAB(const BiCGSTAB& other)
        : Parent(other),
          m_tolerance(other.m_tolerance),
          m_maxIterations(other.m_maxIterations),
          m_result(other.m_result)
      {}

      /// @param other Solver whose settings and work storage are moved.
      BiCGSTAB(BiCGSTAB&& other) = default;

      /// @brief Destructor.
      ~BiCGSTAB() override = default;

      /**
       * @brief Sets the relative residual tolerance.
       * @param tolerance Finite, nonnegative tolerance.
       * @return This solver.
       */
      BiCGSTAB& setTolerance(const Real& tolerance)
      {
        if (!std::isfinite(tolerance) || tolerance < 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The tolerance must be finite and nonnegative." << Alert::Raise;
        }
        m_tolerance = tolerance;
        return *this;
      }

      /**
       * @brief Sets the maximum number of completed BiCGSTAB steps.
       * @param count Iteration limit; zero disables the recurrence.
       * @return This solver.
       */
      BiCGSTAB& setMaxIterations(size_t count)
      {
        m_maxIterations = count;
        return *this;
      }

      /**
       * @brief Solves a backend system with the current solution as initial guess.
       * @param system System to solve. A mismatched-size guess is replaced by zero.
       */
      void solve(LinearSystemType& system) override
      {
        const auto& matrix = system.getOperator();
        const size_t n = matrix.getRowCount();
        if (matrix.getColumnCount() != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "BiCGSTAB requires a square matrix." << Alert::Raise;
        }
        if (static_cast<size_t>(system.getVector().size()) != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The right-hand side size does not match the matrix." << Alert::Raise;
        }
        m_result = {false, 0, std::numeric_limits<Real>::infinity()};
        if (n == 0)
        {
          system.getSolution().resize(0);
          m_result = {true, 0, 0};
          return;
        }
        if constexpr (requires { system.getSolution().getView(); })
        {
          if (static_cast<size_t>(system.getSolution().size()) != n)
            system.getSolution().setZero(n);
        }
        const ExecutionSpace execution;
        if constexpr (!requires { system.getVector().getView(); })
          Krylov::reserve(m_b, n, "Rodin BiCGSTAB RHS");
        if constexpr (!requires { system.getSolution().getView(); })
          Krylov::reserve(m_x, n, "Rodin BiCGSTAB solution");
        Krylov::upload(system.getVector(), m_b);
        if (static_cast<size_t>(system.getSolution().size()) == n)
          Krylov::upload(system.getSolution(), m_x);
        else
          ::Kokkos::deep_copy(m_x, Real(0));
        Krylov::reserve(m_inverse, n, "Rodin BiCGSTAB inverse diagonal");
        Krylov::reserve(m_r, n, "Rodin BiCGSTAB residual");
        Krylov::reserve(m_shadow, n, "Rodin BiCGSTAB shadow residual");
        Krylov::reserve(m_p, n, "Rodin BiCGSTAB direction");
        Krylov::reserve(m_v, n, "Rodin BiCGSTAB direction product");
        Krylov::reserve(m_y, n, "Rodin BiCGSTAB preconditioned direction");
        Krylov::reserve(m_s, n, "Rodin BiCGSTAB intermediate residual");
        Krylov::reserve(m_z, n, "Rodin BiCGSTAB preconditioned residual");
        Krylov::reserve(m_t, n, "Rodin BiCGSTAB residual product");
        Krylov::invertDiagonal(execution, matrix, m_inverse);
        iterate(
          execution, Krylov::wrap(matrix), m_maxIterations ? *m_maxIterations : 2 * n);
        ::Rodin::Kokkos::Internal::checkLaunches("Rodin BiCGSTAB solve");
        Krylov::download(m_x, system.getSolution());
      }

      /// @return Whether the last solve certified the requested relative residual.
      Boolean success() const
      {
        return m_result.converged;
      }

      /// @return Number of completed steps of the last solve, including restarts.
      size_t getIterationNumber() const
      {
        return m_result.iterations;
      }

      /// @return True relative residual of the last solve, or infinity if non-finite.
      Real getError() const
      {
        return m_result.error;
      }

      /// @return New solver with copied settings and independent work storage.
      BiCGSTAB* copy() const noexcept override
      {
        return new BiCGSTAB(*this);
      }

      /// @brief Forms the next direction from the previous step.
      struct Direction
      {
        WorkVector r, p, v;
        Real beta, omega;

        /// @param i Coefficient index.
        KOKKOS_INLINE_FUNCTION void operator()(size_t i) const
        {
          p(i) = r(i) + beta * (p(i) - omega * v(i));
        }
      };

      /// @brief Forms the intermediate residual and reduces its squared norm.
      struct Intermediate
      {
        WorkVector r, v, s;
        Real alpha;

        /// @param i Coefficient index. @param norm Squared norm accumulator.
        KOKKOS_INLINE_FUNCTION void operator()(size_t i, Real& norm) const
        {
          const Real value = r(i) - alpha * v(i);
          s(i) = value;
          norm += value * value;
        }
      };

      /// @brief Updates the iterate and residual and reduces the squared norm.
      struct Update
      {
        WorkVector x, y, z, s, t, r;
        Real alpha, omega;

        /// @param i Coefficient index. @param norm Squared norm accumulator.
        KOKKOS_INLINE_FUNCTION void operator()(size_t i, Real& norm) const
        {
          x(i) += alpha * y(i) + omega * z(i);
          const Real value = s(i) - omega * t(i);
          r(i) = value;
          norm += value * value;
        }
      };

      /// @brief Commits the alpha step when the intermediate residual converges.
      struct AlphaUpdate
      {
        WorkVector x, y;
        Real alpha;

        /// @param i Coefficient index.
        KOKKOS_INLINE_FUNCTION void operator()(size_t i) const
        {
          x(i) += alpha * y(i);
        }
      };

    private:
      /**
       * @brief Recomputes the true residual and updates the convergence diagnostic.
       * @param execution Execution instance.
       * @param matrix Native CSR matrix.
       * @param product SpMV descriptors retained for this solve.
       * @param rhsNorm2 Squared norm of the nonzero right-hand side.
       */
      void residual(const ExecutionSpace& execution,
        const typename Krylov::CrsMatrix& matrix, typename Krylov::SpMV& product,
        Real rhsNorm2)
      {
        ::Kokkos::deep_copy(execution, m_r, m_b);
        ::KokkosSparse::spmv(
          execution, &product, "N", Real(-1), matrix, m_x, Real(1), m_r);
        ::Rodin::Kokkos::Internal::checkLaunches("Rodin BiCGSTAB residual product");
        const Real norm2 = Krylov::dot(execution, m_r, m_r);
        m_result.error = std::isfinite(norm2) ? std::sqrt(norm2 / rhsNorm2)
                                              : std::numeric_limits<Real>::infinity();
        m_result.converged =
          std::isfinite(m_result.error) && m_result.error <= m_tolerance;
      }

      /**
       * @brief Executes the right-preconditioned BiCGSTAB recurrence.
       * @param execution Execution instance.
       * @param matrix Native CSR matrix.
       * @param limit Maximum number of completed steps.
       */
      void iterate(const ExecutionSpace& execution,
        const typename Krylov::CrsMatrix& matrix, size_t limit)
      {
        typename Krylov::SpMV product;
        const Real rhsNorm2 = Krylov::dot(execution, m_b, m_b);
        if (rhsNorm2 == 0)
        {
          ::Kokkos::deep_copy(execution, m_x, Real(0));
          m_result = {true, 0, 0};
          execution.fence("Rodin BiCGSTAB zero RHS");
          return;
        }
        if (!std::isfinite(rhsNorm2))
          return;
        residual(execution, matrix, product, rhsNorm2);
        if (m_result.converged || !std::isfinite(m_result.error))
          return;
        ::Kokkos::deep_copy(execution, m_shadow, m_r);
        const size_t n = m_r.extent(0);
        Real rhoOld = 1, alpha = 1, omega = 1;
        bool first = true;
        while (m_result.iterations < limit)
        {
          const Real rho = Krylov::dot(execution, m_shadow, m_r);
          if (!std::isfinite(rho) || rho == 0)
            break;
          if (first)
            ::Kokkos::deep_copy(execution, m_p, m_r);
          else
          {
            const Real beta = (rho / rhoOld) * (alpha / omega);
            if (!std::isfinite(beta))
              break;
            ::Kokkos::parallel_for("Rodin BiCGSTAB direction", Policy(execution, 0, n),
              Direction{m_r, m_p, m_v, beta, omega});
          }
          Krylov::scale(execution, m_inverse, m_p, m_y, Real(1), Real(0));
          ::KokkosSparse::spmv(
            execution, &product, "N", Real(1), matrix, m_y, Real(0), m_v);
          ::Rodin::Kokkos::Internal::checkLaunches("Rodin BiCGSTAB direction product");
          const Real divisor = Krylov::dot(execution, m_shadow, m_v);
          if (!std::isfinite(divisor) || divisor == 0)
            break;
          alpha = rho / divisor;
          if (!std::isfinite(alpha))
            break;
          Real norm2 = 0;
          ::Kokkos::parallel_reduce("Rodin BiCGSTAB intermediate",
            Policy(execution, 0, n), Intermediate{m_r, m_v, m_s, alpha}, norm2);
          ::Rodin::Kokkos::Internal::checkLaunches("Rodin BiCGSTAB intermediate");
          if (!std::isfinite(norm2))
            break;
          if (std::sqrt(norm2 / rhsNorm2) <= m_tolerance)
          {
            ::Kokkos::parallel_for("Rodin BiCGSTAB alpha update", Policy(execution, 0, n),
              AlphaUpdate{m_x, m_y, alpha});
            ::Rodin::Kokkos::Internal::checkLaunches("Rodin BiCGSTAB alpha update");
            ++m_result.iterations;
          }
          else
          {
            Krylov::scale(execution, m_inverse, m_s, m_z, Real(1), Real(0));
            ::KokkosSparse::spmv(
              execution, &product, "N", Real(1), matrix, m_z, Real(0), m_t);
            ::Rodin::Kokkos::Internal::checkLaunches("Rodin BiCGSTAB residual product");
            const Real tt = Krylov::dot(execution, m_t, m_t);
            const Real ts = Krylov::dot(execution, m_t, m_s);
            if (!std::isfinite(tt) || !std::isfinite(ts) || tt == 0)
              break;
            omega = ts / tt;
            if (!std::isfinite(omega))
              break;
            ::Kokkos::parallel_reduce("Rodin BiCGSTAB update", Policy(execution, 0, n),
              Update{m_x, m_y, m_z, m_s, m_t, m_r, alpha, omega}, norm2);
            ::Rodin::Kokkos::Internal::checkLaunches("Rodin BiCGSTAB update");
            ++m_result.iterations;
            if (std::isfinite(norm2) && std::sqrt(norm2 / rhsNorm2) > m_tolerance &&
              omega != 0)
            {
              rhoOld = rho;
              first = false;
              continue;
            }
          }
          residual(execution, matrix, product, rhsNorm2);
          if (m_result.converged || !std::isfinite(m_result.error) || omega == 0)
            return;
          // A fresh shadow and direction restart from the verified residual.
          ::Kokkos::deep_copy(execution, m_shadow, m_r);
          rhoOld = alpha = omega = 1;
          first = true;
        }
        residual(execution, matrix, product, rhsNorm2);
      }

      Real m_tolerance = std::numeric_limits<Real>::epsilon();
      Optional<size_t> m_maxIterations;
      Internal::KokkosKrylovResult m_result;
      WorkVector m_inverse, m_b, m_x, m_r, m_shadow, m_p, m_v, m_y, m_s, m_z, m_t;
  };
}

#endif
