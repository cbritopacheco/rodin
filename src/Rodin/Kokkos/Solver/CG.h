/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_CG_H
#define RODIN_KOKKOS_SOLVER_CG_H

/**
 * @file
 * @brief Conjugate gradient on the device for a Kokkos::SparseMatrix.
 */

#include <cmath>
#include <limits>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Solver/CG.h"

#include "Rodin/Kokkos/LinearSystem.h"

#include "Krylov.h"

namespace Rodin::Solver::Internal
{
  /// @brief Writes @f$ p = D^{-1} r @f$ and sums @f$ r \cdot p @f$.
  template <class Backend>
  struct KokkosCGPrecondition
  {
    /// @brief Inverse of the diagonal of the matrix.
    typename KokkosKrylov<Backend>::Vector inverse;

    /// @brief Residual.
    typename KokkosKrylov<Backend>::Vector r;

    /// @brief Preconditioned residual, written as the first search direction.
    typename KokkosKrylov<Backend>::Vector p;

    /**
     * @brief Preconditions one entry and accumulates its contribution to
     * @f$ r \cdot p @f$.
     * @param i Entry index.
     * @param sum Reduction accumulator.
     */
    KOKKOS_INLINE_FUNCTION void operator()(const size_t i, Real& sum) const
    {
      const Real v = inverse(i) * r(i);
      p(i) = v;
      sum += r(i) * v;
    }
  };

  /**
   * @brief Updates the iterate and the residual of a conjugate gradient step,
   * writes the preconditioned residual and sums @f$ |r|^2 @f$ and
   * @f$ r \cdot z @f$.
   */
  template <class Backend>
  struct KokkosCGUpdate
  {
    /// @brief Iterate, updated in place.
    typename KokkosKrylov<Backend>::Vector x;

    /// @brief Search direction.
    typename KokkosKrylov<Backend>::Vector p;

    /// @brief Residual, updated in place.
    typename KokkosKrylov<Backend>::Vector r;

    /// @brief Product @f$ q = A p @f$.
    typename KokkosKrylov<Backend>::Vector q;

    /// @brief Inverse of the diagonal of the matrix.
    typename KokkosKrylov<Backend>::Vector inverse;

    /// @brief Preconditioned residual, written by the step.
    typename KokkosKrylov<Backend>::Vector z;

    /// @brief Step length along the search direction.
    Real alpha;

    /**
     * @brief Updates one entry and accumulates its contributions to
     * @f$ |r|^2 @f$ and @f$ r \cdot z @f$.
     * @param k Entry index.
     * @param rr Accumulator of @f$ |r|^2 @f$.
     * @param dz Accumulator of @f$ r \cdot z @f$.
     */
    KOKKOS_INLINE_FUNCTION void operator()(const size_t k, Real& rr, Real& dz) const
    {
      x(k) += alpha * p(k);
      const Real rk = r(k) - alpha * q(k);
      r(k) = rk;
      const Real zk = inverse(k) * rk;
      z(k) = zk;
      rr += rk * rk;
      dz += rk * zk;
    }
  };

  /// @brief Writes the next search direction @f$ p = z + \beta p @f$.
  template <class Backend>
  struct KokkosCGDirection
  {
    /// @brief Search direction, updated in place.
    typename KokkosKrylov<Backend>::Vector p;

    /// @brief Preconditioned residual.
    typename KokkosKrylov<Backend>::Vector z;

    /// @brief Ratio of the successive values of @f$ r \cdot z @f$.
    Real beta;

    /**
     * @brief Updates one entry of the search direction.
     * @param k Entry index.
     */
    KOKKOS_INLINE_FUNCTION void operator()(const size_t k) const
    {
      p(k) = z(k) + beta * p(k);
    }
  };

  /**
   * @brief Jacobi preconditioned conjugate gradient, following the iteration
   * and the stopping test of Eigen's ConjugateGradient.
   * @param execution Execution space instance.
   * @param A Matrix.
   * @param inverse Inverse of the diagonal of @f$ A @f$.
   * @param b Right-hand side.
   * @param x Initial guess, then solution.
   * @param r Work vector.
   * @param z Work vector.
   * @param p Work vector.
   * @param q Work vector.
   * @param tolerance Relative tolerance on @f$ \|b - Ax\| / \|b\| @f$.
   * @param maxIterations Iteration limit.
   * @returns Convergence flag, iteration count and relative residual.
   */
  template <class Backend>
  KokkosKrylovResult kokkosCG(
    const typename KokkosKrylov<Backend>::ExecutionSpace& execution,
    const typename KokkosKrylov<Backend>::CrsMatrix& A,
    const typename KokkosKrylov<Backend>::Vector& inverse,
    const typename KokkosKrylov<Backend>::Vector& b,
    const typename KokkosKrylov<Backend>::Vector& x,
    const typename KokkosKrylov<Backend>::Vector& r,
    const typename KokkosKrylov<Backend>::Vector& z,
    const typename KokkosKrylov<Backend>::Vector& p,
    const typename KokkosKrylov<Backend>::Vector& q, Real tolerance, size_t maxIterations)
  {
    using Krylov = KokkosKrylov<Backend>;
    using ExecutionSpace = typename Krylov::ExecutionSpace;
    using Policy = ::Kokkos::RangePolicy<ExecutionSpace>;
    const size_t n = b.extent(0);
    KokkosKrylovResult res;
    typename Krylov::SpMV product;

    const Real rhsNorm2 = Krylov::dot(execution, b, b);
    if (rhsNorm2 == Real(0))
    {
      ::Kokkos::deep_copy(execution, x, Real(0));
      res.converged = true;
      return res;
    }
    const Real threshold =
      std::max(tolerance * tolerance * rhsNorm2, std::numeric_limits<Real>::min());

    // r = b - A x
    ::Kokkos::deep_copy(execution, r, b);
    ::KokkosSparse::spmv(execution, &product, "N", Real(-1), A, x, Real(1), r);
    Real residualNorm2 = Krylov::dot(execution, r, r);
    if (residualNorm2 < threshold)
    {
      res.converged = true;
      res.error = std::sqrt(residualNorm2 / rhsNorm2);
      return res;
    }

    // p = D^{-1} r
    Real absNew = 0;
    ::Kokkos::parallel_reduce("Rodin::Solver::CG::precondition", Policy(execution, 0, n),
      KokkosCGPrecondition<Backend>{inverse, r, p}, absNew);

    size_t i = 0;
    while (i < maxIterations)
    {
      ::KokkosSparse::spmv(execution, &product, "N", Real(1), A, p, Real(0), q);
      const Real alpha = absNew / Krylov::dot(execution, p, q);

      // x += alpha p, r -= alpha q, z = D^{-1} r, with |r|^2 and r.z.
      Real rz = 0;
      ::Kokkos::parallel_reduce("Rodin::Solver::CG::update", Policy(execution, 0, n),
        KokkosCGUpdate<Backend>{x, p, r, q, inverse, z, alpha},
        ::Kokkos::Sum<Real>(residualNorm2), ::Kokkos::Sum<Real>(rz));
      if (residualNorm2 < threshold)
        break;

      const Real beta = rz / absNew;
      absNew = rz;
      ::Kokkos::parallel_for("Rodin::Solver::CG::direction", Policy(execution, 0, n),
        KokkosCGDirection<Backend>{p, z, beta});
      ++i;
    }
    res.iterations = i;
    res.error = std::sqrt(residualNorm2 / rhsNorm2);
    res.converged = res.error <= tolerance;
    return res;
  }
}

namespace Rodin::Solver
{
  /**
   * @ingroup CGSpecializations
   * @brief Jacobi preconditioned conjugate gradient on the device, for a
   * system whose matrix is a Kokkos::SparseMatrix.
   *
   * The iteration, the defaults and the stopping test are those of the
   * sparse Eigen CG: the solve stops when @f$ \|b - Ax\| < \text{tol}
   * \|b\| @f$, the default tolerance is the machine epsilon and the default
   * iteration limit is twice the size. The matrix is read in place by
   * Kokkos-Kernels; the right-hand side and the solution are copied between
   * the host and the device at every solve, and the solution is used as the
   * initial guess when its size matches.
   *
   * @tparam Backend Kokkos backend of the matrix.
   */
  template <class Backend>
  class CG<
    Math::LinearSystem<::Rodin::Kokkos::SparseMatrix<Real, Backend>, Math::Vector<Real>>>
    final : public LinearSolverBase<Math::LinearSystem<
              ::Rodin::Kokkos::SparseMatrix<Real, Backend>, Math::Vector<Real>>>
  {
    using Krylov = Internal::KokkosKrylov<Backend>;

    public:
      /// @brief Scalar type.
      using ScalarType = Real;

      /// @brief Host vector type.
      using VectorType = Math::Vector<Real>;

      /// @brief Device matrix type.
      using OperatorType = ::Rodin::Kokkos::SparseMatrix<Real, Backend>;

      /// @brief Linear system type.
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      /// @brief Problem type.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      /// @brief Parent class type.
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /**
       * @brief Constructs the solver of a problem.
       * @param pb Problem to solve.
       */
      CG(ProblemBaseType& pb)
        : Parent(pb)
      {}

      /**
       * @brief Copies the settings, not the device vectors.
       * @param other Solver to copy.
       */
      CG(const CG& other)
        : Parent(other),
          m_tolerance(other.m_tolerance),
          m_maxIterations(other.m_maxIterations),
          m_result(other.m_result)
      {}

      /**
       * @brief Moves the solver.
       * @param other Solver to move.
       */
      CG(CG&& other) = default;

      /// @brief Destructor.
      ~CG() = default;

      /**
       * @brief Sets the relative tolerance.
       * @param tol Tolerance on @f$ \|b - Ax\| / \|b\| @f$.
       * @returns This solver.
       */
      CG& setTolerance(const Real& tol)
      {
        m_tolerance = tol;
        return *this;
      }

      /**
       * @brief Sets the iteration limit.
       * @param maxIt Maximum number of iterations.
       * @returns This solver.
       */
      CG& setMaxIterations(size_t maxIt)
      {
        m_maxIterations = maxIt;
        return *this;
      }

      /**
       * @brief Solves the system on the device.
       * @param[in,out] axb System, whose solution is the initial guess when
       * its size matches.
       */
      void solve(LinearSystemType& axb) override
      {
        const auto& matrix = axb.getOperator();
        const size_t n = matrix.getRowCount();
        if (matrix.getColumnCount() != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "CG requires a square matrix." << Alert::Raise;
        }
        if (static_cast<size_t>(axb.getVector().size()) != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The right-hand side size does not match the matrix." << Alert::Raise;
        }
        const typename Krylov::ExecutionSpace execution;
        Krylov::reserve(m_inverse, n, "Rodin::Solver::CG::inverse");
        Krylov::reserve(m_b, n, "Rodin::Solver::CG::b");
        Krylov::reserve(m_x, n, "Rodin::Solver::CG::x");
        Krylov::reserve(m_r, n, "Rodin::Solver::CG::r");
        Krylov::reserve(m_z, n, "Rodin::Solver::CG::z");
        Krylov::reserve(m_p, n, "Rodin::Solver::CG::p");
        Krylov::reserve(m_q, n, "Rodin::Solver::CG::q");
        Krylov::upload(axb.getVector(), m_b);
        if (static_cast<size_t>(axb.getSolution().size()) == n)
          Krylov::upload(axb.getSolution(), m_x);
        else
          ::Kokkos::deep_copy(m_x, Real(0));
        if (n == 0)
        {
          m_result = {true, 0, 0};
          axb.getSolution().resize(0);
          return;
        }
        Krylov::invertDiagonal(execution, matrix, m_inverse);
        const auto A = Krylov::wrap(matrix);
        m_result = Internal::kokkosCG<Backend>(execution, A, m_inverse, m_b, m_x, m_r,
          m_z, m_p, m_q, m_tolerance, m_maxIterations ? *m_maxIterations : 2 * n);
        Krylov::download(m_x, axb.getSolution());
      }

      /**
       * @brief Whether the last solve reached the tolerance.
       * @returns True if the last solve converged.
       */
      Boolean success() const
      {
        return m_result.converged;
      }

      /**
       * @brief Number of iterations of the last solve.
       * @returns Iteration count of the last solve.
       */
      size_t getIterationNumber() const
      {
        return m_result.iterations;
      }

      /**
       * @brief Relative residual of the last solve.
       * @returns Relative residual estimated by the last solve.
       */
      Real getError() const
      {
        return m_result.error;
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new CG instance.
       */
      CG* copy() const noexcept override
      {
        return new CG(*this);
      }

    private:
      Real m_tolerance = std::numeric_limits<Real>::epsilon();
      Optional<size_t> m_maxIterations;
      Internal::KokkosKrylovResult m_result;
      typename Krylov::Vector m_inverse, m_b, m_x, m_r, m_z, m_p, m_q;
  };
}

#endif
