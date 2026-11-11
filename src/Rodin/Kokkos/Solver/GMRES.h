/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_GMRES_H
#define RODIN_KOKKOS_SOLVER_GMRES_H

/**
 * @file
 * @brief Restarted GMRES on the device for a Kokkos::Math::SparseMatrix.
 */

#include <limits>

#include <KokkosKernels_Handle.hpp>
#include <KokkosSparse_gmres.hpp>
#include <KokkosSparse_Preconditioner.hpp>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Solver/GMRES.h"

#include "Rodin/Kokkos/Math/LinearSystem.h"

#include "Krylov.h"

namespace Rodin::Solver::Internal
{
  /**
   * @brief Jacobi preconditioner of the Kokkos-Kernels GMRES.
   * @tparam Backend Kokkos backend of the matrix.
   */
  template <class Backend>
  class KokkosJacobi final : public ::KokkosSparse::Experimental::Preconditioner<
                               typename KokkosKrylov<Backend>::CrsMatrix>
  {
    public:
      /// @brief Krylov storage type.
      using Krylov = KokkosKrylov<Backend>;

      /// @brief Kokkos device of the vectors.
      using Device = typename Krylov::Device;

      /**
       * @brief Preconditioner holding the inverse of a diagonal.
       * @param inverse Inverse of the diagonal.
       */
      explicit KokkosJacobi(const typename Krylov::Vector& inverse)
        : m_inverse(inverse)
      {}

      /**
       * @brief Computes @f$ Y = \beta Y + \alpha D^{-1} X @f$.
       * @param X Input vector.
       * @param Y Output vector.
       * @param trans Transposition mode, ignored since the diagonal equals
       * its transpose.
       * @param alpha Coefficient of @f$ D^{-1} X @f$.
       * @param beta Coefficient of @f$ Y @f$.
       */
      void apply(const ::Kokkos::View<const Real*, Device>& X,
        const ::Kokkos::View<Real*, Device>& Y, [[maybe_unused]] const char trans[] = "N",
        Real alpha = Real(1), Real beta = Real(0)) const override
      {
        Krylov::scale(typename Krylov::ExecutionSpace(), m_inverse, X, Y, alpha, beta);
      }

      /// @brief Has no parameter.
      void setParameters() override {}

      /// @brief Has no structure to set up.
      void initialize() override {}

      /**
       * @brief Always initialized.
       * @returns True.
       */
      bool isInitialized() const override
      {
        return true;
      }

      /// @brief The diagonal is computed by the solver.
      void compute() override {}

      /**
       * @brief Always computed.
       * @returns True.
       */
      bool isComputed() const override
      {
        return true;
      }

      /**
       * @brief A diagonal is its own transpose.
       * @returns True.
       */
      bool hasTransposeApply() const override
      {
        return true;
      }

    private:
      typename Krylov::Vector m_inverse;
  };
}

namespace Rodin::Solver
{
  /**
   * @ingroup GMRESSpecializations
   * @brief Restarted GMRES of Kokkos-Kernels on the device, right
   * preconditioned by Jacobi, for a system whose matrix is a
   * Kokkos::Math::SparseMatrix.
   *
   * The defaults are those of the sparse Eigen GMRES: the machine epsilon as
   * tolerance on @f$ \|b - Ax\| / \|b\| @f$, twice the size as iteration
   * limit and a restart of 30. The iteration limit is rounded up to a whole
   * number of restart cycles. The matrix is read in place; the right-hand
   * side and solution are borrowed directly for backend vectors, while
   * explicit host vectors are staged. The solution supplies the initial guess
   * when its size matches.
   *
   * @tparam Backend Kokkos backend of the matrix.
   */
  template <class Backend, class VectorType_>
  class GMRES<::Rodin::Math::LinearSystem<
    ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>, VectorType_>>
    final : public LinearSolverBase<::Rodin::Math::LinearSystem<
              ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>, VectorType_>>
  {
    using Krylov = Internal::KokkosKrylov<Backend>;

    using Handle =
      ::KokkosKernels::Experimental::KokkosKernelsHandle<::Rodin::Math::SparseIndex,
        ::Rodin::Math::SparseIndex, Real, typename Krylov::ExecutionSpace,
      typename Krylov::MemorySpace, typename Krylov::MemorySpace>;

    public:
      /// @brief Scalar type.
      using ScalarType = Real;

      /// @brief System coefficient vector type.
      using VectorType = VectorType_;

      /// @brief Device matrix type.
      using OperatorType = ::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>;

      /// @brief Linear system type.
      using LinearSystemType = ::Rodin::Math::LinearSystem<OperatorType, VectorType>;

      /// @brief Problem type.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      /// @brief Parent class type.
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /**
       * @brief Constructs the solver of a problem.
       * @param problem Problem to solve.
       */
      explicit GMRES(ProblemBaseType& problem)
        : Parent(problem)
      {}

      /**
       * @brief Copies the settings, not the device vectors.
       * @param other Solver to copy.
       */
      GMRES(const GMRES& other)
        : Parent(other),
          m_tolerance(other.m_tolerance),
          m_maxIterations(other.m_maxIterations),
          m_restart(other.m_restart),
          m_result(other.m_result)
      {}

      /**
       * @brief Moves the solver.
       * @param other Solver to move.
       */
      GMRES(GMRES&& other) = default;

      /// @brief Destructor.
      ~GMRES() = default;

      /**
       * @brief Sets the relative tolerance.
       * @param tolerance Tolerance on @f$ \|b - Ax\| / \|b\| @f$.
       * @returns This solver.
       */
      GMRES& setTolerance(const Real& tolerance)
      {
        m_tolerance = tolerance;
        return *this;
      }

      /**
       * @brief Sets the iteration limit, rounded up to whole restart cycles.
       * @param count Maximum number of iterations.
       * @returns This solver.
       */
      GMRES& setMaxIterations(size_t count)
      {
        m_maxIterations = count;
        return *this;
      }

      /**
       * @brief Sets the restart dimension.
       * @param count Krylov dimension between restarts.
       * @returns This solver.
       */
      GMRES& setRestart(size_t count)
      {
        m_restart = count;
        return *this;
      }

      /**
       * @brief Solves the system on the device.
       * @param[in,out] system System, whose solution is the initial guess when
       * its size matches.
       */
      void solve(LinearSystemType& system) override
      {
        const auto& matrix = system.getOperator();
        const size_t n = matrix.getRowCount();
        if (matrix.getColumnCount() != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "GMRES requires a square matrix." << Alert::Raise;
        }
        if (static_cast<size_t>(system.getVector().size()) != n)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The right-hand side size does not match the matrix." << Alert::Raise;
        }
        if (m_restart == 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The restart dimension must be positive." << Alert::Raise;
        }
        if (n == 0)
        {
          m_result = {true, 0, 0};
          system.getSolution().resize(0);
          return;
        }
        if constexpr (requires { system.getSolution().getView(); })
        {
          if (static_cast<size_t>(system.getSolution().size()) != n)
            system.getSolution().setZero(n);
        }
        const typename Krylov::ExecutionSpace execution;
        Krylov::reserve(m_inverse, n, "Rodin::Solver::GMRES::inverse");
        if constexpr (!requires { system.getVector().getView(); })
          Krylov::reserve(m_b, n, "Rodin::Solver::GMRES::b");
        if constexpr (!requires { system.getSolution().getView(); })
          Krylov::reserve(m_x, n, "Rodin::Solver::GMRES::x");
        Krylov::upload(system.getVector(), m_b);
        if (static_cast<size_t>(system.getSolution().size()) == n)
          Krylov::upload(system.getSolution(), m_x);
        else
          ::Kokkos::deep_copy(m_x, Real(0));
        Krylov::invertDiagonal(execution, matrix, m_inverse);
        auto A = Krylov::wrap(matrix);

        const size_t restart = std::min(m_restart, n);
        const size_t limit = m_maxIterations ? *m_maxIterations : 2 * n;
        const size_t cycles = std::max<size_t>(1, (limit + restart - 1) / restart);
        Handle handle;
        handle.create_gmres_handle(static_cast<::Rodin::Math::SparseIndex>(restart),
          m_tolerance, static_cast<::Rodin::Math::SparseIndex>(cycles));
        Internal::KokkosJacobi<Backend> jacobi(m_inverse);
        ::KokkosSparse::Experimental::gmres(&handle, A, m_b, m_x, &jacobi);
        execution.fence();
        const auto* gmres = handle.get_gmres_handle();
        m_result.converged = gmres->get_conv_flag_val() == Handle::GMRESHandleType::Conv;
        m_result.iterations = static_cast<size_t>(gmres->get_num_iters());
        m_result.error = static_cast<Real>(gmres->get_end_rel_res());
        Krylov::download(m_x, system.getSolution());
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
       * @returns Iteration count reported by the GMRES handle.
       */
      size_t getIterationNumber() const
      {
        return m_result.iterations;
      }

      /**
       * @brief Relative residual of the last solve.
       * @returns Final relative residual reported by the GMRES handle.
       */
      Real getError() const
      {
        return m_result.error;
      }

      /**
       * @brief Creates a copy of this solver.
       * @returns Pointer to a new GMRES instance.
       */
      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }

    private:
      Real m_tolerance = std::numeric_limits<Real>::epsilon();
      Optional<size_t> m_maxIterations;
      size_t m_restart = 30;
      Internal::KokkosKrylovResult m_result;
      typename Krylov::Vector m_inverse, m_b, m_x;
  };
}

#endif
