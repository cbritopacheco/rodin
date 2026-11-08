/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_KRYLOV_H
#define RODIN_KOKKOS_SOLVER_KRYLOV_H

/**
 * @file
 * @brief Device storage and kernels shared by the Kokkos-Kernels Krylov
 * solvers.
 */

#include <string>

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_spmv.hpp>

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"

#include "Rodin/Kokkos/SparseMatrix.h"

namespace Rodin::Solver::Internal
{
  /**
   * @brief Device vectors and kernels of a Jacobi preconditioned Krylov
   * solver on a Kokkos::SparseMatrix.
   *
   * The matrix is wrapped without copy. The right-hand side and the solution
   * are copied between the host vectors of the system and device vectors,
   * which are kept from one solve to the next while the size is unchanged.
   *
   * @tparam Backend Kokkos backend of the matrix.
   */
  template <class Backend>
  struct KokkosKrylov
  {
    /// @brief Matrix type.
    using MatrixType = ::Rodin::Kokkos::SparseMatrix<Real, Backend>;

    /// @brief Execution space of the kernels.
    using ExecutionSpace = typename MatrixType::ExecutionSpace;

    /// @brief Memory space of the vectors.
    using MemorySpace = typename MatrixType::MemorySpace;

    /// @brief Kokkos device of the vectors.
    using Device = ::Kokkos::Device<ExecutionSpace, MemorySpace>;

    /// @brief Device vector.
    using Vector = ::Kokkos::View<Real*, Device>;

    /// @brief Kokkos-Kernels view of the matrix.
    using CrsMatrix = ::KokkosSparse::CrsMatrix<Real, Math::SparseIndex, Device, void,
      Math::SparseIndex>;

    /**
     * @brief Kokkos-Kernels SpMV handle, which keeps the cuSPARSE descriptors
     * of one matrix between products.
     */
    using SpMV = ::KokkosSparse::SPMVHandle<Device, CrsMatrix, Vector, Vector>;

    /**
     * @brief Wraps the arrays of a matrix.
     * @param matrix Square matrix.
     * @returns Kokkos-Kernels matrix sharing the arrays.
     */
    static CrsMatrix wrap(const MatrixType& matrix)
    {
      const auto n = static_cast<Math::SparseIndex>(matrix.getRowCount());
      return CrsMatrix("Rodin::Solver::Krylov::A", n, n,
        static_cast<Math::SparseIndex>(matrix.getNonZeroCount()), matrix.getValues(),
        matrix.getRowStarts(), matrix.getColumns());
    }

    /**
     * @brief Resizes a vector, keeping it when the size is unchanged.
     * @param v Vector.
     * @param n Size.
     * @param label Label of a new allocation.
     */
    static void reserve(Vector& v, size_t n, const char* label)
    {
      if (v.extent(0) != n)
      {
        v = Vector(
          ::Kokkos::view_alloc(::Kokkos::WithoutInitializing, std::string(label)), n);
      }
    }

    /**
     * @brief Copies a host vector to the device.
     * @param host Host vector.
     * @param device Device vector of the same size.
     */
    static void upload(const Math::Vector<Real>& host, const Vector& device)
    {
      const ::Kokkos::View<const Real*, ::Kokkos::HostSpace, ::Kokkos::MemoryUnmanaged>
        h(host.data(), static_cast<size_t>(host.size()));
      ::Kokkos::deep_copy(device, h);
    }

    /**
     * @brief Copies a device vector to the host.
     * @param device Device vector.
     * @param host Host vector, resized.
     */
    static void download(const Vector& device, Math::Vector<Real>& host)
    {
      host.resize(static_cast<Eigen::Index>(device.extent(0)));
      const ::Kokkos::View<Real*, ::Kokkos::HostSpace, ::Kokkos::MemoryUnmanaged> h(
        host.data(), device.extent(0));
      ::Kokkos::deep_copy(h, device);
    }

    /// @brief Writes the inverse of the diagonal coefficient of each row.
    struct InvertDiagonal
    {
      /// @brief CSR row offsets of the matrix.
      ::Kokkos::View<Math::SparseIndex*, MemorySpace> starts;

      /// @brief CSR column indices of the matrix.
      ::Kokkos::View<Math::SparseIndex*, MemorySpace> columns;

      /// @brief CSR values of the matrix.
      ::Kokkos::View<Real*, MemorySpace> values;

      /// @brief Inverse diagonal, set to one where the diagonal is zero or absent.
      Vector inverse;

      /**
       * @brief Inverts the diagonal coefficient of one row.
       * @param i Row index.
       */
      KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
      {
        Real d = 0;
        for (auto k = starts(i); k < starts(i + 1); ++k)
        {
          if (static_cast<size_t>(columns(k)) == i)
            d = values(k);
        }
        inverse(i) = d != Real(0) ? Real(1) / d : Real(1);
      }
    };

    /// @brief Sums the products of the coefficients of two vectors.
    struct Dot
    {
      /// @brief First vector.
      Vector a;

      /// @brief Second vector.
      Vector b;

      /**
       * @brief Accumulates the product of one pair of coefficients.
       * @param i Coefficient index.
       * @param sum Reduction accumulator.
       */
      KOKKOS_INLINE_FUNCTION void operator()(const size_t i, Real& sum) const
      {
        sum += a(i) * b(i);
      }
    };

    /// @brief Writes @f$ y = \beta y + \alpha D^{-1} x @f$ coefficientwise.
    struct Scale
    {
      /// @brief Inverse diagonal @f$ D^{-1} @f$.
      Vector inverse;

      /// @brief Input vector @f$ x @f$.
      ::Kokkos::View<const Real*, Device> x;

      /// @brief Output vector @f$ y @f$.
      Vector y;

      /// @brief Factor @f$ \alpha @f$ of @f$ D^{-1} x @f$.
      Real alpha;

      /// @brief Factor @f$ \beta @f$ of @f$ y @f$, zero to overwrite @f$ y @f$.
      Real beta;

      /**
       * @brief Writes one coefficient of @f$ y @f$.
       * @param i Coefficient index.
       */
      KOKKOS_INLINE_FUNCTION void operator()(const size_t i) const
      {
        const Real v = alpha * inverse(i) * x(i);
        y(i) = beta == Real(0) ? v : beta * y(i) + v;
      }
    };

    /**
     * @brief Inverts the diagonal of a matrix, as Eigen's
     * DiagonalPreconditioner: a zero or missing diagonal coefficient gives one.
     * @param execution Execution space instance.
     * @param matrix Square matrix.
     * @param inverse Inverse of the diagonal, of the size of the matrix.
     */
    static void invertDiagonal(
      const ExecutionSpace& execution, const MatrixType& matrix, const Vector& inverse)
    {
      const auto starts = matrix.getRowStarts();
      const auto columns = matrix.getColumns();
      const auto values = matrix.getValues();
      ::Kokkos::parallel_for("Rodin::Solver::Krylov::invertDiagonal",
        ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, inverse.extent(0)),
        InvertDiagonal{starts, columns, values, inverse});
    }

    /**
     * @brief Dot product.
     * @param execution Execution space instance.
     * @param a First vector.
     * @param b Second vector.
     * @returns @f$ a \cdot b @f$, after the kernel completed.
     */
    static Real dot(const ExecutionSpace& execution, const Vector& a, const Vector& b)
    {
      Real res = 0;
      ::Kokkos::parallel_reduce("Rodin::Solver::Krylov::dot",
        ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, a.extent(0)), Dot{a, b},
        res);
      return res;
    }

    /**
     * @brief Applies the inverse diagonal: @f$ y = \beta y + \alpha D^{-1} x @f$.
     * @param execution Execution space instance.
     * @param inverse Inverse of the diagonal.
     * @param x Input vector.
     * @param y Output vector.
     * @param alpha Coefficient of @f$ D^{-1} x @f$.
     * @param beta Coefficient of @f$ y @f$, which is not read when zero.
     */
    static void scale(const ExecutionSpace& execution, const Vector& inverse,
      const ::Kokkos::View<const Real*, Device>& x, const Vector& y, Real alpha,
      Real beta)
    {
      ::Kokkos::parallel_for("Rodin::Solver::Krylov::scale",
        ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, y.extent(0)),
        Scale{inverse, x, y, alpha, beta});
    }
  };

  /// @brief Outcome of a Krylov solve on the device.
  struct KokkosKrylovResult
  {
    /// @brief Whether the tolerance was reached.
    bool converged = false;

    /// @brief Number of iterations.
    size_t iterations = 0;

    /// @brief Relative residual estimated by the method.
    Real error = 0;
  };
}

#endif
