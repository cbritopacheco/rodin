/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file RocSolverRF.h
 * @brief Sparse direct solver on AMD GPUs by rocSOLVER refactorization.
 *
 * rocSOLVER has no sparse factorization of its own. It refactorizes on the
 * GPU a matrix whose pivot order and factor pattern are known, as
 * @f$ PAQ = LU @f$ for a general matrix and as @f$ Q^TAQ = LL^T @f$ for a
 * symmetric positive definite one, and solves with the factors.
 *
 * The pivot order comes from the host: UMFPACK factorizes a general matrix
 * once with partial pivoting and returns @f$ P @f$ and @f$ Q @f$, and a
 * symmetric positive definite matrix is ordered by AMD. The pattern of the
 * factors is then computed symbolically from these permutations, without
 * dropping the coefficients that happen to vanish, so that every later
 * matrix with the same pattern is refactorized exactly.
 *
 * rocSOLVER reads the matrix in CSR format, whereas Math::SparseMatrix
 * stores it as CSC. A host matrix is therefore reordered as CSR on the host
 * before every upload. A Kokkos::Math::SparseMatrix is stored as CSR, and
 * rocSOLVER reads it in place when it is in device memory.
 *
 * @note This solver requires rocSOLVER, UMFPACK and RODIN_USE_ROCSOLVER.
 * Link Rodin::Kokkos and keep the Kokkos runtime initialized during use.
 */
#ifndef RODIN_KOKKOS_SOLVER_ROCSOLVERRF_H
#define RODIN_KOKKOS_SOLVER_ROCSOLVERRF_H

#include "Rodin/Configure.h"

#include <limits>
#include <Kokkos_Core.hpp>
#include "Rodin/Kokkos/Utility.h"

namespace Rodin::Kokkos::Solver
{
  /**
   * @brief Certifies that all factor pivots are finite and nonzero.
   *
   * ## Architecture
   * A Kokkos reduction counts invalid gathered pivots in their memory space.
   * Only that scalar count returns to the host. The same kernel runs with a
   * host execution space in tests and with HIP in rocSOLVER. No tolerance or
   * norm is used, so every finite nonzero pivot is accepted.
   */
  struct RocSolverRFPivots
  {
    /**
     * @brief Tests a pivot without introducing a numerical tolerance.
     * @param value Pivot value.
     * @return Whether the pivot is finite and nonzero.
     */
    KOKKOS_INLINE_FUNCTION static bool valid(double value)
    {
      constexpr double Maximum = std::numeric_limits<double>::max();
      return value != 0 && value >= -Maximum && value <= Maximum;
    }

    /**
     * @brief Checks every pivot on the supplied execution space.
     * @tparam Execution Execution space accessing the pivot view.
     * @tparam View One-dimensional view of factor pivots.
     * @param execution Execution instance on which the pivots are ready.
     * @param values Gathered pivots, which may be empty.
     * @return Whether every pivot is finite and nonzero, after reduction completion.
     */
    template <class Execution, class View>
    static bool check(const Execution& execution, const View& values)
    {
      static_assert(::Kokkos::SpaceAccessibility<Execution,
        typename View::memory_space>::accessible);
      size_t invalid = 0;
      ::Kokkos::parallel_reduce("Rodin rocSOLVER check pivots",
        ::Kokkos::RangePolicy<Execution, ::Kokkos::IndexType<size_t>>(
          execution, 0, values.extent(0)),
        Check<View>{values}, invalid);
      ::Rodin::Kokkos::Internal::checkLaunches("Rodin rocSOLVER check pivots");
      return invalid == 0;
    }

    /// @brief Counts invalid pivots without copying their values to the host.
    template <class View>
    struct Check
    {
      /// @brief Device pivots to inspect.
      View values;

      /**
       * @param i Pivot index.
       * @param invalid Number of invalid pivots in the reduction.
       */
      KOKKOS_INLINE_FUNCTION void operator()(size_t i, size_t& invalid) const
      {
        if (!valid(values(i)))
          ++invalid;
      }
    };
  };
}

#ifdef RODIN_USE_ROCSOLVER

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/OrderingMethods>

#include <hip/hip_runtime_api.h>
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>
#include <rocsparse/rocsparse.h>
#include <umfpack.h>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Vector.h"

#include "CuDSSOperator.h"
#include "CuDSSVector.h"
#include "Rodin/Solver/ForwardDecls.h"
#include "Rodin/Solver/Info.h"
#include "Rodin/Solver/LinearSolver.h"
#include "CompressedPattern.h"
#include "Rodin/Kokkos/Math/SparseMatrix.h"

namespace Rodin::FormLanguage
{
  template <class LinearSystem>
  struct Traits<Solver::RocSolverRF<LinearSystem>>
  {
    /// @brief Linear system type.
    using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup RocSolverRFSpecializations RocSolverRF Template Specializations
   * @brief Template specializations of the RocSolverRF class.
   * @see @ref RocSolverRF
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for RocSolverRF.
   */
  template <class LinearSystem>
  RocSolverRF(Variational::ProblemBase<LinearSystem>&) -> RocSolverRF<LinearSystem>;

  /**
   * @ingroup RocSolverRFSpecializations
   * @brief GPU direct solver for real sparse square systems on AMD GPUs.
   *
   * The first factorization of a sparsity pattern fixes the pivot order on
   * the host, computes the pattern of the factors, and hands both to the
   * rocSOLVER analysis. Every numeric factorization, the first included, is
   * then a rocSOLVER refactorization on the device, followed by device
   * solves. The matrix, the right-hand side and the solution are copied to
   * the device, and the solution is copied back. A matrix already in device
   * memory, such as a Kokkos::Math::SparseMatrix on the HIP backend, is read in
   * place and only copied to the host for the first factorization.
   *
   * The analysis, the pivot order included, is reused across explicit
   * @ref factorize calls while the sparsity pattern is unchanged; only the
   * values are copied again. A matrix with another pattern, such as a
   * tangent whose assembly drops other zero contributions, is analyzed
   * again. The pivots are therefore those of the first matrix of a pattern:
   * a later matrix whose values call for other pivots needs the analysis to
   * be cleared. The numeric factorization is retained after @ref factorize,
   * so several right-hand sides are solved without refactorization. The
   * ordinary @ref solve entry point builds a factorization only when none is
   * retained. Call @ref factorize explicitly after changing the matrix.
   *
   * Symmetry is never inferred from the entries. A matrix is factorized as
   * @f$ LL^T @f$ only after @ref setSymmetric, and only its lower triangle
   * is then read.
   *
   * @par Architecture
   * 1. Copy the compressed arrays of a host matrix to the device once per
   *    sparsity pattern, reordered as CSR, and the values once per
   *    factorization. Hand the arrays of a device matrix to rocSOLVER at
   *    every factorization, as the assembly may have moved them.
   * 2. Fix the pivot order on the host once per sparsity pattern, by
   *    UMFPACK for @f$ LU @f$ and by AMD for @f$ LL^T @f$, and compute the
   *    exact pattern of the factors from it.
   * 3. Refactorize on the device for every numeric factorization of the
   *    pattern. Kokkos graphs are retained and compared exactly on the backend;
   *    only a scalar mismatch result returns to the host.
   * 4. Retain the numeric factorization for explicit repeated solves.
   *
   * The solver runs on the current HIP device, on the default stream. Every
   * call returns once the device work it started has completed.
   *
   * The refactorization does not pivot, and rocSOLVER neither reports a zero
   * pivot nor fails on one: its triangular solves skip the division and
   * return a wrong solution. The factorization therefore gathers the
   * diagonal of the factors after every refactorization, checks it on the
   * device, and fails on a zero or non-finite pivot. The solve fails on a
   * non-finite solution. A
   * merely small pivot goes unnoticed, and a caller that needs to know
   * checks the residual. The Cholesky refactorization takes the square root
   * of the absolute value of every pivot, so an indefinite matrix declared
   * Symmetry::PositiveDefinite yields a wrong solution without failing.
   *
   * Info::status holds the @c rocblas_status of a failed rocSOLVER call, the
   * @c rocsparse_status of a failed pivot gather, the UMFPACK status of a
   * failed host factorization, and @ref Breakdown on a zero or non-finite
   * pivot or a non-finite solution. It is zero on success.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref RocSolverRF "RocSolverRF<Math::LinearSystem<Matrix, Math::Vector<Real>>>" | Real double-precision sparse systems whose matrix has a CuDSSOperator. |
   *
   * @tparam Matrix Math::SparseMatrix<Real> or Kokkos::Math::SparseMatrix<Real, Backend>.
   */
  template <class Matrix, class VectorType_>
  class RocSolverRF<Math::LinearSystem<Matrix, VectorType_>> final
    : public LinearSolverBase<Math::LinearSystem<Matrix, VectorType_>>
  {
    /// @brief Access to the compressed arrays of the matrix.
    using Arrays = CuDSSOperator<Matrix>;

    static_assert(!(Arrays::Resident && Arrays::Transposed),
      "A device matrix must hold the CSR arrays of the matrix itself.");
    static_assert(sizeof(Math::SparseIndex) == sizeof(rocblas_int));
    static_assert(std::is_same_v<Real, double>);

    public:
      using ScalarType = Real;
      using VectorType = VectorType_;
      using OperatorType = Matrix;
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /// @brief Factorization stage, shared by the factorization solvers.
      using Factorization = Rodin::Solver::Factorization;

      /**
       * @brief Info::status of a zero or non-finite pivot, or of a solution
       * that is not finite.
       */
      static constexpr Integer Breakdown = std::numeric_limits<Integer>::max();

      /// @brief Matrix symmetry declared to rocSOLVER.
      enum class Symmetry
      {
        /// @brief No symmetry; the matrix is factorized as @f$ LU @f$.
        Unsymmetric,
        /// @brief Symmetric positive definite; factorized as @f$ LL^T @f$.
        PositiveDefinite
      };

      /**
       * @brief Device array owned by the solver.
       * @tparam T Element type.
       */
      template <class T>
      struct Array
      {
        Array() = default;
        Array(const Array&) = delete;
        Array& operator=(const Array&) = delete;

        ~Array()
        {
          release();
        }

        /// @brief Device pointer, or null when nothing is allocated.
        T* data = nullptr;

        /// @brief Number of elements allocated.
        size_t size = 0;

        /**
         * @brief Allocates @p count elements, discarding the previous ones.
         * @returns The HIP status of the allocation.
         */
        hipError_t resize(size_t count) noexcept
        {
          if (count == size)
            return hipSuccess;
          release();
          if (count == 0)
            return hipSuccess;
          const hipError_t status = hipMalloc(&data, count * sizeof(T));
          if (status == hipSuccess)
            size = count;
          else
            data = nullptr;
          return status;
        }

        void release() noexcept
        {
          if (data)
            static_cast<void>(hipFree(data));
          data = nullptr;
          size = 0;
        }
      };

      /// @brief rocSOLVER objects and device arrays retained between solves.
      struct Resources
      {
        Resources() = default;
        Resources(const Resources&) = delete;
        Resources& operator=(const Resources&) = delete;

        /// @brief rocBLAS library handle, or null before the first use.
        rocblas_handle handle = nullptr;

        /// @brief rocSPARSE library handle, or null before the first use.
        rocsparse_handle sparse = nullptr;

        /// @brief Refactorization metadata, valid while @ref symbolic is true.
        rocsolver_rfinfo rfinfo = nullptr;

        /// @brief Row starts of the CSR arrays, copied from a host matrix.
        Array<Math::SparseIndex> outer;

        /// @brief Column indices of the CSR arrays, copied from a host matrix.
        Array<Math::SparseIndex> inner;

        /// @brief Values, copied from a host matrix.
        Array<ScalarType> values;

        /// @brief Row starts of the bundled factors.
        Array<Math::SparseIndex> factorStarts;

        /// @brief Column indices of the bundled factors.
        Array<Math::SparseIndex> factorIndices;

        /// @brief Values of the bundled factors.
        Array<ScalarType> factorValues;

        /// @brief Original row of every pivot row.
        Array<Math::SparseIndex> rowPivots;

        /// @brief Original column of every pivot column.
        Array<Math::SparseIndex> columnPivots;

        /// @brief Right-hand side, overwritten by the solution.
        Array<ScalarType> b;

        /// @brief Position in the bundled factors of every diagonal entry.
        Array<Math::SparseIndex> diagonal;

        /// @brief Diagonal of the factors, gathered to check the pivots.
        Array<ScalarType> pivots;

        /**
         * @brief Position in the arrays of a CSC host matrix of every CSR
         * coefficient; empty for CSR arrays.
         */
        std::vector<Math::SparseIndex> gather;

        /// @brief Values in CSR order, staged on the host for the upload.
        std::vector<ScalarType> staging;

        /// @brief Whether a reusable analysis is retained.
        bool symbolic = false;

        /// @brief Whether a reusable numeric factorization is retained.
        bool numeric = false;

        /// @brief Number of rows of the analyzed matrix.
        Index rows = 0;

        /// @brief Number of nonzeros of the analyzed matrix.
        Index nonzeros = 0;

        /**
         * @brief Owned snapshot of the analyzed graph, in the memory space
         * of the matrix.
         */
        typename Arrays::Pattern pattern;

        /// @brief Number of nonzeros of the bundled factors.
        Index factorNonzeros = 0;

        /**
         * @brief Releases a retained factorization stage.
         *
         * Clearing Factorization::Numeric only marks the factors as
         * unusable, because rocSOLVER overwrites them in place at the next
         * refactorization. Clearing Factorization::Symbolic also releases
         * the analysis and the device arrays, keeping the handle.
         */
        void clear(Factorization factorization) noexcept
        {
          numeric = false;
          if (factorization == Factorization::Numeric)
            return;
          symbolic = false;
          if (rfinfo)
            rocsolver_destroy_rfinfo(rfinfo);
          rfinfo = nullptr;
          outer.release();
          inner.release();
          values.release();
          factorStarts.release();
          factorIndices.release();
          factorValues.release();
          rowPivots.release();
          columnPivots.release();
          diagonal.release();
          pivots.release();
          b.release();
          gather.clear();
          gather.shrink_to_fit();
          staging.clear();
          staging.shrink_to_fit();
          pattern.clear();
          rows = 0;
          nonzeros = 0;
          factorNonzeros = 0;
        }

        ~Resources()
        {
          clear(Factorization::Symbolic);
          if (handle)
            rocblas_destroy_handle(handle);
          if (sparse)
            rocsparse_destroy_handle(sparse);
        }
      };

      /// @brief Constructs a rocSOLVER solver for the given problem.
      RocSolverRF(ProblemBaseType& pb)
        : Parent(pb)
      {}

      RocSolverRF(const RocSolverRF& other)
        : Parent(other),
          m_symmetry(other.m_symmetry)
      {}

      RocSolverRF(RocSolverRF&& other) noexcept
        : Parent(std::move(other)),
          m_symmetry(other.m_symmetry)
      {}

      ~RocSolverRF() override = default;

      /**
       * @brief Declares the symmetry of the system matrix.
       *
       * The caller guarantees that the matrix is symmetric positive definite
       * for Symmetry::PositiveDefinite, since rocSOLVER does not detect a
       * negative pivot. Symmetry selects the factorization itself and
       * therefore releases any retained factorization.
       *
       * @param symmetry Symmetry of the matrices to be factorized.
       * @returns Reference to this solver.
       */
      RocSolverRF& setSymmetric(Symmetry symmetry) noexcept
      {
        if (m_symmetry != symmetry)
        {
          m_symmetry = symmetry;
          clear(Factorization::Symbolic);
          m_info = Info{};
        }
        return *this;
      }

      /// @brief Returns the declared symmetry of the system matrix.
      Symmetry getSymmetric() const noexcept
      {
        return m_symmetry;
      }

      /**
       * @brief Analyzes and numerically factorizes the system matrix.
       *
       * Reuses the existing analysis and replaces any previous numeric
       * factorization. A matrix whose sparsity pattern differs from the
       * analyzed one is analyzed again.
       */
      void factorize(LinearSystemType& axb) override
      {
        auto& matrix = axb.getOperator();
        if (Arrays::getRowCount(matrix) != Arrays::getColumnCount(matrix))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "rocSOLVER requires a square matrix." << Alert::Raise;
        }

        Arrays::prepare(matrix);
        initialize();

        auto& res = m_resources;
        const Index rows = static_cast<Index>(Arrays::getRowCount(matrix));
        const Index nonzeros = static_cast<Index>(Arrays::getNonZeroCount(matrix));
        // The assemblies list only the nonzero contributions, so the pattern
        // of a tangent may change with the iterate. A new pattern is
        // analyzed again.
        if (res.symbolic && !hasPattern(matrix, rows, nonzeros))
          res.clear(Factorization::Symbolic);

        res.numeric = false;
        if (!res.symbolic)
        {
          if (!analyze(matrix, rows, nonzeros))
          {
            res.clear(Factorization::Symbolic);
            return;
          }
        }
        else if constexpr (!Arrays::Resident)
        {
          upload(matrix);
        }

        Math::SparseIndex* starts = res.outer.data;
        Math::SparseIndex* indices = res.inner.data;
        ScalarType* values = res.values.data;
        if constexpr (Arrays::Resident)
        {
          // rocSOLVER reads the arrays without writing them.
          starts = Arrays::getStarts(matrix);
          indices = Arrays::getIndices(matrix);
          values = Arrays::getValues(matrix);
        }

        rocblas_status status;
        if (m_symmetry == Symmetry::PositiveDefinite)
        {
          status =
            rocsolver_dcsrrf_refactchol(res.handle, rows, nonzeros, starts, indices,
              values, res.factorNonzeros, res.factorStarts.data, res.factorIndices.data,
              res.factorValues.data, res.columnPivots.data, res.rfinfo);
        }
        else
        {
          status = rocsolver_dcsrrf_refactlu(res.handle, rows, nonzeros, starts, indices,
            values, res.factorNonzeros, res.factorStarts.data, res.factorIndices.data,
            res.factorValues.data, res.rowPivots.data, res.columnPivots.data, res.rfinfo);
        }
        if (!complete(status, Factorization::Symbolic))
          return;
        if (!checkPivots())
          return;
        res.numeric = true;
        record(0, Factorization::Numeric);
      }

      /**
       * @brief Returns the outcome of the most recent operation.
       *
       * Updated by every factorization and every solve, so a caller reads it
       * after the call it wants to check.
       */
      const Info& getInfo() const noexcept
      {
        return m_info;
      }

      /**
       * @brief Checks whether the most recent operation succeeded.
       * @returns true if the solver succeeded, false otherwise.
       */
      Boolean success() const noexcept
      {
        return m_info.success;
      }

      /**
       * @brief Returns the retained rocSOLVER resources.
       *
       * A true @ref Resources::numeric denotes a reusable numeric
       * factorization, and a true @ref Resources::symbolic denotes a
       * reusable analysis.
       */
      const Resources& getResources() const noexcept
      {
        return m_resources;
      }

      /**
       * @brief Releases a retained factorization stage.
       *
       * Clearing Factorization::Numeric keeps the analysis, so a later
       * @ref factorize call skips it. Clearing Factorization::Symbolic
       * releases both stages and the device arrays.
       */
      void clear(Factorization factorization) noexcept
      {
        m_resources.clear(factorization);
        if (m_resources.symbolic)
          m_info.factorization = Factorization::Symbolic;
        else
          m_info.factorization.reset();
      }

      /**
       * @brief Solves the linear system with rocSOLVER.
       *
       * Reuses the retained numeric factorization and builds one only when
       * none is retained. The caller must call @ref factorize after changing
       * the matrix. A failed solve leaves the solution vector untouched.
       *
       * @param[in,out] axb System whose solution vector receives @f$ A^{-1}b @f$.
       */
      void solve(LinearSystemType& axb) override
      {
        if (static_cast<size_t>(axb.getVector().size()) !=
          Arrays::getRowCount(axb.getOperator()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The right-hand side size does not match the matrix." << Alert::Raise;
        }

        if (!m_resources.numeric)
        {
          factorize(axb);
          if (!success())
            return;
        }

        auto& res = m_resources;
        copy(res.b.data, axb.getVector().data(), res.rows,
          CuDSSVector<VectorType>::Resident ? hipMemcpyDeviceToDevice
                                                      : hipMemcpyHostToDevice);
        const rocblas_status status =
          rocsolver_dcsrrf_solve(res.handle, res.rows, 1, res.factorNonzeros,
            res.factorStarts.data, res.factorIndices.data, res.factorValues.data,
            res.rowPivots.data, res.columnPivots.data, res.b.data, res.rows, res.rfinfo);
        if (!complete(status, Factorization::Numeric))
          return;
        VectorType x(res.rows);
        copy(x.data(), res.b.data, res.rows,
          CuDSSVector<VectorType>::Resident ? hipMemcpyDeviceToDevice
                                                      : hipMemcpyDeviceToHost);
        if (!x.allFinite())
        {
          record(Breakdown, Factorization::Numeric);
          return;
        }
        axb.getSolution() = std::move(x);
        record(0, Factorization::Numeric);
      }

      RocSolverRF* copy() const noexcept override
      {
        return new RocSolverRF(*this);
      }

    private:
      /// @brief Creates the rocBLAS handle on first use.
      void initialize()
      {
        auto& res = m_resources;
        if (!res.handle)
          check(rocblas_create_handle(&res.handle), "create the rocBLAS handle");
        if (!res.sparse)
        {
          const rocsparse_status status = rocsparse_create_handle(&res.sparse);
          if (status != rocsparse_status_success)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "rocSPARSE failed to create its handle: "
              << rocsparse_get_status_name(status) << "." << Alert::Raise;
          }
        }
      }

      /**
       * @brief Gathers the diagonal of the factors and checks every pivot.
       *
       * The diagonal holds the pivots of @f$ U @f$ for @f$ LU @f$ and of
       * @f$ L @f$ for @f$ LL^T @f$.
       *
       * @returns Whether every pivot is finite and nonzero.
       */
      bool checkPivots()
      {
        auto& res = m_resources;
        const rocsparse_status status = rocsparse_dgthr(res.sparse,
          static_cast<rocsparse_int>(res.rows), res.factorValues.data, res.pivots.data,
          res.diagonal.data, rocsparse_index_base_zero);
        if (status != rocsparse_status_success)
        {
          record(static_cast<Integer>(status), Factorization::Symbolic);
          return false;
        }
        hipStream_t stream = nullptr;
        const rocsparse_status streamStatus = rocsparse_get_stream(res.sparse, &stream);
        if (streamStatus != rocsparse_status_success)
        {
          record(static_cast<Integer>(streamStatus), Factorization::Symbolic);
          return false;
        }
        using Pivots = ::Kokkos::View<const ScalarType*, ::Kokkos::HIPSpace,
          ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>;
        const ::Kokkos::HIP execution(stream);
        if (!::Rodin::Kokkos::Solver::RocSolverRFPivots::check(
              execution, Pivots(res.pivots.data, res.rows)))
        {
          record(Breakdown, Factorization::Symbolic);
          return false;
        }
        return true;
      }

      /**
       * @brief Fixes the pivot order and the pattern of the factors, and
       * runs the rocSOLVER analysis.
       * @param matrix Non-const owner borrowed for the vendor input arrays.
       * @param rows Number of rows.
       * @param nonzeros Number of stored coefficients.
       * @returns Whether the analysis succeeded.
       */
      bool analyze(OperatorType& matrix, Index rows, Index nonzeros)
      {
        auto& res = m_resources;

        // The host needs the arrays, wherever they are.
        std::vector<Math::SparseIndex> hostStarts;
        std::vector<Math::SparseIndex> hostIndices;
        std::vector<ScalarType> hostValues;
        const Math::SparseIndex* starts = Arrays::getStarts(matrix);
        const Math::SparseIndex* indices = Arrays::getIndices(matrix);
        const ScalarType* values = Arrays::getValues(matrix);
        if constexpr (Arrays::Resident)
        {
          hostStarts.resize(rows + 1);
          hostIndices.resize(nonzeros);
          hostValues.resize(nonzeros);
          copy(hostStarts.data(), starts, rows + 1, hipMemcpyDeviceToHost);
          copy(hostIndices.data(), indices, nonzeros, hipMemcpyDeviceToHost);
          copy(hostValues.data(), values, nonzeros, hipMemcpyDeviceToHost);
          starts = hostStarts.data();
          indices = hostIndices.data();
          values = hostValues.data();
        }

        // CSR arrays for rocSOLVER and CSC arrays for UMFPACK.
        ::Rodin::Kokkos::Solver::CompressedPattern other;
        std::vector<Math::SparseIndex> positions;
        other = ::Rodin::Kokkos::Solver::CompressedPattern::Transpose(
          rows, starts, indices, positions);
        std::vector<ScalarType> otherValues(nonzeros);
        for (Index k = 0; k < nonzeros; ++k)
          otherValues[k] = values[positions[k]];
        const Math::SparseIndex* csrStarts = starts;
        const Math::SparseIndex* csrIndices = indices;
        const Math::SparseIndex* cscStarts = other.starts.data();
        const Math::SparseIndex* cscIndices = other.indices.data();
        const ScalarType* cscValues = otherValues.data();
        if constexpr (Arrays::Transposed)
        {
          std::swap(csrStarts, cscStarts);
          std::swap(csrIndices, cscIndices);
          cscValues = values;
          res.gather = std::move(positions);
        }

        // Pivot order and pattern of the factors.
        std::vector<Math::SparseIndex> p(rows);
        std::vector<Math::SparseIndex> q(rows);
        ::Rodin::Kokkos::Solver::CompressedPattern factors;
        if (m_symmetry == Symmetry::PositiveDefinite)
        {
          // AMD orders the pattern of A + A^T, which it builds from a matrix
          // it owns rather than from a view.
          const Eigen::SparseMatrix<ScalarType, Eigen::ColMajor, Math::SparseIndex> a =
            Eigen::Map<
              const Eigen::SparseMatrix<ScalarType, Eigen::ColMajor, Math::SparseIndex>>(
              rows, rows, nonzeros, cscStarts, cscIndices, cscValues);
          Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, Math::SparseIndex>
            ordering;
          Eigen::AMDOrdering<Math::SparseIndex>()(a, ordering);
          for (Index k = 0; k < rows; ++k)
            q[k] = ordering.indices()[k];
          p = q;
          factors = ::Rodin::Kokkos::Solver::CompressedPattern::Cholesky(
            rows, cscStarts, cscIndices, q);
        }
        else
        {
          const int status = pivot(rows, cscStarts, cscIndices, cscValues, p, q);
          if (status != UMFPACK_OK)
          {
            record(static_cast<Integer>(status), Optional<Factorization>{});
            return false;
          }
          factors = ::Rodin::Kokkos::Solver::CompressedPattern::LU(
            rows, cscStarts, cscIndices, p, q);
        }

        // Device arrays.
        const Index factorNonzeros = static_cast<Index>(factors.indices.size());
        check(res.b.resize(rows), "allocate the right-hand side on the device");
        check(res.factorStarts.resize(rows + 1), "allocate the factors on the device");
        check(
          res.factorIndices.resize(factorNonzeros), "allocate the factors on the device");
        check(
          res.factorValues.resize(factorNonzeros), "allocate the factors on the device");
        check(res.rowPivots.resize(rows), "allocate the pivots on the device");
        check(res.columnPivots.resize(rows), "allocate the pivots on the device");
        copy(
          res.factorStarts.data, factors.starts.data(), rows + 1, hipMemcpyHostToDevice);
        copy(res.factorIndices.data, factors.indices.data(), factorNonzeros,
          hipMemcpyHostToDevice);
        check(hipMemset(res.factorValues.data, 0, factorNonzeros * sizeof(ScalarType)),
          "clear the factors");
        copy(res.rowPivots.data, p.data(), rows, hipMemcpyHostToDevice);
        copy(res.columnPivots.data, q.data(), rows, hipMemcpyHostToDevice);
        res.rows = rows;
        res.nonzeros = nonzeros;
        res.factorNonzeros = factorNonzeros;

        // Every row of the factors holds its diagonal entry.
        std::vector<Math::SparseIndex> diagonal(rows);
        for (Index i = 0; i < rows; ++i)
        {
          const auto first = factors.indices.begin() + factors.starts[i];
          const auto last = factors.indices.begin() + factors.starts[i + 1];
          diagonal[i] = static_cast<Math::SparseIndex>(
            std::lower_bound(first, last, static_cast<Math::SparseIndex>(i)) -
            factors.indices.begin());
        }
        check(res.diagonal.resize(rows), "allocate the pivots on the device");
        check(res.pivots.resize(rows), "allocate the pivots on the device");
        copy(res.diagonal.data, diagonal.data(), rows, hipMemcpyHostToDevice);

        Math::SparseIndex* deviceStarts;
        Math::SparseIndex* deviceIndices;
        ScalarType* deviceValues;
        if constexpr (Arrays::Resident)
        {
          deviceStarts = Arrays::getStarts(matrix);
          deviceIndices = Arrays::getIndices(matrix);
          deviceValues = Arrays::getValues(matrix);
        }
        else
        {
          check(res.outer.resize(rows + 1), "allocate the matrix on the device");
          check(res.inner.resize(nonzeros), "allocate the matrix on the device");
          check(res.values.resize(nonzeros), "allocate the matrix on the device");
          copy(res.outer.data, csrStarts, rows + 1, hipMemcpyHostToDevice);
          copy(res.inner.data, csrIndices, nonzeros, hipMemcpyHostToDevice);
          upload(matrix);
          deviceStarts = res.outer.data;
          deviceIndices = res.inner.data;
          deviceValues = res.values.data;
        }

        check(rocsolver_create_rfinfo(&res.rfinfo, res.handle),
          "create the rocSOLVER metadata");
        check(rocsolver_set_rfinfo_mode(res.rfinfo,
                m_symmetry == Symmetry::PositiveDefinite ? rocsolver_rfinfo_mode_cholesky
                                                         : rocsolver_rfinfo_mode_lu),
          "set the rocSOLVER factorization mode");
        const rocblas_status status = rocsolver_dcsrrf_analysis(res.handle, rows, 1,
          nonzeros, deviceStarts, deviceIndices, deviceValues, factorNonzeros,
          res.factorStarts.data, res.factorIndices.data, res.factorValues.data,
          res.rowPivots.data, res.columnPivots.data, res.b.data, rows, res.rfinfo);
        if (!complete(status, Optional<Factorization>{}))
          return false;
        res.symbolic = true;
        res.pattern.assign(matrix);
        return true;
      }

      /// @brief Returns whether the matrix has the analyzed sparsity pattern.
      bool hasPattern(const OperatorType& matrix, Index rows, Index nonzeros) const
      {
        const auto& res = m_resources;
        return rows == res.rows && nonzeros == res.nonzeros &&
          res.pattern.matches(matrix);
      }

      /**
       * @brief Factorizes the CSC matrix once with UMFPACK for its pivots.
       * @returns The UMFPACK status, with a singular matrix as a failure.
       */
      static int pivot(Index rows, const Math::SparseIndex* starts,
        const Math::SparseIndex* indices, const ScalarType* values,
        std::vector<Math::SparseIndex>& p, std::vector<Math::SparseIndex>& q)
      {
        static_assert(sizeof(Math::SparseIndex) == sizeof(int));
        const int n = static_cast<int>(rows);
        double control[UMFPACK_CONTROL];
        umfpack_di_defaults(control);
        void* symbolic = nullptr;
        void* numeric = nullptr;
        int status =
          umfpack_di_symbolic(n, n, starts, indices, values, &symbolic, control, nullptr);
        if (status == UMFPACK_OK)
        {
          status = umfpack_di_numeric(
            starts, indices, values, symbolic, &numeric, control, nullptr);
        }
        if (symbolic)
          umfpack_di_free_symbolic(&symbolic);
        if (status == UMFPACK_OK)
        {
          status = umfpack_di_get_numeric(nullptr, nullptr, nullptr, nullptr, nullptr,
            nullptr, p.data(), q.data(), nullptr, nullptr, nullptr, numeric);
        }
        if (numeric)
          umfpack_di_free_numeric(&numeric);
        return status;
      }

      /// @brief Copies the values of a host matrix to the device, as CSR.
      void upload(const OperatorType& matrix)
      {
        auto& res = m_resources;
        const ScalarType* values = Arrays::getValues(matrix);
        if (!res.gather.empty())
        {
          res.staging.resize(res.nonzeros);
          for (Index k = 0; k < res.nonzeros; ++k)
            res.staging[k] = values[res.gather[k]];
          values = res.staging.data();
        }
        copy(res.values.data, values, res.nonzeros, hipMemcpyHostToDevice);
      }

      template <class T>
      void copy(T* destination, const T* source, Index count, hipMemcpyKind kind)
      {
        check(hipMemcpy(destination, source, count * sizeof(T), kind), "copy the data");
      }

      /**
       * @brief Waits for a rocSOLVER call and records its failure.
       * @param status Status of the call.
       * @param survived Factorization stage retained if the call failed.
       * @returns Whether the call succeeded.
       */
      bool complete(rocblas_status status, Optional<Factorization> survived)
      {
        if (status == rocblas_status_success)
          check(hipStreamSynchronize(nullptr), "wait for the device");
        else
        {
          record(static_cast<Integer>(status), survived);
          return false;
        }
        return true;
      }

      void check(hipError_t status, StringView operation) const
      {
        if (status != hipSuccess)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "HIP failed to " << operation << ": " << hipGetErrorString(status) << "."
            << Alert::Raise;
        }
      }

      void check(rocblas_status status, StringView operation) const
      {
        if (status != rocblas_status_success)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "rocSOLVER failed to " << operation << ": "
            << rocblas_status_to_string(status) << "." << Alert::Raise;
        }
      }

      /**
       * @brief Records the status and the stage that survived it.
       *
       * A failing analysis, factorization or solve is an outcome rather than
       * a defect in the caller, so it is reported through getInfo() instead
       * of raising.
       */
      void record(Integer status, Optional<Factorization> factorization)
      {
        m_info.status = status;
        m_info.success = status == 0;
        m_info.factorization = factorization;
        if (!m_info.success && factorization != Factorization::Numeric)
          m_resources.numeric = false;
      }

      /// @brief Outcome of the most recent operation.
      Info m_info;

      Symmetry m_symmetry = Symmetry::Unsymmetric;
      Resources m_resources;
  };
}

#endif
#endif
