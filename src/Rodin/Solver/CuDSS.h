/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CuDSS.h
 * @brief NVIDIA cuDSS sparse direct solver wrapper.
 *
 * cuDSS factorizes a sparse matrix on an NVIDIA GPU, as @f$ LU @f$ for a
 * general matrix and as @f$ LDL^T @f$ or @f$ LL^T @f$ for a symmetric one.
 *
 * cuDSS reads the matrix in CSR format, whereas Math::SparseMatrix stores
 * it as CSC. The CSC arrays of @f$ A @f$ are the CSR arrays of @f$ A^T @f$,
 * so they are copied to the device unchanged and cuDSS is asked to solve
 * with the transpose of the matrix it holds. A symmetric matrix needs no
 * transpose. A Kokkos::SparseMatrix is stored as CSR, and cuDSS reads it in
 * place when it is in device memory.
 *
 * @note This solver requires cuDSS and RODIN_USE_CUDSS.
 */
#ifndef RODIN_SOLVER_CUDSS_H
#define RODIN_SOLVER_CUDSS_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_CUDSS

#include <cstdint>
#include <limits>

#include <cuda_runtime.h>
#include <cudss.h>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Vector.h"

#include "ForwardDecls.h"
#include "Info.h"
#include "LinearSolver.h"

namespace Rodin::FormLanguage
{
  template <class LinearSystem>
  struct Traits<Solver::CuDSS<LinearSystem>>
  {
    /// @brief Linear system type.
    using LinearSystemType = LinearSystem;
  };
}

namespace Rodin::Solver::Internal
{
  /// @brief Host CSC arrays, read as the CSR arrays of the transpose.
  template <>
  struct CuDSSOperator<Math::SparseMatrix<Real>>
  {
    /// @brief Matrix type.
    using MatrixType = Math::SparseMatrix<Real>;

    /// @brief Whether the arrays are those of the transpose.
    static constexpr bool Transposed = true;

    /// @brief Whether the arrays are in device memory.
    static constexpr bool Resident = false;

    /// @brief Compresses the matrix.
    static void prepare(MatrixType& matrix)
    {
      matrix.makeCompressed();
    }

    /// @brief Number of rows.
    static size_t getRowCount(const MatrixType& matrix)
    {
      return static_cast<size_t>(matrix.rows());
    }

    /// @brief Number of columns.
    static size_t getColumnCount(const MatrixType& matrix)
    {
      return static_cast<size_t>(matrix.cols());
    }

    /// @brief Number of stored coefficients.
    static size_t getNonZeroCount(const MatrixType& matrix)
    {
      return static_cast<size_t>(matrix.nonZeros());
    }

    /// @brief Start of every compressed row, followed by the coefficient count.
    static const Math::SparseIndex* getStarts(const MatrixType& matrix)
    {
      return matrix.outerIndexPtr();
    }

    /// @brief Index of every coefficient in its compressed row.
    static const Math::SparseIndex* getIndices(const MatrixType& matrix)
    {
      return matrix.innerIndexPtr();
    }

    /// @brief Value of every coefficient.
    static const Real* getValues(const MatrixType& matrix)
    {
      return matrix.valuePtr();
    }
  };
}

namespace Rodin::Solver
{
  /**
   * @defgroup CuDSSSpecializations CuDSS Template Specializations
   * @brief Template specializations of the CuDSS class.
   * @see @ref CuDSS
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD (Class Template Argument Deduction) guide for CuDSS.
   */
  template <class LinearSystem>
  CuDSS(Variational::ProblemBase<LinearSystem>&) -> CuDSS<LinearSystem>;

  /**
   * @ingroup CuDSSSpecializations
   * @brief GPU direct solver for real sparse square systems.
   *
   * The matrix, the right-hand side and the solution are copied to the
   * device, factorized and solved there, and the solution is copied back. A
   * matrix already in device memory, such as a Kokkos::SparseMatrix on the
   * CUDA backend, is read in place and only the vectors are copied.
   * Symbolic analysis (reordering and symbolic factorization) is reused
   * across explicit @ref factorize calls, which assume that the sparsity
   * pattern is unchanged; only the values are copied again. The numeric
   * factorization is retained after @ref factorize, so several right-hand
   * sides are solved without refactorization. The ordinary @ref solve entry
   * point builds a factorization only when none is retained. Call
   * @ref factorize explicitly after changing matrix values, and clear the
   * symbolic factorization after changing the sparsity pattern.
   *
   * Symmetry is never inferred from the entries. A matrix is factorized as
   * symmetric only after @ref setSymmetric.
   *
   * @par Architecture
   * 1. Copy the compressed arrays of a host matrix to the device once per
   *    sparsity pattern, and the values once per factorization. Hand the
   *    arrays of a device matrix to cuDSS at every factorization, as the
   *    assembly may have moved them.
   * 2. Hand CSC arrays to cuDSS as the CSR arrays of @f$ A^T @f$ and solve
   *    with the transpose, so that no conversion is needed.
   * 3. Reuse symbolic analysis until explicitly cleared.
   * 4. Retain the numeric factorization for explicit repeated solves.
   *
   * The solver runs on the current CUDA device, on the default stream. Every
   * call returns once the device work it started has completed.
   *
   * cuDSS replaces a pivot that is too small by a small perturbation instead
   * of stopping, so a singular matrix is not always reported as a failure.
   * A caller that needs to know checks the residual.
   *
   * Info::status holds the @c cudssStatus_t of a failed cuDSS call. When the
   * calls succeed, it holds cuDSS's @c CUDSS_DATA_INFO, which is zero on
   * success and positive when the factorization met a zero pivot.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref CuDSS "CuDSS<Math::LinearSystem<Matrix, Math::Vector<Real>>>" | Real double-precision sparse systems whose matrix has an Internal::CuDSSOperator. |
   *
   * @tparam Matrix Math::SparseMatrix<Real> or Kokkos::SparseMatrix<Real, Backend>.
   */
  template <class Matrix>
  class CuDSS<Math::LinearSystem<Matrix, Math::Vector<Real>>> final
    : public LinearSolverBase<Math::LinearSystem<Matrix, Math::Vector<Real>>>
  {
    /// @brief Access to the compressed arrays of the matrix.
    using Arrays = Internal::CuDSSOperator<Matrix>;

    public:
      using ScalarType = Real;
      using VectorType = Math::Vector<ScalarType>;
      using OperatorType = Matrix;
      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      using Parent = LinearSolverBase<LinearSystemType>;

      using Parent::solve;

      /// @brief Factorization stage, shared by the factorization solvers.
      using Factorization = Rodin::Solver::Factorization;

      /// @brief Matrix symmetry declared to cuDSS.
      enum class Symmetry
      {
        /// @brief No symmetry; the matrix is factorized as @f$ LU @f$.
        Unsymmetric,
        /// @brief Symmetric, possibly indefinite; factorized as @f$ LDL^T @f$.
        General,
        /// @brief Symmetric positive definite; factorized as @f$ LL^T @f$.
        PositiveDefinite
      };

      /**
       * @brief Device array owned by the solver.
       * @tparam T Element type.
       */
      template <class T>
      struct DeviceArray
      {
        DeviceArray() = default;
        DeviceArray(const DeviceArray&) = delete;
        DeviceArray& operator=(const DeviceArray&) = delete;

        ~DeviceArray()
        {
          release();
        }

        /// @brief Device pointer, or null when nothing is allocated.
        T* data = nullptr;

        /// @brief Number of elements allocated.
        size_t size = 0;

        /**
         * @brief Allocates @p count elements, discarding the previous ones.
         * @returns The CUDA status of the allocation.
         */
        cudaError_t resize(size_t count) noexcept
        {
          if (count == size)
            return cudaSuccess;
          release();
          if (count == 0)
            return cudaSuccess;
          const cudaError_t status = cudaMalloc(&data, count * sizeof(T));
          if (status == cudaSuccess)
            size = count;
          else
            data = nullptr;
          return status;
        }

        void release() noexcept
        {
          if (data)
            cudaFree(data);
          data = nullptr;
          size = 0;
        }
      };

      /// @brief cuDSS objects and device arrays retained between solves.
      struct Resources
      {
        Resources() = default;
        Resources(const Resources&) = delete;
        Resources& operator=(const Resources&) = delete;

        /// @brief cuDSS library handle, or null before the first use.
        cudssHandle_t handle = nullptr;

        /// @brief Solver settings.
        cudssConfig_t config = nullptr;

        /// @brief Factorization data, valid while @ref symbolic is true.
        cudssData_t data = nullptr;

        /// @brief Matrix descriptor over the device CSR arrays.
        cudssMatrix_t matrix = nullptr;

        /// @brief Dense descriptor over the device right-hand side.
        cudssMatrix_t rhs = nullptr;

        /// @brief Dense descriptor over the device solution.
        cudssMatrix_t solution = nullptr;

        /// @brief Row starts of the CSR arrays, copied from a host matrix.
        DeviceArray<Math::SparseIndex> outer;

        /// @brief Column indices of the CSR arrays, copied from a host matrix.
        DeviceArray<Math::SparseIndex> inner;

        /// @brief Values, copied from a host matrix.
        DeviceArray<ScalarType> values;

        /// @brief Right-hand side on the device.
        DeviceArray<ScalarType> b;

        /// @brief Solution on the device.
        DeviceArray<ScalarType> x;

        /// @brief Whether reusable symbolic analysis is retained.
        bool symbolic = false;

        /// @brief Whether a reusable numeric factorization is retained.
        bool numeric = false;

        /// @brief Number of rows of the analyzed matrix.
        Index rows = 0;

        /// @brief Number of nonzeros of the analyzed matrix.
        Index nonzeros = 0;

        /**
         * @brief Releases a retained factorization stage.
         *
         * Clearing Factorization::Numeric only marks the factors as
         * unusable, because cuDSS overwrites them in place at the next
         * factorization. Clearing Factorization::Symbolic also releases the
         * factorization data and the device arrays, keeping the handle.
         */
        void clear(Factorization factorization) noexcept
        {
          numeric = false;
          if (factorization == Factorization::Numeric)
            return;
          symbolic = false;
          destroyDescriptors();
          if (data)
            cudssDataDestroy(handle, data);
          data = nullptr;
          outer.release();
          inner.release();
          values.release();
          b.release();
          x.release();
          rows = 0;
          nonzeros = 0;
        }

        void destroyDescriptors() noexcept
        {
          if (matrix)
            cudssMatrixDestroy(matrix);
          if (rhs)
            cudssMatrixDestroy(rhs);
          if (solution)
            cudssMatrixDestroy(solution);
          matrix = nullptr;
          rhs = nullptr;
          solution = nullptr;
        }

        ~Resources()
        {
          clear(Factorization::Symbolic);
          if (config)
            cudssConfigDestroy(config);
          if (handle)
            cudssDestroy(handle);
        }
      };

      /// @brief Constructs a cuDSS solver for the given problem.
      CuDSS(ProblemBaseType& pb)
        : Parent(pb)
      {}

      CuDSS(const CuDSS& other)
        : Parent(other),
          m_symmetry(other.m_symmetry),
          m_refinementSteps(other.m_refinementSteps)
      {}

      CuDSS(CuDSS&& other) noexcept
        : Parent(std::move(other)),
          m_symmetry(other.m_symmetry),
          m_refinementSteps(other.m_refinementSteps)
      {}

      ~CuDSS() override = default;

      /**
       * @brief Declares the symmetry of the system matrix.
       *
       * The caller guarantees that the matrix is symmetric for
       * Symmetry::General and Symmetry::PositiveDefinite. Symmetry selects
       * the factorization itself and therefore releases any retained
       * factorization.
       *
       * @param symmetry Symmetry of the matrices to be factorized.
       * @returns Reference to this solver.
       */
      CuDSS& setSymmetric(Symmetry symmetry) noexcept
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
       * @brief Sets the number of iterative refinement steps of each solve.
       *
       * cuDSS refines the solution on the device with the retained factors.
       * The default is zero.
       *
       * @param steps Number of refinement steps.
       * @returns Reference to this solver.
       */
      CuDSS& setRefinementSteps(size_t steps)
      {
        if (steps > static_cast<size_t>(std::numeric_limits<int>::max()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The refinement step count exceeds the cuDSS integer range."
            << Alert::Raise;
        }
        m_refinementSteps = steps;
        return *this;
      }

      /// @brief Returns the configured number of refinement steps.
      size_t getRefinementSteps() const noexcept
      {
        return m_refinementSteps;
      }

      /**
       * @brief Analyzes and numerically factorizes the system matrix.
       *
       * Reuses existing symbolic analysis and replaces any previous numeric
       * factorization. The matrix sparsity pattern must remain unchanged
       * until the symbolic factorization is cleared.
       */
      void factorize(LinearSystemType& axb)
      {
        auto& matrix = axb.getOperator();
        if (Arrays::getRowCount(matrix) != Arrays::getColumnCount(matrix))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "cuDSS requires a square matrix." << Alert::Raise;
        }

        Arrays::prepare(matrix);
        initialize();

        auto& res = m_resources;
        const Index rows = static_cast<Index>(Arrays::getRowCount(matrix));
        const Index nonzeros = static_cast<Index>(Arrays::getNonZeroCount(matrix));
        if (res.symbolic && (rows != res.rows || nonzeros != res.nonzeros))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "The matrix sparsity pattern changed since the symbolic analysis."
            << Alert::Raise;
        }

        res.numeric = false;
        const Math::SparseIndex* starts = res.outer.data;
        const Math::SparseIndex* indices = res.inner.data;
        const ScalarType* values = res.values.data;
        if constexpr (Arrays::Resident)
        {
          if (!res.symbolic)
            allocate(rows);
          starts = Arrays::getStarts(matrix);
          indices = Arrays::getIndices(matrix);
          values = Arrays::getValues(matrix);
        }
        else
        {
          if (!res.symbolic)
          {
            allocate(rows);
            check(res.outer.resize(rows + 1), "allocate the matrix on the device");
            check(res.inner.resize(nonzeros), "allocate the matrix on the device");
            check(res.values.resize(nonzeros), "allocate the matrix on the device");
            upload(res.outer.data, Arrays::getStarts(matrix), rows + 1);
            upload(res.inner.data, Arrays::getIndices(matrix), nonzeros);
            starts = res.outer.data;
            indices = res.inner.data;
            values = res.values.data;
          }
          upload(res.values.data, Arrays::getValues(matrix), nonzeros);
        }

        configure();
        if (!res.symbolic)
        {
          if (!describe(rows, nonzeros, starts, indices, values) ||
            !execute(CUDSS_PHASE_ANALYSIS, Optional<Factorization>{}))
          {
            res.clear(Factorization::Symbolic);
            return;
          }
          res.symbolic = true;
          res.rows = rows;
          res.nonzeros = nonzeros;
        }
        else if constexpr (Arrays::Resident)
        {
          check(cudssMatrixSetCsrPointers(res.matrix, starts, nullptr, indices, values),
            "update the matrix arrays");
        }

        if (!execute(CUDSS_PHASE_FACTORIZATION, Factorization::Symbolic))
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
       * @brief Returns the retained cuDSS resources.
       *
       * A true @ref Resources::numeric denotes a reusable numeric
       * factorization, and a true @ref Resources::symbolic denotes reusable
       * symbolic analysis.
       */
      const Resources& getResources() const noexcept
      {
        return m_resources;
      }

      /**
       * @brief Releases a retained factorization stage.
       *
       * Clearing Factorization::Numeric keeps the symbolic analysis, so a
       * later @ref factorize call skips it. Clearing Factorization::Symbolic
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
       * @brief Solves the linear system with cuDSS.
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
        configure();
        upload(res.b.data, axb.getVector().data(), res.rows);
        if (!execute(CUDSS_PHASE_SOLVE, Factorization::Numeric))
          return;
        VectorType& x = axb.getSolution();
        x.resize(res.rows);
        check(cudaMemcpy(x.data(), res.x.data, res.rows * sizeof(ScalarType),
                cudaMemcpyDeviceToHost),
          "copy the solution to the host");
        record(0, Factorization::Numeric);
      }

      CuDSS* copy() const noexcept override
      {
        return new CuDSS(*this);
      }

    private:
      /// @brief Creates the cuDSS handle and configuration on first use.
      void initialize()
      {
        auto& res = m_resources;
        if (!res.handle)
          check(cudssCreate(&res.handle), "create the cuDSS handle");
        if (!res.config)
          check(cudssConfigCreate(&res.config), "create the cuDSS configuration");
      }

      /// @brief Allocates the vectors on the device.
      void allocate(Index rows)
      {
        auto& res = m_resources;
        check(res.b.resize(rows), "allocate the right-hand side on the device");
        check(res.x.resize(rows), "allocate the solution on the device");
      }

      template <class T>
      void upload(T* device, const T* host, Index count)
      {
        check(cudaMemcpy(device, host, count * sizeof(T), cudaMemcpyHostToDevice),
          "copy to the device");
      }

      /// @brief Creates the descriptors over the device arrays.
      bool describe(Index rows, Index nonzeros, const Math::SparseIndex* starts,
        const Math::SparseIndex* indices, const ScalarType* values)
      {
        auto& res = m_resources;
        res.destroyDescriptors();
        cudssMatrixType_t type = CUDSS_MTYPE_GENERAL;
        if (m_symmetry == Symmetry::General)
          type = CUDSS_MTYPE_SYMMETRIC;
        else if (m_symmetry == Symmetry::PositiveDefinite)
          type = CUDSS_MTYPE_SPD;
        static_assert(sizeof(Math::SparseIndex) == 4);
        // cuDSS reads the arrays without writing them.
        cudssStatus_t status = cudssMatrixCreateCsr(&res.matrix, rows, rows, nonzeros,
          const_cast<Math::SparseIndex*>(starts), nullptr,
          const_cast<Math::SparseIndex*>(indices), const_cast<ScalarType*>(values),
          CUDA_R_32I, CUDA_R_32I, CUDA_R_64F, type, CUDSS_MVIEW_FULL, CUDSS_BASE_ZERO);
        if (status == CUDSS_STATUS_SUCCESS)
        {
          status = cudssMatrixCreateDn(
            &res.rhs, rows, 1, rows, res.b.data, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR);
        }
        if (status == CUDSS_STATUS_SUCCESS)
        {
          status = cudssMatrixCreateDn(
            &res.solution, rows, 1, rows, res.x.data, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR);
        }
        if (status == CUDSS_STATUS_SUCCESS && !res.data)
          status = cudssDataCreate(res.handle, &res.data);
        if (status != CUDSS_STATUS_SUCCESS)
        {
          record(static_cast<Integer>(status), Optional<Factorization>{});
          return false;
        }
        return true;
      }

      void configure()
      {
        auto& res = m_resources;
        // CSC arrays describe A^T, so an unsymmetric system is solved with
        // its transpose. A symmetric matrix is its own transpose.
        const cudssSolveMode_t mode =
          Arrays::Transposed && m_symmetry == Symmetry::Unsymmetric
          ? CUDSS_SOLVE_MODE_TRANSPOSE
          : CUDSS_SOLVE_MODE_NORMAL;
        check(cudssConfigSet(res.config, CUDSS_CONFIG_SOLVE_MODE, &mode, sizeof(mode)),
          "set the solve mode");
        const int steps = static_cast<int>(m_refinementSteps);
        check(cudssConfigSet(res.config, CUDSS_CONFIG_IR_N_STEPS, &steps, sizeof(steps)),
          "set the refinement steps");
      }

      /**
       * @brief Runs a cuDSS phase and records its outcome.
       * @param phase cuDSS phase.
       * @param survived Factorization stage retained if the phase fails.
       * @returns Whether the phase succeeded.
       */
      bool execute(int phase, Optional<Factorization> survived)
      {
        auto& res = m_resources;
        cudssStatus_t status = cudssExecute(
          res.handle, phase, res.config, res.data, res.matrix, res.solution, res.rhs);
        if (status == CUDSS_STATUS_SUCCESS)
          check(cudaStreamSynchronize(nullptr), "wait for the device");
        int info = 0;
        size_t written = 0;
        if (status == CUDSS_STATUS_SUCCESS)
        {
          status = cudssDataGet(
            res.handle, res.data, CUDSS_DATA_INFO, &info, sizeof(info), &written);
        }
        if (status != CUDSS_STATUS_SUCCESS)
        {
          record(static_cast<Integer>(status), survived);
          return false;
        }
        if (info != 0)
        {
          record(static_cast<Integer>(info), survived);
          return false;
        }
        return true;
      }

      void check(cudaError_t status, StringView operation) const
      {
        if (status != cudaSuccess)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "CUDA failed to " << operation << ": " << cudaGetErrorString(status) << "."
            << Alert::Raise;
        }
      }

      void check(cudssStatus_t status, StringView operation) const
      {
        if (status != CUDSS_STATUS_SUCCESS)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "cuDSS failed to " << operation << " with status "
            << static_cast<Integer>(status) << "." << Alert::Raise;
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
      size_t m_refinementSteps = 0;
      Resources m_resources;
  };
}

#endif
#endif
