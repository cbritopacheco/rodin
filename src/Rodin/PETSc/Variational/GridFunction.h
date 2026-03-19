/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_GRIDFUNCTION_H

/**
 * @file GridFunction.h
 * @brief PETSc specialization of finite element grid functions.
 *
 * This file provides the partial specialization of
 * @ref Rodin::Variational::GridFunction for the PETSc vector backend
 * (`::Vec`). It enables the Rodin finite element framework to store
 * degrees of freedom inside PETSc distributed vectors, supporting both
 * purely local (sequential) and MPI-parallel computations.
 *
 * ## Mathematical Foundation
 *
 * A grid function represents a discrete finite element approximation
 * @f$ u_h \in V_h @f$ defined by its expansion in the finite element
 * basis:
 * @f[
 *   u_h(x) = \sum_{i=1}^{N} u_i \,\phi_i(x)
 * @f]
 * where:
 * - @f$ u_i @f$ are the degrees of freedom (DOF coefficients)
 * - @f$ \phi_i @f$ are the finite element basis functions
 * - @f$ N @f$ is the total number of DOFs in the finite element space
 *
 * In the PETSc specialization the coefficient vector @f$ \mathbf{u} =
 * (u_1, \ldots, u_N) @f$ is stored in a PETSc `Vec` object.  In MPI
 * mode this vector is partitioned across MPI ranks and augmented with a
 * ghost layer so that every rank can read DOFs shared with adjacent mesh
 * partitions without explicit point-to-point communication.
 *
 * ## Ghost-Vector Management
 *
 * In MPI mode the PETSc vector is created with `VecMPISetGhost()`.
 * The ghost indices are computed from the shard's local DOFs that map to
 * global indices outside the local ownership range @f$ [\text{begin},
 * \text{end}) @f$.  Two synchronization patterns are used:
 *
 * | Direction          | PETSc call                                 | Purpose |
 * |--------------------|--------------------------------------------|---------|
 * | owner → ghost      | `VecGhostUpdateBegin/End(INSERT, FORWARD)` | Read remote DOFs before element-level evaluation |
 * | ghost → owner      | `VecGhostUpdateBegin/End(INSERT, REVERSE)` | Scatter locally-modified ghost values back to owners |
 *
 * The `acquire()` / `flush()` / `release()` methods encapsulate
 * these patterns.
 *
 * ## Key Features
 *
 * - **Dual context support**: Compiles for both `Context::Local`
 *   (sequential `PETSC_COMM_SELF`) and `Context::MPI` (distributed) meshes.
 * - **Ghost DOF mapping**: A bidirectional map (`GhostBimap`) translates
 *   between global ghost DOF indices and their local offsets in the
 *   ghosted vector.
 * - **Lazy array access**: The `operator[]` methods lazily `acquire()` the
 *   PETSc array on first use, and the `flush()` method restores it and
 *   triggers the appropriate ghost update.
 * - **Arithmetic operations**: Overloaded `+=`, `-=`, `*=`, `/=` for both
 *   scalar and grid-function operands, delegating to optimized PETSc
 *   routines (`VecShift`, `VecScale`, `VecAXPY`, `VecPointwiseMult`, …).
 * - **Point evaluation**: `operator()` and `getValue()` evaluate
 *   @f$ u_h @f$ at arbitrary geometric points using the finite element
 *   interpolation machinery, with MPI-aware shard lookup.
 * - **Sub-vector import**: `setData()` copies a contiguous slice from
 *   another PETSc vector into this grid function, useful for extracting
 *   sub-solutions from block linear systems.
 * - **Projection**: `project()` assigns DOF values by evaluating a
 *   user-supplied function over a filtered mesh region.
 *
 * ## Usage Example
 *
 * ```cpp
 * // Sequential usage
 * Mesh mesh;
 * mesh = mesh.UniformGrid(Polytope::Type::Triangle, {16, 16});
 * P1 Vh(mesh);
 * PETSc::Variational::GridFunction u(Vh); // zero-initialised
 *
 * // Assign from a lambda
 * u = [](const Point& p) { return std::sin(p.x()) * std::cos(p.y()); };
 *
 * // Arithmetic
 * u *= 2.0;
 *
 * // Evaluate at a point
 * Geometry::Point pt = ...;
 * PetscScalar val = u(pt);
 * ```
 *
 * @see Rodin::Variational::GridFunctionBase,
 *      Rodin::PETSc::Variational::TrialFunction,
 *      Rodin::PETSc::Variational::TestFunction,
 *      Rodin::PETSc::Math::LinearSystem
 */

#include <petsc.h>
#include <petscmacros.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <petscvec.h>
#include <utility>

#include "Rodin/Types.h"
#include "Rodin/Context/Local.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/GridFunction.h"

#include "Rodin/PETSc/Math/Vector.h"

#ifdef RODIN_USE_MPI
#include "Rodin/MPI/Context.h"
#endif

namespace Rodin::Variational
{
  /**
   * @ingroup GridFunctionSpecializations
   * @brief Grid function specialization storing DOF coefficients in a PETSc
   *        vector (`::Vec`).
   *
   * This partial specialization replaces the default Eigen-based coefficient
   * storage with a PETSc `Vec`, making it possible to assemble and solve
   * variational problems entirely within the PETSc ecosystem.  The class
   * supports both sequential (`Context::Local`) and distributed
   * (`Context::MPI`) finite element spaces.
   *
   * ## Ownership and Ghost Model (MPI)
   *
   * When the mesh context is `Context::MPI`, the constructor partitions the
   * global DOF vector across MPI ranks.  Each rank *owns* a contiguous
   * range @f$[\texttt{m\_begin}, \texttt{m\_end})@f$ and additionally
   * maintains *ghost* copies of DOFs needed for element-level operations
   * but owned by neighbouring ranks.  The ghost indices are recorded in
   * a `GhostBimap` so that `operator[]` can transparently translate a
   * global DOF index to the correct local position in the ghosted vector.
   *
   * ## Array Access Protocol
   *
   * PETSc requires explicit lock/unlock of the underlying raw array via
   * `VecGetArray*` / `VecRestoreArray*`.  This class manages the lock
   * state through two helper structs (`ArrayRead` and `ArrayWrite`) and
   * exposes three public methods:
   *
   * | Method      | Direction       | PETSc calls |
   * |-------------|-----------------|-------------|
   * | `acquire()` | begin access    | `VecGhostUpdate(FORWARD)`, `VecGetLocalForm`, `VecGetArrayWrite` |
   * | `flush()`   | end access      | `VecRestoreArrayWrite`, `VecGhostRestoreLocalForm`, `VecGhostUpdate(REVERSE)` |
   * | `release()` | destroy handles | Restores any acquired arrays, then `VecDestroy` |
   *
   * The const overloads use `VecGetArrayRead` / `VecRestoreArrayRead`.
   *
   * @tparam FES Finite element space type (e.g. `P1<Real, Mesh<Context::Local>>`
   *             or `P1<Real, Mesh<Context::MPI>>`).
   *
   * @see Rodin::Variational::GridFunctionBase,
   *      Rodin::PETSc::Variational::GridFunction (convenience wrapper)
   */
  template <class FES>
  class GridFunction<FES, ::Vec>
    : public GridFunctionBase<GridFunction<FES, ::Vec>>
  {
    public:
      /**
       * @brief Bidirectional mapping between global ghost DOF indices and
       *        local ghost offsets within the PETSc ghosted vector.
       *
       * In MPI mode, the finite element space may require DOF values that
       * are owned by remote ranks.  These *ghost* DOFs are appended after
       * the locally-owned entries in the PETSc ghosted vector.  The
       * `GhostBimap` maintains both directions of the mapping:
       *
       * - `left[offset]  → global` : given a local ghost offset, retrieve
       *   the corresponding global DOF index.
       * - `right[global] → offset` : given a global DOF index, find its
       *   position among the local ghost entries.
       *
       * This is used by `operator[]` to convert an arbitrary global DOF
       * index to the correct position in the local (ghosted) array.
       */
      struct GhostBimap
      {
        std::vector<PetscInt> left;           ///< Ghost offset → global DOF index.
        FlatMap<PetscInt, PetscInt> right;     ///< Global DOF index → ghost offset.
      };

      /**
       * @brief State tracking for mutable (write) array access to the PETSc
       *        vector.
       *
       * Stores the raw pointer returned by `VecGetArrayWrite()` as well as
       * the ghosted local form obtained from `VecGhostGetLocalForm()` in
       * MPI mode.  The `acquired` flag prevents redundant lock calls.
       */
      struct ArrayWrite
      {
        Boolean acquired = false;            ///< `true` while the array is locked for writing.
        ::Vec ghost = PETSC_NULLPTR;         ///< Local form of the ghosted vector (MPI only).
        PetscScalar* raw = PETSC_NULLPTR;    ///< Raw pointer to the mutable DOF array.
      };

      /**
       * @brief State tracking for read-only array access to the PETSc
       *        vector.
       *
       * Analogous to @ref ArrayWrite, but uses `VecGetArrayRead()` so that
       * multiple const accessors can coexist safely.
       */
      struct ArrayRead
      {
        Boolean acquired = false;              ///< `true` while the array is locked for reading.
        ::Vec ghost = PETSC_NULLPTR;           ///< Local form of the ghosted vector (MPI only).
        const PetscScalar* raw = PETSC_NULLPTR; ///< Raw pointer to the read-only DOF array.
      };

      /// @brief Finite element space type, e.g. @ref Rodin::Variational::P1.
      using FESType =
        FES;

      /// @brief Underlying PETSc vector data type (`::Vec`), used to store the
      ///        DOF coefficient vector @f$ \mathbf{u} @f$.
      using DataType =
        ::Vec;

      /// @brief Scalar type of each DOF coefficient (`PetscScalar`); determines
      ///        the arithmetic precision of all vector operations.
      using ScalarType =
        PetscScalar;

      /// @brief Range type of the finite element space (e.g. `PetscScalar`
      ///        for scalar FE spaces, or `Math::Vector<PetscScalar>` for
      ///        vector-valued spaces).
      using RangeType =
        typename FormLanguage::Traits<FESType>::RangeType;

      /// @brief Mesh type associated with the finite element space (e.g.
      ///        `Geometry::Mesh<Context::Local>` or `Geometry::Mesh<Context::MPI>`).
      using FESMeshType =
        typename FormLanguage::Traits<FESType>::MeshType;

      /// @brief Context type of the finite element space mesh; either
      ///        @ref Rodin::Context::Local for sequential problems or
      ///        @ref Rodin::Context::MPI for distributed problems.
      using FESMeshContextType =
        typename FormLanguage::Traits<FESMeshType>::ContextType;

      /// @brief Parent CRTP base class providing the generic
      ///        GridFunctionBase interface (projection, interpolation, I/O).
      using Parent =
        GridFunctionBase<GridFunction<FESType, DataType>>;

      using Parent::operator=;

      using Parent::project;

      using Parent::min;

      using Parent::max;

      static_assert(
          std::is_same_v<FESMeshContextType, Context::Local> ||
          std::is_same_v<FESMeshContextType, Context::MPI>);

      /**
       * @brief Constructs a zero-initialised grid function on the given
       *        finite element space.
       *
       * Allocates a PETSc vector whose size matches the number of DOFs in
       * @p fes and sets all entries to zero.
       *
       * ### Local mode (`Context::Local`)
       * Creates a sequential vector on `PETSC_COMM_SELF` with
       * `VecCreate` → `VecSetSizes` → `VecSetFromOptions` → `VecSetUp` →
       * `VecZeroEntries`.
       *
       * ### MPI mode (`Context::MPI`)
       * 1. Queries the ownership range from the finite element space via
       *    `fes.getOwnershipRange(m_begin, m_end)`.
       * 2. Iterates over local shard DOFs to identify ghost indices —
       *    those whose global index falls outside @f$[\text{m\_begin},
       *    \text{m\_end})@f$ — and populates the `GhostBimap`.
       * 3. Creates a distributed vector with `VecCreate` → `VecSetSizes`
       *    → `VecSetFromOptions` → `VecMPISetGhost` → `VecZeroEntries`.
       *
       * @param[in] fes The finite element space that defines the DOF layout,
       *                mesh association, and (in MPI mode) the parallel
       *                partitioning.
       */
      GridFunction(const FESType& fes)
        : Parent(fes),
          m_read{.acquired = false, .raw = PETSC_NULLPTR},
          m_write{.acquired = false, .raw = PETSC_NULLPTR}
      {
        PetscErrorCode ierr;
        auto& data = this->getData();
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          ierr = VecCreate(PETSC_COMM_SELF, &data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetSizes(data, fes.getSize(), PETSC_DECIDE);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetFromOptions(data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetUp(data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecZeroEntries(data);
          assert(ierr == PETSC_SUCCESS);
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          const auto& mesh = fes.getMesh();
          const auto& ctx = mesh.getContext();
          const auto& comm = ctx.getCommunicator();
          const size_t globalSize = fes.getSize();
          const auto& shard = fes.getShard();
          fes.getOwnershipRange(m_begin, m_end);
          const size_t ownedSize = m_end - m_begin;

          Index offset = 0;
          for (size_t i = 0; i < shard.getSize(); ++i)
          {
            const Index global = fes.getGlobalIndex(i);
            if (global < m_begin || global >= m_end)
            {
              auto [it, inserted] = m_ghosts.right.emplace(global, offset);
              assert(inserted);
              offset++;
            }
          }

          m_ghosts.left.resize(m_ghosts.right.size());
          for (const auto& [global, offset] : m_ghosts.right)
          {
            assert(offset >= 0);
            assert(offset < m_ghosts.left.size());
            m_ghosts.left[offset] = global;
          }

          ierr = VecCreate(comm, &data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetSizes(data, ownedSize, globalSize);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecSetFromOptions(data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecMPISetGhost(data, m_ghosts.left.size(), m_ghosts.left.data());
          assert(ierr == PETSC_SUCCESS);

          ierr = VecZeroEntries(data);
          assert(ierr == PETSC_SUCCESS);
        }
        else
        {
          assert(false);
        }
        (void) ierr;
      }

      /**
       * @brief Deep-copy constructor.
       *
       * Duplicates the PETSc vector from @p other with `VecDuplicate` +
       * `VecCopy` so that this instance owns an independent copy of the
       * DOF data.  Ghost mappings and ownership ranges are also copied.
       * Any acquired array state in @p other is **not** propagated; the
       * new object starts with no arrays acquired.
       *
       * @param[in] other Grid function to copy from.
       */
      GridFunction(const GridFunction& other)
        : Parent(other.getFiniteElementSpace()),
          m_begin(other.m_begin),
          m_end(other.m_end),
          m_ghosts(other.m_ghosts),
          m_read{.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR},
          m_write{.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR}
      {
        PetscErrorCode ierr = VecDuplicate(other.m_data, &m_data);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecCopy(other.m_data, m_data);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
      }

      /**
       * @brief Move constructor.
       *
       * Transfers ownership of the PETSc vector handle from @p other using
       * `std::exchange`, leaving @p other in a valid but empty state (null
       * vector, no array access).  No PETSc API calls are made.
       *
       * @param[in] other Grid function to move from.
       */
      GridFunction(GridFunction&& other) noexcept
        : Parent(std::move(other)),
          m_data(std::exchange(other.m_data, PETSC_NULLPTR)),
          m_begin(std::exchange(other.m_begin, 0)),
          m_end(std::exchange(other.m_end, 0)),
          m_ghosts(std::move(other.m_ghosts)),
          m_read{.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR},
          m_write{.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR}
      {
        other.m_read = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};
        other.m_write = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};
      }

      /**
       * @brief Copy assignment operator.
       *
       * Releases any existing PETSc handles via `release()`, then deep-copies
       * the vector from @p other with `VecDuplicate` + `VecCopy`.  Ghost
       * mappings, ownership ranges, and assembly state are also copied.
       *
       * @pre Both grid functions must reference the same finite element space.
       *
       * @param[in] other Grid function to copy from.
       * @returns Reference to `*this`.
       */
      GridFunction& operator=(const GridFunction& other)
      {
        if (this == &other)
          return *this;

        assert(&this->getFiniteElementSpace() == &other.getFiniteElementSpace());

        this->release();

        m_begin = other.m_begin;
        m_end = other.m_end;
        m_ghosts = other.m_ghosts;

        m_read = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};
        m_write = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};

        PetscErrorCode ierr = VecDuplicate(other.m_data, &m_data);
        assert(ierr == PETSC_SUCCESS);

        ierr = VecCopy(other.m_data, m_data);
        assert(ierr == PETSC_SUCCESS);

        (void) ierr;
        return *this;
      }

      /**
       * @brief Move assignment operator.
       *
       * Releases any existing PETSc handles via `release()`, then transfers
       * ownership from @p other.
       *
       * @pre Both grid functions must reference the same finite element space.
       *
       * @param[in] other Grid function to move from.
       * @returns Reference to `*this`.
       */
      GridFunction& operator=(GridFunction&& other) noexcept
      {
        if (this == &other)
          return *this;

        assert(&this->getFiniteElementSpace() == &other.getFiniteElementSpace());

        this->release();

        m_data = std::exchange(other.m_data, PETSC_NULLPTR);
        m_begin = std::exchange(other.m_begin, 0);
        m_end = std::exchange(other.m_end, 0);
        m_ghosts = std::move(other.m_ghosts);

        m_read = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};
        m_write = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};

        other.m_read = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};
        other.m_write = {.acquired = false, .ghost = PETSC_NULLPTR, .raw = PETSC_NULLPTR};

        return *this;
      }

      /// @brief Destructor; calls `release()` to restore any acquired arrays
      ///        and destroy the underlying PETSc vector.
      virtual ~GridFunction()
      {
        this->release();
      }

      /**
       * @brief Finds the minimum DOF value in the grid function.
       *
       * Flushes any pending write access, then delegates to `VecMin`.
       * In MPI mode the result is the global minimum across all ranks.
       *
       * @param[out] idx On return, the global index of the minimum value.
       * @returns The minimum scalar value @f$ \min_i u_i @f$.
       */
      constexpr
      ScalarType min(Index& idx) const
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ScalarType res;
        ierr = VecMin(data, idx, &res);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return res;
      }

      /**
       * @brief Finds the maximum DOF value in the grid function.
       *
       * Flushes any pending write access, then delegates to `VecMax`.
       * In MPI mode the result is the global maximum across all ranks.
       *
       * @param[out] idx On return, the global index of the maximum value.
       * @returns The maximum scalar value @f$ \max_i u_i @f$.
       */
      constexpr
      ScalarType max(Index& idx) const
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ScalarType res;
        ierr = VecMax(data, idx, &res);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return res;
      }

      /**
       * @brief Provides mutable element access by global DOF index.
       *
       * Lazily acquires write access to the underlying PETSc array on first
       * call.  In `Context::Local` mode the index maps directly to the
       * sequential array.  In `Context::MPI` mode, owned indices are
       * translated to @f$ \texttt{global} - \texttt{m\_begin} @f$ while
       * ghost indices are looked up in the `GhostBimap`.
       *
       * @note The caller must call `flush()` after all modifications to
       *       propagate ghost values back to their owners.
       *
       * @param[in] global Global DOF index.
       * @returns Mutable reference to the DOF coefficient @f$ u_i @f$.
       */
      ScalarType& operator[](Index global)
      {
        this->acquire();

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          return m_write.raw[global];
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          PetscInt local;
          if (m_begin <= global && global < m_end)
            local = global - m_begin;
          else
            local = m_end - m_begin + m_ghosts.right.at(global);
          return m_write.raw[local];
        }
        else
        {
          assert(false);
        }
      }

      /**
       * @brief Provides read-only element access by global DOF index.
       *
       * Lazily acquires read-only access to the underlying PETSc array.
       * The index translation logic is identical to the mutable overload.
       *
       * @param[in] global Global DOF index.
       * @returns Const reference to the DOF coefficient @f$ u_i @f$.
       */
      const ScalarType& operator[](Index global) const
      {
        this->acquire();

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          return m_read.raw[global];
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          PetscInt local;
          if (m_begin <= global && global < m_end)
            local = global - m_begin;
          else
            local = m_end - m_begin + m_ghosts.right.at(global);
          return m_read.raw[local];
        }
        else
        {
          assert(false);
        }
      }

      /**
       * @brief Adds a scalar to every DOF value: @f$ u_i \leftarrow u_i + c @f$.
       *
       * Delegates to `VecShift(data, rhs)`.  Only valid when the range type
       * is a scalar.
       *
       * @param[in] rhs Scalar value @f$ c @f$ to add.
       * @returns Reference to `*this`.
       */
      GridFunction& operator+=(const ScalarType& rhs)
      {
        static_assert(std::is_same_v<RangeType, ScalarType>);
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecShift(data, rhs);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      /**
       * @brief Subtracts a scalar from every DOF value: @f$ u_i \leftarrow u_i - c @f$.
       *
       * Delegates to `VecShift(data, -rhs)`.
       *
       * @param[in] rhs Scalar value @f$ c @f$ to subtract.
       * @returns Reference to `*this`.
       */
      GridFunction& operator-=(const ScalarType& rhs)
      {
        static_assert(std::is_same_v<RangeType, ScalarType>);
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecShift(data, -rhs);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      /**
       * @brief Scales every DOF value by a scalar: @f$ u_i \leftarrow c \, u_i @f$.
       *
       * Delegates to `VecScale(data, rhs)`.
       *
       * @param[in] rhs Scalar multiplier @f$ c @f$.
       * @returns Reference to `*this`.
       */
      GridFunction& operator*=(const ScalarType& rhs)
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecScale(data, rhs);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      /**
       * @brief Divides every DOF value by a scalar: @f$ u_i \leftarrow u_i / c @f$.
       *
       * Delegates to `VecScale(data, 1/rhs)`.
       *
       * @param[in] rhs Scalar divisor @f$ c @f$  (must be non-zero).
       * @returns Reference to `*this`.
       */
      GridFunction& operator/=(const ScalarType& rhs)
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecScale(data, 1.0 / rhs);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return static_cast<GridFunction&>(*this);
      }

      /**
       * @brief Component-wise addition of another grid function:
       *        @f$ u_i \leftarrow u_i + v_i @f$.
       *
       * Both grid functions must live on the same finite element space.
       * Delegates to `VecAXPY(data, 1.0, rhs.getData())`.
       *
       * @pre `&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace()`
       *
       * @param[in] rhs Grid function @f$ v @f$ to add.
       * @returns Reference to `*this`.
       */
      GridFunction& operator+=(const GridFunction& rhs)
      {
        assert(&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace());
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecAXPY(data, 1.0, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      /**
       * @brief Component-wise subtraction of another grid function:
       *        @f$ u_i \leftarrow u_i - v_i @f$.
       *
       * Delegates to `VecAXPY(data, -1.0, rhs.getData())`.
       *
       * @pre `&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace()`
       *
       * @param[in] rhs Grid function @f$ v @f$ to subtract.
       * @returns Reference to `*this`.
       */
      GridFunction& operator-=(const GridFunction& rhs)
      {
        assert(&this->getFiniteElementSpace() == &rhs.getFiniteElementSpace());
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecAXPY(data, -1.0, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      /**
       * @brief Component-wise (Hadamard) multiplication by another grid
       *        function: @f$ u_i \leftarrow u_i \cdot v_i @f$.
       *
       * Delegates to `VecPointwiseMult(data, data, rhs.getData())`.
       *
       * @param[in] rhs Grid function @f$ v @f$ to multiply by.
       * @returns Reference to `*this`.
       */
      GridFunction& operator*=(const GridFunction& rhs)
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecPointwiseMult(data, data, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      /**
       * @brief Component-wise division by another grid function:
       *        @f$ u_i \leftarrow u_i / v_i @f$.
       *
       * Delegates to `VecPointwiseDivide(data, data, rhs.getData())`.
       *
       * @param[in] rhs Grid function @f$ v @f$ to divide by (entries must
       *            be non-zero).
       * @returns Reference to `*this`.
       */
      GridFunction& operator/=(const GridFunction& rhs)
      {
        this->flush();
        PetscErrorCode ierr;
        auto& data = this->getData();
        ierr = VecPointwiseDivide(data, data, rhs.getData());
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
        return *this;
      }

      /**
       * @brief Copies a contiguous slice of DOF values from another PETSc
       *        vector into this grid function.
       *
       * This is typically used to extract a sub-solution from a block
       * system's global solution vector.  In local mode the slice
       * @f$[\text{offset}, \text{offset}+N)@f$ is extracted with
       * `VecGetSubVector` / `VecCopy`.  In MPI mode only the locally-owned
       * portion is copied, and ghost values are refreshed by a
       * `VecGhostUpdate(INSERT, FORWARD)`.
       *
       * @param[in] data   Source PETSc vector from which to read.
       * @param[in] offset Global index offset into @p data (default 0).
       * @returns Reference to `*this`.
       */
      GridFunction& setData(const DataType& data, size_t offset = 0)
      {
        this->flush();

        PetscErrorCode ierr;

        // Global size checks
        PetscInt nDst = 0, nSrc = 0;
        ierr = VecGetSize(m_data, &nDst); assert(ierr == PETSC_SUCCESS);
        ierr = VecGetSize(data,   &nSrc); assert(ierr == PETSC_SUCCESS);

        const PetscInt off = static_cast<PetscInt>(offset);
        assert(off >= 0);
        assert(off + nDst <= nSrc);

        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          // Sequential: simple local slice
          IS is = nullptr;
          ierr = ISCreateStride(PETSC_COMM_SELF, nDst, off, 1, &is);
          assert(ierr == PETSC_SUCCESS);

          Vec sub = nullptr;
          ierr = VecGetSubVector(data, is, &sub);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecCopy(sub, m_data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecRestoreSubVector(data, is, &sub);
          assert(ierr == PETSC_SUCCESS);

          ierr = ISDestroy(&is);
          assert(ierr == PETSC_SUCCESS);
        }
        else
        {
          static_assert(std::is_same_v<FESMeshContextType, Context::MPI>);

          // MPI: copy only the owned range of m_data on this rank
          PetscInt rb = 0, re = 0;
          ierr = VecGetOwnershipRange(m_data, &rb, &re);
          assert(ierr == PETSC_SUCCESS);

          const PetscInt nLocal = re - rb;
          assert(nLocal >= 0);

          // Indices in 'data' we want are shifted by 'off'
          const PetscInt start = off + rb;

          // Important: IS local size must be nLocal, not nDst (global)
          MPI_Comm comm;
          PetscObjectGetComm((PetscObject) m_data, &comm);

          IS is = nullptr;
          ierr = ISCreateStride(comm, nLocal, start, 1, &is);
          assert(ierr == PETSC_SUCCESS);

          Vec sub = nullptr;
          ierr = VecGetSubVector(data, is, &sub);
          assert(ierr == PETSC_SUCCESS);

          // Now sub and m_data have identical local layouts -> VecCopy is valid
          ierr = VecCopy(sub, m_data);
          assert(ierr == PETSC_SUCCESS);

          ierr = VecRestoreSubVector(data, is, &sub);
          assert(ierr == PETSC_SUCCESS);

          ierr = ISDestroy(&is);
          assert(ierr == PETSC_SUCCESS);

          // Refresh ghosts from owned values
          ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_FORWARD);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_FORWARD);
          assert(ierr == PETSC_SUCCESS);
        }

        (void)ierr;
        return *this;
      }

      /**
       * @brief Projects a user-supplied function onto the grid function
       *        over a filtered mesh region.
       *
       * Iterates over the polytopes in @p region, applies the predicate
       * @p pred to each polytope, and for those that pass, sets the
       * corresponding DOF values by evaluating the function @p v.
       *
       * The method calls `acquire()` before the loop and `flush()` after
       * it, so the caller does not need to manage array access manually.
       *
       * @tparam Function Callable type compatible with the finite element
       *                  space's projection interface.
       * @tparam Pred     Unary predicate `bool(const Polytope&)` filtering
       *                  polytopes within @p region.
       * @param[in] region Mesh region to iterate (Cells, Faces, Boundary,
       *                   or Interface).
       * @param[in] v      Function to project.
       * @param[in] pred   Predicate selecting polytopes within @p region.
       * @returns Reference to `*this`.
       */
      template <class Function, class Pred>
      GridFunction& project(const Geometry::Region& region, const Function& v, const Pred& pred)
      {
        const auto& fes = this->getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        Geometry::PolytopeIterator it;
        switch (region)
        {
          case Geometry::Region::Cells:
            it = mesh.getCell();
            break;
          case Geometry::Region::Faces:
            it = mesh.getFace();
            break;
          case Geometry::Region::Boundary:
            it = mesh.getBoundary();
            break;
          case Geometry::Region::Interface:
            it = mesh.getInterface();
            break;
        }

        this->acquire();

        while (it)
        {
          const auto& polytope = *it;
          if (pred(polytope))
            this->project({ polytope.getDimension(), polytope.getIndex() }, v);
          ++it;
        }

        this->flush();

        return *this;
      }

      /**
       * @brief Acquires mutable (write) access to the underlying PETSc array.
       *
       * This is the entry point of the array access protocol.  If the
       * array has not already been acquired (`m_write.acquired == false`),
       * the method:
       *
       * ### Local mode
       * - Calls `VecGetArrayWrite(m_data, &raw)`.
       *
       * ### MPI mode
       * 1. Performs a **forward** ghost update
       *    (`VecGhostUpdateBegin/End(INSERT_VALUES, SCATTER_FORWARD)`)
       *    so that ghost entries reflect the latest owner values.
       * 2. Obtains the local (ghosted) form via `VecGhostGetLocalForm`.
       * 3. Locks the local form with `VecGetArrayWrite`.
       *
       * Subsequent calls while the array is already acquired are no-ops.
       *
       * @returns Reference to `*this`.
       */
      GridFunction& acquire()
      {
        PetscErrorCode ierr;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (!m_write.acquired)
          {
            ierr = VecGetArrayWrite(m_data, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = true;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (!m_write.acquired)
          {
            ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_FORWARD);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_FORWARD);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostGetLocalForm(m_data, &m_write.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGetArrayWrite(m_write.ghost, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = true;
          }
        }
        else
        {
          assert(false);
        }
        (void) ierr;
        return *this;
      }

      /**
       * @brief Acquires read-only access to the underlying PETSc array.
       *
       * Behaves like the mutable `acquire()` but uses `VecGetArrayRead`
       * instead, permitting concurrent read access from multiple threads
       * or call sites.  In MPI mode the forward ghost update is still
       * performed to ensure ghost values are up-to-date.
       *
       * @returns Const reference to `*this`.
       */
      const GridFunction& acquire() const
      {
        PetscErrorCode ierr;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (!m_read.acquired)
          {
            ierr = VecGetArrayRead(m_data, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = true;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (!m_read.acquired)
          {
            ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_FORWARD);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_FORWARD);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostGetLocalForm(m_data, &m_read.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGetArrayRead(m_read.ghost, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = true;
          }
        }
        else
        {
          assert(false);
        }
        (void) ierr;
        return *this;
      }

      /**
       * @brief Releases mutable (write) access and scatters locally-modified
       *        ghost values back to their owners.
       *
       * ### Local mode
       * - Calls `VecRestoreArrayWrite(m_data, &raw)`.
       *
       * ### MPI mode
       * 1. Restores the array with `VecRestoreArrayWrite(ghost, &raw)`.
       * 2. Returns the local form via `VecGhostRestoreLocalForm`.
       * 3. Performs a **reverse** ghost update
       *    (`VecGhostUpdateBegin/End(INSERT_VALUES, SCATTER_REVERSE)`)
       *    to propagate ghost modifications back to owners.
       *
       * After this call, `m_write.acquired` is `false` and subsequent
       * reads of the owned portion will reflect any changes that were
       * made through ghost entries.
       *
       * @returns Reference to `*this`.
       */
      GridFunction& flush()
      {
        PetscErrorCode ierr;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (m_write.acquired)
          {
            ierr = VecRestoreArrayWrite(m_data, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = false;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (m_write.acquired)
          {
            ierr = VecRestoreArrayWrite(m_write.ghost, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostRestoreLocalForm(m_data, &m_write.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            m_write.acquired = false;
          }
        }
        else
        {
          assert(false);
        }
        (void) ierr;
        return *this;
      }

      /**
       * @brief Releases read-only access to the PETSc array.
       *
       * Calls `VecRestoreArrayRead` (and, in MPI mode,
       * `VecGhostRestoreLocalForm`) to unlock the array.  No ghost
       * scatter is performed because the array was accessed read-only.
       *
       * @returns Const reference to `*this`.
       */
      const GridFunction& flush() const
      {
        PetscErrorCode ierr;
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          if (m_read.acquired)
          {
            ierr = VecRestoreArrayRead(m_data, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = false;
          }
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          if (m_read.acquired)
          {
            ierr = VecRestoreArrayRead(m_read.ghost, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostRestoreLocalForm(m_data, &m_read.ghost);
            assert(ierr == PETSC_SUCCESS);

            m_read.acquired = false;
          }
        }
        else
        {
          assert(false);
        }
        (void) ierr;
        return *this;
      }

      /// @brief Returns a mutable reference to the raw PETSc `Vec` handle,
      ///        e.g. for passing to PETSc API functions directly.
      auto& getData()
      {
        return m_data;
      }

      /// @brief Returns a read-only reference to the raw PETSc `Vec` handle.
      const DataType& getData() const
      {
        return m_data;
      }

      /// @brief Returns the current read-access state, including the
      ///        acquired flag and the raw pointer (may be null if not acquired).
      const ArrayRead& getArrayRead() const
      {
        return m_read;
      }

      /// @brief Returns the current write-access state, including the
      ///        acquired flag and the raw pointer (may be null if not acquired).
      const ArrayWrite& getArrayWrite() const
      {
        return m_write;
      }

      /// @brief Returns the polynomial order of the finite element space on
      ///        the given polytope, if known.
      ///
      /// For general PETSc grid functions the order is not stored explicitly
      /// and `std::nullopt` is returned.
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const
      {
        return std::nullopt;
      }

      /**
       * @brief Evaluates the grid function at a geometric point.
       *
       * Convenience operator that delegates to `getValue(p)`.
       *
       * @param[in] p The evaluation point @f$ x @f$, carrying the
       *              polytope it belongs to.
       * @returns The interpolated value @f$ u_h(x) @f$.
       */
      decltype(auto) operator()(const Geometry::Point& p) const
      {
        return this->getValue(p);
      }

      /**
       * @brief Evaluates the grid function @f$ u_h @f$ at a geometric point.
       *
       * The evaluation strategy depends on the mesh context:
       *
       * ### Local mode (`Context::Local`)
       * Delegates directly to the parent `GridFunctionBase::getValue()`, which
       * performs the standard finite element interpolation
       * @f$ u_h(x) = \sum_i u_i \phi_i(x) @f$.
       *
       * ### MPI mode (`Context::MPI`)
       * The point may lie on the local shard, on its parent mesh, or on a
       * submesh.  The method checks, in order:
       * 1. Whether the point's polytope belongs to the shard (direct match).
       * 2. Whether the shard geometrically includes the point
       *    (`shard.inclusion(p)`).
       * 3. If the shard is a submesh, whether the point can be restricted to
       *    it (`submesh.restriction(p)`).
       *
       * If none of these succeed, an exception is raised.
       *
       * @param[in] p The evaluation point @f$ x @f$.
       * @returns The interpolated value @f$ u_h(x) @f$.
       *
       * @throws Alert::MemberFunctionException if the point does not belong
       *         to the local shard in MPI mode, or if the mesh context type
       *         is unsupported.
       */
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
        {
          return Parent::getValue(p);
        }
        else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
        {
          static thread_local RangeType s_out;

          const auto& fes = this->getFiniteElementSpace();
          const auto& mesh = fes.getMesh();
          const auto& shard = mesh.getShard();
          const auto& polytope = p.getPolytope();
          const auto& polytopeMesh = polytope.getMesh();

          if (polytopeMesh == shard)
          {
            this->interpolate(s_out, p);
            return s_out;
          }

          if (const auto inclusion = shard.inclusion(p))
          {
            this->interpolate(s_out, *inclusion);
            return s_out;
          }

          if (shard.isSubMesh())
          {
            const auto& submesh = shard.asSubMesh();
            if (const auto restriction = submesh.restriction(p))
            {
              this->interpolate(s_out, *restriction);
              return s_out;
            }
          }

          Alert::MemberFunctionException(*this, __func__)
            << "Point does not belong to the PETSc GridFunction shard."
            << Alert::Raise;

          assert(false);
          return s_out;
        }
        else
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Unsupported mesh context type for PETSc GridFunction."
            << Alert::Raise;
          return Parent::getValue(p);
        }
      }

      /**
       * @brief Releases all PETSc array handles and destroys the vector.
       *
       * Performs a full teardown of the PETSc resources owned by this grid
       * function:
       * 1. If a read array is acquired, calls `VecRestoreArrayRead` (and
       *    `VecGhostRestoreLocalForm` in MPI mode).
       * 2. If a write array is acquired, calls `VecRestoreArrayWrite`,
       *    `VecGhostRestoreLocalForm`, and a reverse ghost scatter.
       * 3. Destroys the PETSc vector with `VecDestroy`.
       *
       * After this call the grid function is in an empty (null-vector)
       * state.  Called automatically by the destructor.  Can also be called
       * explicitly to free resources before the object goes out of scope.
       */
      void release()
      {
        PetscErrorCode ierr;

        if (m_read.acquired)
        {
          assert(m_read.raw);
          if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
          {
            ierr = VecRestoreArrayRead(m_data, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);
          }
          else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
          {
            ierr = VecRestoreArrayRead(m_read.ghost, &m_read.raw);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostRestoreLocalForm(m_data, &m_read.ghost);
            assert(ierr == PETSC_SUCCESS);
          }

          m_read.acquired = false;
          m_read.raw = PETSC_NULLPTR;
          m_read.ghost = PETSC_NULLPTR;
        }

        if (m_write.acquired)
        {
          assert(m_write.raw);
          if constexpr (std::is_same_v<FESMeshContextType, Context::Local>)
          {
            ierr = VecRestoreArrayWrite(m_data, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);
          }
          else if constexpr (std::is_same_v<FESMeshContextType, Context::MPI>)
          {
            ierr = VecRestoreArrayWrite(m_write.ghost, &m_write.raw);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostRestoreLocalForm(m_data, &m_write.ghost);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateBegin(m_data, INSERT_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);

            ierr = VecGhostUpdateEnd(m_data, INSERT_VALUES, SCATTER_REVERSE);
            assert(ierr == PETSC_SUCCESS);
          }

          m_write.acquired = false;
          m_write.raw = PETSC_NULLPTR;
          m_write.ghost = PETSC_NULLPTR;
        }

        if (m_data)
        {
          ierr = VecDestroy(&m_data);
          assert(ierr == PETSC_SUCCESS);
          m_data = PETSC_NULLPTR;
        }

        (void) ierr;
      }

    private:
      DataType m_data;          ///< Underlying PETSc `Vec` storing the DOF coefficient vector @f$ \mathbf{u} @f$.
      size_t m_begin, m_end;    ///< Owned DOF range @f$[\texttt{m\_begin}, \texttt{m\_end})@f$ (MPI mode); unused in Local mode.
      GhostBimap m_ghosts;      ///< Bidirectional ghost DOF index mapping (MPI mode); empty in Local mode.

      mutable ArrayRead m_read;   ///< Read-access state (mutable to allow const `acquire`/`flush`).
      mutable ArrayWrite m_write; ///< Write-access state (mutable to allow const `acquire`/`flush`).
  };
}

namespace Rodin::PETSc::Variational
{
  /**
   * @brief Convenience wrapper that inherits all functionality from
   *        @ref Rodin::Variational::GridFunction<FES, ::Vec>.
   *
   * This thin derived class lives in the `Rodin::PETSc::Variational`
   * namespace and pulls in all base-class constructors, operators, and
   * methods via `using` declarations.  Its sole purpose is to enable
   * cleaner user-facing syntax:
   *
   * ```cpp
   * PETSc::Variational::GridFunction u(Vh);
   * ```
   *
   * instead of the fully-qualified specialization name.
   *
   * @tparam FES Finite element space type (deduced by CTAD).
   *
   * @see Rodin::Variational::GridFunction<FES, ::Vec>
   */
  template <class FES>
  class GridFunction
    : public Rodin::Variational::GridFunction<FES, ::Vec>
  {
    public:
      /// @brief Parent specialization type.
      using Parent = Rodin::Variational::GridFunction<FES, ::Vec>;
      using Parent::Parent;       ///< Inherit all constructors.
      using Parent::operator[];   ///< Inherit element-access operators.
      using Parent::operator();   ///< Inherit point-evaluation operator.
      using Parent::getValue;     ///< Inherit getValue.
      using Parent::operator=;    ///< Inherit assignment operators.
      using Parent::operator+=;   ///< Inherit addition operators.
      using Parent::operator-=;   ///< Inherit subtraction operators.
      using Parent::operator*=;   ///< Inherit multiplication operators.
      using Parent::operator/=;   ///< Inherit division operators.
  };

  /**
   * @ingroup RodinCTAD
   * @brief Deduction guide for PETSc::Variational::GridFunction.
   */
  template <class FES>
  GridFunction(const FES&) -> GridFunction<FES>;
}

namespace Rodin::FormLanguage
{
  /**
   * @brief Form-language traits specialization for PETSc grid functions.
   *
   * Provides the `FESType` and `DataType` aliases that the Rodin form
   * language uses to determine the finite element space and storage
   * backend associated with a given grid function type.
   *
   * @tparam FES Finite element space type.
   */
  template <class FES>
  struct Traits<PETSc::Variational::GridFunction<FES>>
  {
    /// @brief Finite element space type associated with this grid function.
    using FESType = FES;
    /// @brief Data storage type (PETSc vector `::Vec`).
    using DataType = ::Vec;
  };
}

#endif
