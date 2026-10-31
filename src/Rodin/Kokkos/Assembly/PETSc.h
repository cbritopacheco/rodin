/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_ASSEMBLY_PETSC_H
#define RODIN_KOKKOS_ASSEMBLY_PETSC_H

/**
 * @file
 * @brief Device assembly of PETSc problems on resident Kokkos meshes.
 *
 * The integrators write their entries straight into coordinate lists in
 * device memory, which PETSc consumes through MatSetPreallocationCOO and
 * MatSetValuesCOO. The matrix and the vectors are created with the PETSc
 * Kokkos types, so a solve on the device never copies the system to the
 * host.
 */

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

#include <petscmat.h>
#include <petscvec.h>

#include "Rodin/PETSc/Assembly/VectorSetup.h"
#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/PETSc/Variational/Problem.h"

#include "Rodin/Kokkos/Assembly/Kokkos.h"

#ifndef PETSC_HAVE_KOKKOS
#error "Device assembly of PETSc problems requires PETSc configured with Kokkos."
#endif

namespace Rodin::Assembly::Internal
{
  /**
   * @brief Sparsity pattern last handed to a PETSc matrix.
   *
   * PETSc analyses the coordinate list in MatSetPreallocationCOO, which is
   * much more expensive than MatSetValuesCOO. The pattern is kept to skip
   * that analysis when the same matrix is assembled again with the same
   * coordinates.
   *
   * @tparam MemorySpace Kokkos memory space of the coordinates.
   */
  template <class MemorySpace>
  struct COOPattern
  {
    /// @brief PETSc identifier of the matrix, or -1 before the first assembly.
    PetscInt64 matrix = -1;

    /// @brief Row coordinates.
    ::Kokkos::View<PetscInt*, MemorySpace> rows;

    /// @brief Column coordinates.
    ::Kokkos::View<PetscInt*, MemorySpace> cols;
  };

  /**
   * @brief PETSc linear system assembled in device memory, with value
   * constraints eliminated.
   *
   * The elimination is the one of ConstrainedSystem: rows of fixed degrees
   * of freedom become identity rows holding the fixed value, and the
   * columns of fixed degrees of freedom move to the right-hand side.
   * Eliminated entries keep their coordinates with a zero value, and every
   * diagonal entry is listed, so the pattern of the matrix does not depend
   * on the constraints.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class DeviceCOOSystem
  {
    public:
      /// @brief Mesh context.
      using ContextType = ::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>;

      /// @brief Execution space of the kernels.
      using ExecutionSpace = typename ContextType::ExecutionSpaceType;

      /// @brief Memory space of the system.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Device integrator triplets.
      using TripletsType = ::Rodin::Kokkos::Variational::DeviceTriplets<MemorySpace>;

      /// @brief Device integrator entries.
      using EntriesType = ::Rodin::Kokkos::Variational::DeviceEntries<MemorySpace>;

      /// @brief Pattern cache type.
      using PatternType = COOPattern<MemorySpace>;

      static_assert(std::is_same_v<PetscScalar, Real>,
        "Device assembly of PETSc problems requires real PETSc scalars.");

      /**
       * @brief Creates an empty system.
       * @param execution Execution space instance of the mesh.
       * @param rows Number of rows, the size of the test spaces.
       * @param cols Number of columns, the size of the trial spaces.
       */
      DeviceCOOSystem(const ExecutionSpace& execution, size_t rows, size_t cols)
        : m_execution(execution),
          m_rows(rows),
          m_cols(cols),
          m_fixed(std::max(rows, cols), 0),
          m_values(std::max(rows, cols), Real(0))
      {
        if (std::max(rows, cols) >
          static_cast<size_t>(std::numeric_limits<PetscInt>::max()))
        {
          Alert::Exception() << "The system is too large for the PETSc index type."
                             << Alert::Raise;
        }
      }

      /**
       * @brief Assembles a Dirichlet condition and fixes its degrees of
       * freedom.
       * @param dbc Value Dirichlet condition.
       * @param offset Offset of the constrained space in the system.
       */
      void fix(Variational::DirichletBCBase<Real>& dbc, size_t offset)
      {
        using ValueDOFs = typename Variational::DirichletBCBase<Real>::ValueDOFs;
        dbc.assemble();
        const auto* dofs = std::get_if<ValueDOFs>(&dbc.getDOFs());
        if (!dofs)
        {
          Alert::Exception()
            << "Identification conditions are not supported on Kokkos meshes."
            << Alert::Raise;
        }
        for (const auto& [local, value] : *dofs)
        {
          m_fixed[offset + local] = 1;
          m_values[offset + local] = value;
        }
      }

      /// @brief Whether a degree of freedom is fixed.
      bool hasFixed() const
      {
        return std::find(m_fixed.begin(), m_fixed.end(), 1) != m_fixed.end();
      }

      /**
       * @brief Integrates a local integrator on the device.
       * @param lbfi Local integrator, possibly negated.
       * @param row Offset of its test space. @param col Offset of its trial space.
       */
      void addMatrix(const Variational::LocalBilinearFormIntegratorBase<Real>& lbfi,
        Index row, Index col)
      {
        const auto [integrator, sign] = getDeviceIntegrator<Backend>(lbfi);
        m_matrices.push_back({integrator->assemble(), row, col, sign});
      }

      /**
       * @brief Integrates a linear integrator on the device.
       * @param lfi Integrator, possibly negated.
       * @param row Offset of its test space.
       * @param sign Sign of the contributions.
       */
      void addVector(
        const Variational::LinearFormIntegratorBase<Real>& lfi, Index row, Real sign)
      {
        const auto [integrator, s] = getDeviceIntegrator<Backend>(lfi, sign);
        m_vectors.push_back({integrator->assemble(), row, s});
      }

      /**
       * @brief Eliminates the constraints and fills the PETSc objects.
       * @param pattern Pattern last handed to the matrix.
       * @param A Matrix, filled when not null.
       * @param b Right-hand side, filled when not null.
       */
      void finish(PatternType& pattern, ::Mat A, ::Vec b)
      {
        const size_t n = std::max(m_rows, m_cols);
        const size_t diagonal = std::min(m_rows, m_cols);
        const auto fixed = upload(m_fixed, "Rodin fixed DOFs");
        const auto values = upload(m_values, "Rodin fixed values");
        ::Kokkos::View<Real*, MemorySpace> rhs(
          ::Kokkos::view_alloc(m_execution, std::string("Rodin right-hand side")), n);

        size_t count = diagonal;
        for (const auto& m : m_matrices)
          count += m.triplets.values.extent(0);
        ::Kokkos::View<PetscInt*, MemorySpace> rows(
          ::Kokkos::view_alloc(
            m_execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO rows")),
          count);
        ::Kokkos::View<PetscInt*, MemorySpace> cols(
          ::Kokkos::view_alloc(
            m_execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO cols")),
          count);
        ::Kokkos::View<PetscScalar*, MemorySpace> coo(
          ::Kokkos::view_alloc(
            m_execution, ::Kokkos::WithoutInitializing, std::string("Rodin COO values")),
          count);

        size_t offset = 0;
        for (const auto& m : m_matrices)
        {
          scatterMatrix(m_execution, m.triplets, m.row, m.col, m.sign, offset, fixed,
            values, rows, cols, coo, rhs);
          offset += m.triplets.values.extent(0);
        }
        scatterDiagonal(m_execution, offset, diagonal, fixed, rows, cols, coo);
        for (const auto& v : m_vectors)
          scatterVector(m_execution, v.entries, v.row, v.sign, fixed, rhs);
        imposeFixed(m_execution, m_rows, fixed, values, rhs);
        m_execution.fence("Rodin PETSc system");

        if (A)
          setMatrix(pattern, A, rows, cols, coo);
        if (b)
          setVector(b, rhs);
      }

      /**
       * @brief Writes the entries of a local integrator.
       *
       * Public as the enclosing function of device lambdas.
       */
      template <class Fixed, class Values, class Indices, class Scalars, class RHS>
      static void scatterMatrix(const ExecutionSpace& execution,
        const TripletsType& triplets, Index rowOffset, Index colOffset, Real sign,
        size_t offset, const Fixed& fixed, const Values& values, const Indices& rows,
        const Indices& cols, const Scalars& coo, const RHS& rhs)
      {
        const auto tr = triplets.rows;
        const auto tc = triplets.cols;
        const auto tv = triplets.values;
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::matrix",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, tv.extent(0)),
          KOKKOS_LAMBDA(const size_t e) {
            const Index r = rowOffset + tr(e);
            const Index c = colOffset + tc(e);
            PetscScalar v = sign * tv(e);
            if (fixed(r))
            {
              v = 0;
            }
            else if (fixed(c))
            {
              ::Kokkos::atomic_add(&rhs(r), -v * values(c));
              v = 0;
            }
            rows(offset + e) = static_cast<PetscInt>(r);
            cols(offset + e) = static_cast<PetscInt>(c);
            coo(offset + e) = v;
          });
      }

      /**
       * @brief Lists every diagonal entry, one on fixed rows.
       *
       * Public as the enclosing function of device lambdas.
       */
      template <class Fixed, class Indices, class Scalars>
      static void scatterDiagonal(const ExecutionSpace& execution, size_t offset,
        size_t count, const Fixed& fixed, const Indices& rows, const Indices& cols,
        const Scalars& coo)
      {
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::diagonal",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          KOKKOS_LAMBDA(const size_t i) {
            rows(offset + i) = static_cast<PetscInt>(i);
            cols(offset + i) = static_cast<PetscInt>(i);
            coo(offset + i) = fixed(i) ? PetscScalar(1) : PetscScalar(0);
          });
      }

      /**
       * @brief Adds the entries of a linear integrator to the free rows.
       *
       * Public as the enclosing function of device lambdas.
       */
      template <class Fixed, class RHS>
      static void scatterVector(const ExecutionSpace& execution,
        const EntriesType& entries, Index rowOffset, Real sign, const Fixed& fixed,
        const RHS& rhs)
      {
        const auto indices = entries.indices;
        const auto values = entries.values;
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::vector",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, values.extent(0)),
          KOKKOS_LAMBDA(const size_t e) {
            const Index i = rowOffset + indices(e);
            if (!fixed(i))
              ::Kokkos::atomic_add(&rhs(i), sign * values(e));
          });
      }

      /**
       * @brief Writes the fixed values on the fixed rows.
       *
       * Public as the enclosing function of device lambdas.
       */
      template <class Fixed, class Values, class RHS>
      static void imposeFixed(const ExecutionSpace& execution, size_t count,
        const Fixed& fixed, const Values& values, const RHS& rhs)
      {
        ::Kokkos::parallel_for(
          "Rodin::Assembly::Kokkos::PETSc::fixed",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, count),
          KOKKOS_LAMBDA(const size_t i) {
            if (fixed(i))
              rhs(i) = values(i);
          });
      }

      /**
       * @brief Counts the coordinates that differ between two patterns.
       *
       * Public as the enclosing function of device lambdas.
       */
      template <class Indices>
      static size_t countChanges(const ExecutionSpace& execution, const Indices& rows,
        const Indices& cols, const Indices& previousRows, const Indices& previousCols)
      {
        size_t changes = 0;
        ::Kokkos::parallel_reduce(
          "Rodin::Assembly::Kokkos::PETSc::pattern",
          ::Kokkos::RangePolicy<ExecutionSpace>(execution, 0, rows.extent(0)),
          KOKKOS_LAMBDA(const size_t e, size_t& sum) {
            sum += rows(e) != previousRows(e) || cols(e) != previousCols(e);
          },
          changes);
        return changes;
      }

    private:
      struct Matrix
      {
        TripletsType triplets;
        Index row;
        Index col;
        Real sign;
      };

      struct Vector
      {
        EntriesType entries;
        Index row;
        Real sign;
      };

      template <class T>
      ::Kokkos::View<T*, MemorySpace> upload(
        const std::vector<T>& host, const char* label) const
      {
        ::Kokkos::View<T*, MemorySpace> res(
          ::Kokkos::view_alloc(
            m_execution, ::Kokkos::WithoutInitializing, std::string(label)),
          host.size());
        const ::Kokkos::View<const T*, ::Kokkos::HostSpace,
          ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>
          source(host.data(), host.size());
        ::Kokkos::deep_copy(m_execution, res, source);
        return res;
      }

      /**
       * @brief Hands the coordinate list to the matrix.
       *
       * The matrix is laid out as a PETSc Kokkos matrix unless it already has
       * a type. The pattern is analysed again only when the coordinates
       * changed since the last assembly of the same matrix.
       */
      void setMatrix(PatternType& pattern, ::Mat A,
        const ::Kokkos::View<PetscInt*, MemorySpace>& rows,
        const ::Kokkos::View<PetscInt*, MemorySpace>& cols,
        const ::Kokkos::View<PetscScalar*, MemorySpace>& coo) const
      {
        PetscErrorCode ierr;
        PetscInt64 id = -1;
        ierr = PetscObjectGetId(reinterpret_cast<PetscObject>(A), &id);
        assert(ierr == PETSC_SUCCESS);

        const bool reuse = pattern.matrix == id &&
          pattern.rows.extent(0) == rows.extent(0) &&
          countChanges(m_execution, rows, cols, pattern.rows, pattern.cols) == 0;
        if (!reuse)
        {
          PetscBool assembled = PETSC_FALSE;
          ierr = MatAssembled(A, &assembled);
          assert(ierr == PETSC_SUCCESS);
          if (!assembled)
          {
            const PetscInt m = static_cast<PetscInt>(m_rows);
            const PetscInt n = static_cast<PetscInt>(m_cols);
            ierr = MatSetSizes(A, m, n, m, n);
            assert(ierr == PETSC_SUCCESS);
            MatType type = nullptr;
            ierr = MatGetType(A, &type);
            assert(ierr == PETSC_SUCCESS);
            if (!type)
            {
              ierr = MatSetType(A, MATAIJKOKKOS);
              assert(ierr == PETSC_SUCCESS);
            }
          }
          // PETSc may reorder the coordinates it is given, so it gets copies.
          auto hostRows = ::Kokkos::create_mirror(::Kokkos::HostSpace(), rows);
          auto hostCols = ::Kokkos::create_mirror(::Kokkos::HostSpace(), cols);
          ::Kokkos::deep_copy(hostRows, rows);
          ::Kokkos::deep_copy(hostCols, cols);
          ierr = MatSetPreallocationCOO(
            A, static_cast<PetscCount>(rows.extent(0)), hostRows.data(), hostCols.data());
          assert(ierr == PETSC_SUCCESS);
          pattern = {id, rows, cols};
        }

        PetscBool device = PETSC_FALSE;
        ierr = PetscObjectTypeCompareAny(reinterpret_cast<PetscObject>(A), &device,
          MATSEQAIJKOKKOS, MATMPIAIJKOKKOS, MATSEQAIJCUSPARSE, MATMPIAIJCUSPARSE, "");
        assert(ierr == PETSC_SUCCESS);
        if (device ||
          ::Kokkos::SpaceAccessibility<::Kokkos::HostSpace, MemorySpace>::accessible)
        {
          ierr = MatSetValuesCOO(A, coo.data(), INSERT_VALUES);
        }
        else
        {
          const auto host =
            ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), coo);
          ierr = MatSetValuesCOO(A, host.data(), INSERT_VALUES);
        }
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      /**
       * @brief Copies the right-hand side where the vector stores it.
       */
      void setVector(::Vec b, const ::Kokkos::View<Real*, MemorySpace>& rhs) const
      {
        using Unmanaged = ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>;
        const auto source = ::Kokkos::subview(rhs, std::make_pair(size_t(0), m_rows));
        PetscScalar* array = nullptr;
        PetscMemType type;
        PetscErrorCode ierr = VecGetArrayWriteAndMemType(b, &array, &type);
        assert(ierr == PETSC_SUCCESS);
        if (PetscMemTypeHost(type))
        {
          ::Kokkos::deep_copy(
            ::Kokkos::View<PetscScalar*, ::Kokkos::HostSpace, Unmanaged>(array, m_rows),
            source);
        }
        else
        {
          ::Kokkos::deep_copy(
            ::Kokkos::View<PetscScalar*, ::Kokkos::DefaultExecutionSpace::memory_space,
              Unmanaged>(array, m_rows),
            source);
        }
        ierr = VecRestoreArrayWriteAndMemType(b, &array);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      ExecutionSpace m_execution;
      size_t m_rows;
      size_t m_cols;
      std::vector<unsigned char> m_fixed;
      std::vector<Real> m_values;
      std::vector<Matrix> m_matrices;
      std::vector<Vector> m_vectors;
  };

  /**
   * @brief Prepares the PETSc objects of a device assembly.
   * @param axb PETSc linear system.
   * @param rows Number of rows. @param cols Number of columns.
   * @param doVector Whether the right-hand side is assembled.
   */
  inline void prepareDeviceSystem(
    PETSc::Math::LinearSystem& axb, size_t rows, size_t cols, bool doVector)
  {
    PetscErrorCode ierr;
    const PetscInt m = static_cast<PetscInt>(rows);
    const PetscInt n = static_cast<PetscInt>(cols);
    if (doVector)
    {
      ierr = PETSc::Assembly::VectorSetup(axb.getVector())
               .prepare({m, m, VECKOKKOS, false, false});
      assert(ierr == PETSC_SUCCESS);
    }
    // The solution carries the previous iterate as the initial guess.
    ierr = PETSc::Assembly::VectorSetup(axb.getSolution())
             .prepare({n, n, VECKOKKOS, false, false});
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  /**
   * @brief Rejects the PETSc problem terms without a device implementation.
   * @param pb Problem body.
   */
  template <class ProblemBody>
  void requireDevicePETScProblem(const ProblemBody& pb)
  {
    requireDeviceProblem(pb);
    if (pb.getBFs().size() > 0 || pb.getLFs().size() > 0)
    {
      Alert::Exception() << "Preassembled PETSc forms are not supported on Kokkos meshes."
                         << Alert::Raise;
    }
  }
}

namespace Rodin::Assembly
{
  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a single field PETSc problem on a Kokkos mesh.
   *
   * The integrators of the problem, possibly negated, are integrated on the
   * device and handed to PETSc as coordinate lists in device memory. The
   * value Dirichlet conditions are eliminated on the device, as the host
   * assemblies do.
   *
   * Global integrators, preassembled forms, periodic conditions and
   * identification conditions are not supported.
   *
   * @tparam U Trial function on a resident Kokkos mesh.
   * @tparam V Test function on the same mesh.
   */
  template <class U, class V>
  class Kokkos<PETSc::Math::LinearSystem,
    Variational::Problem<PETSc::Math::LinearSystem, U, V>>
    final : public AssemblyBase<PETSc::Math::LinearSystem,
              Variational::Problem<PETSc::Math::LinearSystem, U, V>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType = PETSc::Math::LinearSystem;

      /// @brief Problem type.
      using ProblemType = Variational::Problem<LinearSystemType, U, V>;

      /// @brief Parent class.
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the trial mesh context.
      using BackendType = typename FormLanguage::Traits<typename FormLanguage::Traits<
        typename FormLanguage::Traits<U>::FESType>::MeshType>::ContextType::BackendType;

      /// @brief Device system type.
      using SystemType = Internal::DeviceCOOSystem<BackendType>;

      /// @brief Default constructor.
      Kokkos() = default;

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Kokkos(const Kokkos& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Kokkos(Kokkos&& other)
        : Parent(std::move(other)),
          m_pattern(std::move(other.m_pattern))
      {}

      /**
       * @brief Assigns the assembly policy, forgetting the cached pattern.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other)
      {
        Parent::operator=(other);
        m_pattern = {};
        return *this;
      }

      /**
       * @brief Moves the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Integrates the problem on the device and builds its system.
       * @param axb Assembled linear system.
       * @param input Problem input.
       */
      void execute(LinearSystemType& axb, const InputType& input) override
      {
        assemble(axb, input, true, true);
      }

      /**
       * @brief Assembles one side of the system.
       * @param axb Assembled linear system.
       * @param input Problem input.
       * @param target Side of the system to assemble.
       */
      void execute(
        LinearSystemType& axb, const InputType& input, Variational::AssemblyTarget target)
      {
        const bool lhs = target == Variational::AssemblyTarget::LHS;
        assemble(axb, input, lhs, !lhs);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }

    private:
      void assemble(
        LinearSystemType& axb, const InputType& input, bool doMatrix, bool doVector)
      {
        auto& pb = input.getProblemBody();
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();
        Internal::requireDevicePETScProblem(pb);
        const size_t rows = v.getFiniteElementSpace().getSize();
        const size_t cols = u.getFiniteElementSpace().getSize();
        SystemType system(
          u.getFiniteElementSpace().getMesh().getContext().getExecutionSpace(), rows,
          cols);
        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() == u.getUUID())
            system.fix(dbc, 0);
        }
        // The matrix is still needed on the right-hand side to lift the
        // fixed values.
        if (doMatrix || system.hasFixed())
        {
          for (const auto& lbfi : pb.getLocalBFIs())
            system.addMatrix(lbfi, 0, 0);
        }
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        if (doVector)
        {
          for (const auto& lfi : pb.getLFIs())
            system.addVector(lfi, 0, Real(-1));
        }
        Internal::prepareDeviceSystem(axb, rows, cols, doVector);
        system.finish(m_pattern, doMatrix ? axb.getOperator() : nullptr,
          doVector ? axb.getVector() : nullptr);
      }

      Internal::COOPattern<typename SystemType::MemorySpace> m_pattern;
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a mixed PETSc problem on a Kokkos mesh.
   *
   * Each field occupies a block of the system, at the offsets given by the
   * problem input, and the system is assembled as in the single field case.
   *
   * Global integrators, preassembled forms, periodic conditions and
   * identification conditions are not supported.
   */
  template <class U1, class U2, class U3, class... Us>
  class Kokkos<PETSc::Math::LinearSystem,
    Variational::Problem<PETSc::Math::LinearSystem, U1, U2, U3, Us...>>
    final : public AssemblyBase<PETSc::Math::LinearSystem,
              Variational::Problem<PETSc::Math::LinearSystem, U1, U2, U3, Us...>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType = PETSc::Math::LinearSystem;

      /// @brief Problem type.
      using ProblemType = Variational::Problem<LinearSystemType, U1, U2, U3, Us...>;

      /// @brief Parent class.
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the mesh context of the first field.
      using BackendType = typename FormLanguage::Traits<typename FormLanguage::Traits<
        typename FormLanguage::Traits<U1>::FESType>::MeshType>::ContextType::BackendType;

      /// @brief Device system type.
      using SystemType = Internal::DeviceCOOSystem<BackendType>;

      /// @brief Default constructor.
      Kokkos() = default;

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Kokkos(const Kokkos& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Kokkos(Kokkos&& other)
        : Parent(std::move(other)),
          m_pattern(std::move(other.m_pattern))
      {}

      /**
       * @brief Assigns the assembly policy, forgetting the cached pattern.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other)
      {
        Parent::operator=(other);
        m_pattern = {};
        return *this;
      }

      /**
       * @brief Moves the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Integrates the problem on the device and builds its system.
       * @param axb Assembled linear system.
       * @param input Problem input.
       */
      void execute(LinearSystemType& axb, const InputType& input) override
      {
        assemble(axb, input, true, true);
      }

      /**
       * @brief Assembles one side of the system.
       * @param axb Assembled linear system.
       * @param input Problem input.
       * @param target Side of the system to assemble.
       */
      void execute(
        LinearSystemType& axb, const InputType& input, Variational::AssemblyTarget target)
      {
        const bool lhs = target == Variational::AssemblyTarget::LHS;
        assemble(axb, input, lhs, !lhs);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }

    private:
      void assemble(
        LinearSystemType& axb, const InputType& input, bool doMatrix, bool doVector)
      {
        auto& pb = input.getProblemBody();
        const auto& trialOffsets = input.getTrialOffsets();
        const auto& testOffsets = input.getTestOffsets();
        const auto& trialUUIDMap = input.getTrialUUIDMap();
        const auto& testUUIDMap = input.getTestUUIDMap();
        const auto trialOffset = [&](const auto& uuid) -> Index {
          const auto it = trialUUIDMap.left.find(uuid);
          assert(it != trialUUIDMap.left.end());
          return static_cast<Index>(trialOffsets[it->second]);
        };
        const auto testOffset = [&](const auto& uuid) -> Index {
          const auto it = testUUIDMap.left.find(uuid);
          assert(it != testUUIDMap.left.end());
          return static_cast<Index>(testOffsets[it->second]);
        };

        Internal::requireDevicePETScProblem(pb);
        const size_t rows = input.getTotalTestSize();
        const size_t cols = input.getTotalTrialSize();
        const auto& u1 = input.getTrialFunctions().template get<0>().get();
        SystemType system(
          u1.getFiniteElementSpace().getMesh().getContext().getExecutionSpace(), rows,
          cols);
        for (auto& dbc : pb.getDBCs())
          system.fix(dbc, trialOffset(dbc.getOperand().getUUID()));
        if (doMatrix || system.hasFixed())
        {
          for (const auto& lbfi : pb.getLocalBFIs())
          {
            system.addMatrix(lbfi, testOffset(lbfi.getTestFunction().getUUID()),
              trialOffset(lbfi.getTrialFunction().getUUID()));
          }
        }
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        if (doVector)
        {
          for (const auto& lfi : pb.getLFIs())
            system.addVector(lfi, testOffset(lfi.getTestFunction().getUUID()), Real(-1));
        }
        Internal::prepareDeviceSystem(axb, rows, cols, doVector);
        system.finish(m_pattern, doMatrix ? axb.getOperator() : nullptr,
          doVector ? axb.getVector() : nullptr);
      }

      Internal::COOPattern<typename SystemType::MemorySpace> m_pattern;
  };
}

#endif
