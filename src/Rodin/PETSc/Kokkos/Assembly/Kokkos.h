/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_KOKKOS_ASSEMBLY_KOKKOS_H
#define RODIN_PETSC_KOKKOS_ASSEMBLY_KOKKOS_H

/**
 * @file
 * @brief Kokkos assembly policies targeting PETSc matrices, vectors and systems.
 */

#include "COOAssembly.h"

#ifndef PETSC_HAVE_KOKKOS
#error "Device assembly of PETSc problems requires PETSc configured with Kokkos."
#endif

namespace Rodin::Assembly::Internal
{
  /**
   * @brief Prepares the PETSc objects of a distributed device assembly.
   * @param axb PETSc linear system.
   * @param rows Number of owned rows. @param globalRows Global number of rows.
   * @param cols Number of owned columns. @param globalCols Global number of columns.
   * @param doVector Whether the right-hand side is assembled.
   */
  inline void prepareSystem(PETSc::Math::LinearSystem& axb, size_t rows,
    size_t globalRows, size_t cols, size_t globalCols, bool doVector)
  {
    [[maybe_unused]] PetscErrorCode ierr;
    const PetscInt m = static_cast<PetscInt>(rows);
    const PetscInt n = static_cast<PetscInt>(cols);
    const PetscInt M = static_cast<PetscInt>(globalRows);
    const PetscInt N = static_cast<PetscInt>(globalCols);
    if (doVector)
    {
      ierr = PETSc::Assembly::VectorSetup(axb.getVector())
               .prepare({m, M, VECKOKKOS, false, false});
      assert(ierr == PETSC_SUCCESS);
    }
    // The solution carries the previous iterate as the initial guess.
    ierr = PETSc::Assembly::VectorSetup(axb.getSolution())
             .prepare({n, N, VECKOKKOS, false, false});
    assert(ierr == PETSC_SUCCESS);
  }

  /**
   * @brief Prepares the PETSc objects of a device assembly.
   * @param axb PETSc linear system.
   * @param rows Number of rows. @param cols Number of columns.
   * @param doVector Whether the right-hand side is assembled.
   */
  inline void prepareSystem(
    PETSc::Math::LinearSystem& axb, size_t rows, size_t cols, bool doVector)
  {
    prepareSystem(axb, rows, rows, cols, cols, doVector);
  }
}

namespace Rodin::Assembly
{
  /**
   * @ingroup KokkosAssemblySpecializations
   * @brief Device assembly of a PETSc bilinear form on a Kokkos mesh.
   *
   * The integrators of the form, possibly negated, are integrated on the
   * device and handed to PETSc as a coordinate list in device memory, as in
   * the problem assemblies. The pattern holds the diagonal.
   *
   * @tparam Solution Solution variable type.
   * @tparam TrialFES Trial space on a resident Kokkos mesh.
   * @tparam TestFES Test space on the same resident Kokkos mesh.
   */
  template <class Solution, class TrialFES, class TestFES>
  class Kokkos<::Mat, Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>> final
    : public AssemblyBase<::Mat,
        Variational::BilinearForm<Solution, TrialFES, TestFES, ::Mat>>
  {
    public:
      /// @brief Operator type.
      using OperatorType = ::Mat;

      /// @brief Bilinear form type.
      using BilinearFormType =
        Variational::BilinearForm<Solution, TrialFES, TestFES, OperatorType>;

      /// @brief Parent class.
      using Parent = AssemblyBase<OperatorType, BilinearFormType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the trial mesh context.
      using BackendType = typename FormLanguage::Traits<
        typename FormLanguage::Traits<TrialFES>::MeshType>::ContextType::BackendType;

      /// @brief Device system type.
      using SystemType = Internal::COOAssembly<BackendType>;

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
          m_cache(std::move(other.m_cache)),
          m_coefficients(std::move(other.m_coefficients))
      {}

      /**
       * @brief Assigns the assembly policy, forgetting the cached pattern.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other)
      {
        Parent::operator=(other);
        m_cache = {};
        m_coefficients = {};
        return *this;
      }

      /**
       * @brief Moves the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Integrates the form on the device and fills its matrix.
       * @param res Matrix of the form.
       * @param input Input data.
       */
      void execute(OperatorType& res, const InputType& input) override
      {
        const auto& trial = input.getTrialFES();
        const auto& test = input.getTestFES();
        constexpr bool distributed = Context::IsMPI<typename FormLanguage::Traits<
          typename FormLanguage::Traits<TrialFES>::MeshType>::ContextType>::value;
        // The kernels number the local degrees of freedom, ghosts included.
        size_t rows, cols;
        if constexpr (distributed)
        {
          rows = test.getLocalSize();
          cols = trial.getLocalSize();
        }
        else
        {
          rows = test.getSize();
          cols = trial.getSize();
        }
        SystemType system(trial.getMesh().getContext().getExecutionSpace(), rows, cols);
        if constexpr (distributed)
        {
          system.distribute(
            Internal::COODistribution::of(test), Internal::COODistribution::of(trial));
        }
        for (auto& lbfi : input.getLocalBFIs())
          system.add(lbfi, 0, 0);
        for (auto& gbfi : input.getGlobalBFIs())
          system.add(gbfi, 0, 0);
        system.finish(m_cache, m_coefficients, res, nullptr);
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
      typename SystemType::CacheType m_cache;
      typename SystemType::CoefficientCache m_coefficients;
  };

  /**
   * @ingroup KokkosAssemblySpecializations
   * @brief Device assembly of a PETSc linear form on a Kokkos mesh.
   *
   * The integrators of the form, possibly negated, are integrated on the
   * device and summed into a PETSc Kokkos vector.
   *
   * @tparam FES Test space on a resident Kokkos mesh.
   */
  template <class FES>
  class Kokkos<::Vec, Variational::LinearForm<FES, ::Vec>> final
    : public AssemblyBase<::Vec, Variational::LinearForm<FES, ::Vec>>
  {
    public:
      /// @brief Vector type.
      using VectorType = ::Vec;

      /// @brief Linear form type.
      using LinearFormType = Variational::LinearForm<FES, VectorType>;

      /// @brief Parent class.
      using Parent = AssemblyBase<VectorType, LinearFormType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the test mesh context.
      using BackendType = typename FormLanguage::Traits<
        typename FormLanguage::Traits<FES>::MeshType>::ContextType::BackendType;

      /// @brief Device system type.
      using SystemType = Internal::COOAssembly<BackendType>;

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
          m_cache(std::move(other.m_cache)),
          m_coefficients(std::move(other.m_coefficients))
      {}

      /**
       * @brief Assigns the assembly policy, forgetting the cached buffers.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other)
      {
        Parent::operator=(other);
        m_cache = {};
        m_coefficients = {};
        return *this;
      }

      /**
       * @brief Moves the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Integrates the form on the device and fills its vector.
       * @param res Vector of the form.
       * @param input Input data.
       */
      void execute(VectorType& res, const InputType& input) override
      {
        const auto& fes = input.getFES();
        constexpr bool distributed = Context::IsMPI<typename FormLanguage::Traits<
          typename FormLanguage::Traits<FES>::MeshType>::ContextType>::value;
        PetscInt n = static_cast<PetscInt>(fes.getSize());
        const PetscInt N = n;
        size_t size = fes.getSize();
        if constexpr (distributed)
        {
          // The kernels number the local degrees of freedom, ghosts included.
          size = fes.getLocalSize();
          Index begin = 0, end = 0;
          fes.getOwnershipRange(begin, end);
          n = static_cast<PetscInt>(end - begin);
        }
        SystemType system(fes.getMesh().getContext().getExecutionSpace(), size, 0);
        if constexpr (distributed)
          system.distribute(Internal::COODistribution::of(fes));
        for (auto& lfi : input.getLFIs())
          system.add(lfi, 0, Real(1));
        PetscErrorCode ierr =
          PETSc::Assembly::VectorSetup(res).prepare({n, N, VECKOKKOS, false, false});
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
        system.finish(m_cache, m_coefficients, nullptr, res);
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
      typename SystemType::CacheType m_cache;
      typename SystemType::CoefficientCache m_coefficients;
  };

  /**
   * @ingroup KokkosAssemblySpecializations
   * @brief Device assembly of a single field PETSc problem on a Kokkos mesh.
   *
   * The integrators of the problem, possibly negated, are integrated on the
   * device and handed to PETSc as coordinate lists in device memory. The
   * Dirichlet and periodic conditions are imposed on the device, as the
   * host assemblies do.
   *
   * The entries of the preassembled forms are copied from their PETSc
   * objects and combined with the integrated ones.
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
      using SystemType = Internal::COOAssembly<BackendType>;

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
          m_cache(std::move(other.m_cache)),
          m_coefficients(std::move(other.m_coefficients))
      {}

      /**
       * @brief Assigns the assembly policy, forgetting the cached pattern.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other)
      {
        Parent::operator=(other);
        m_cache = {};
        m_coefficients = {};
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
      /// @brief Whether the mesh is distributed over MPI processes.
      static constexpr bool IsDistributed =
        Context::IsMPI<typename FormLanguage::Traits<typename FormLanguage::Traits<
          typename FormLanguage::Traits<U>::FESType>::MeshType>::ContextType>::value;

      /**
       * @brief Assembles the problem into the linear system.
       * @param axb Linear system to assemble.
       * @param input Problem input.
       * @param doMatrix Whether the matrix is assembled.
       * @param doVector Whether the right-hand side is assembled.
       */
      void assemble(
        LinearSystemType& axb, const InputType& input, bool doMatrix, bool doVector)
      {
        auto& pb = input.getProblemBody();
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();
        const auto& trial = u.getFiniteElementSpace();
        const auto& test = v.getFiniteElementSpace();
        // The kernels number the local degrees of freedom, ghosts included.
        size_t rows, cols;
        if constexpr (IsDistributed)
        {
          rows = test.getLocalSize();
          cols = trial.getLocalSize();
        }
        else
        {
          rows = test.getSize();
          cols = trial.getSize();
        }
        SystemType system(trial.getMesh().getContext().getExecutionSpace(), rows, cols);
        if constexpr (IsDistributed)
        {
          if (&trial.getMesh() != &test.getMesh() || rows != cols)
          {
            Alert::Exception()
              << "A distributed device problem requires the trial and test "
                 "functions to share their space."
              << Alert::Raise;
          }
          system.distribute(Internal::COODistribution::of(trial));
        }
        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() == u.getUUID())
            system.fix(dbc, 0);
        }
        for (auto& pbc : pb.getPBCs())
        {
          if (pbc.getOperand().getUUID() == u.getUUID())
            system.identify(pbc, 0);
        }
        // The matrix is still needed on the right-hand side to lift the
        // constrained values.
        if (doMatrix || system.hasConstraints())
        {
          for (auto& lbfi : pb.getLocalBFIs())
            system.add(lbfi, 0, 0);
          for (auto& gbfi : pb.getGlobalBFIs())
            system.add(gbfi, 0, 0);
          for (const auto& bf : pb.getBFs())
            system.add(bf.getOperator(), 0, 0);
        }
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        if (doVector)
        {
          for (auto& lfi : pb.getLFIs())
            system.add(lfi, 0, Real(-1));
          for (const auto& lf : pb.getLFs())
            system.add(lf.getVector(), 0);
        }
        if constexpr (IsDistributed)
        {
          Index begin = 0, end = 0;
          trial.getOwnershipRange(begin, end);
          const size_t owned = end - begin;
          Internal::prepareSystem(
            axb, owned, test.getSize(), owned, trial.getSize(), doVector);
        }
        else
        {
          Internal::prepareSystem(axb, rows, cols, doVector);
        }
        system.finish(m_cache, m_coefficients, doMatrix ? axb.getOperator() : nullptr,
          doVector ? axb.getVector() : nullptr);
      }

      typename SystemType::CacheType m_cache;
      typename SystemType::CoefficientCache m_coefficients;
  };

  /**
   * @ingroup KokkosAssemblySpecializations
   * @brief Device assembly of a mixed PETSc problem on a Kokkos mesh.
   *
   * Each field occupies a block of the system, at the offsets given by the
   * problem input, and the system is assembled as in the single field case.
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
      using SystemType = Internal::COOAssembly<BackendType>;

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
          m_cache(std::move(other.m_cache)),
          m_coefficients(std::move(other.m_coefficients))
      {}

      /**
       * @brief Assigns the assembly policy, forgetting the cached pattern.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other)
      {
        Parent::operator=(other);
        m_cache = {};
        m_coefficients = {};
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
      /// @brief Whether the mesh of the first field is distributed over MPI processes.
      static constexpr bool IsDistributed =
        Context::IsMPI<typename FormLanguage::Traits<typename FormLanguage::Traits<
          typename FormLanguage::Traits<U1>::FESType>::MeshType>::ContextType>::value;

      /**
       * @brief Assembles the problem into the linear system.
       * @param axb Linear system to assemble.
       * @param input Problem input.
       * @param doMatrix Whether the matrix is assembled.
       * @param doVector Whether the right-hand side is assembled.
       */
      void assemble(
        LinearSystemType& axb, const InputType& input, bool doMatrix, bool doVector)
      {
        auto& pb = input.getProblemBody();
        const auto& us = input.getTrialFunctions();
        const auto& vs = input.getTestFunctions();
        const auto& trialUUIDMap = input.getTrialUUIDMap();
        const auto& testUUIDMap = input.getTestUUIDMap();
        // The kernels number the local degrees of freedom, ghosts included,
        // and the blocks of a distributed system follow each other in the
        // local numbering.
        std::vector<Index> trialOffsets(
          input.getTrialOffsets().begin(), input.getTrialOffsets().end());
        std::vector<Index> testOffsets(
          input.getTestOffsets().begin(), input.getTestOffsets().end());
        size_t rows = input.getTotalTestSize();
        size_t cols = input.getTotalTrialSize();
        Internal::COODistribution distribution;
        if constexpr (IsDistributed)
        {
          std::vector<size_t> sizes(trialOffsets.size());
          std::vector<size_t> globalSizes(trialOffsets.size());
          cols = 0;
          size_t i = 0;
          us.apply([&](const auto& u) {
            const auto& fes = u.get().getFiniteElementSpace();
            trialOffsets[i] = static_cast<Index>(cols);
            sizes[i] = fes.getLocalSize();
            globalSizes[i] = fes.getSize();
            cols += sizes[i];
            distribution.append(fes, static_cast<PetscInt>(input.getTrialOffsets()[i]));
            ++i;
          });
          rows = 0;
          bool same = testOffsets.size() == trialOffsets.size();
          i = 0;
          vs.apply([&](const auto& v) {
            const auto& fes = v.get().getFiniteElementSpace();
            testOffsets[i] = static_cast<Index>(rows);
            rows += fes.getLocalSize();
            same = same && i < sizes.size() && fes.getLocalSize() == sizes[i] &&
              fes.getSize() == globalSizes[i];
            ++i;
          });
          if (!same)
          {
            Alert::Exception()
              << "A distributed device problem requires each test function to "
                 "share the space of the trial function of its block."
              << Alert::Raise;
          }
          distribution.close(us.template get<0>()
                               .get()
                               .getFiniteElementSpace()
                               .getMesh()
                               .getContext()
                               .getCommunicator(),
            static_cast<PetscInt>(input.getTotalTrialSize()));
        }
        const auto trialOffset = [&](const auto& uuid) -> Index {
          const auto it = trialUUIDMap.left.find(uuid);
          assert(it != trialUUIDMap.left.end());
          return trialOffsets[it->second];
        };
        const auto testOffset = [&](const auto& uuid) -> Index {
          const auto it = testUUIDMap.left.find(uuid);
          assert(it != testUUIDMap.left.end());
          return testOffsets[it->second];
        };

        const auto& u1 = us.template get<0>().get();
        SystemType system(
          u1.getFiniteElementSpace().getMesh().getContext().getExecutionSpace(), rows,
          cols);
        if constexpr (IsDistributed)
          system.distribute(std::move(distribution));
        for (auto& dbc : pb.getDBCs())
        {
          const Index offset = trialOffset(dbc.getOperand().getUUID());
          const auto master = dbc.getValueUUID();
          system.fix(dbc, offset, master ? trialOffset(*master) : offset);
        }
        for (auto& pbc : pb.getPBCs())
          system.identify(pbc, trialOffset(pbc.getOperand().getUUID()));
        if (doMatrix || system.hasConstraints())
        {
          for (auto& lbfi : pb.getLocalBFIs())
          {
            system.add(lbfi, testOffset(lbfi.getTestFunction().getUUID()),
              trialOffset(lbfi.getTrialFunction().getUUID()));
          }
          for (auto& gbfi : pb.getGlobalBFIs())
          {
            system.add(gbfi, testOffset(gbfi.getTestFunction().getUUID()),
              trialOffset(gbfi.getTrialFunction().getUUID()));
          }
          for (const auto& bf : pb.getBFs())
          {
            system.add(bf.getOperator(), testOffset(bf.getTestFunction().getUUID()),
              trialOffset(bf.getTrialFunction().getUUID()));
          }
        }
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        if (doVector)
        {
          for (auto& lfi : pb.getLFIs())
            system.add(lfi, testOffset(lfi.getTestFunction().getUUID()), Real(-1));
          for (const auto& lf : pb.getLFs())
            system.add(lf.getVector(), testOffset(lf.getTestFunction().getUUID()));
        }
        if constexpr (IsDistributed)
        {
          const size_t owned = system.getDistribution().getOwnedSize();
          Internal::prepareSystem(axb, owned, input.getTotalTestSize(), owned,
            input.getTotalTrialSize(), doVector);
          axb.setFieldStarts(system.getDistribution().getFieldStarts());
        }
        else
        {
          Internal::prepareSystem(axb, rows, cols, doVector);
        }
        system.finish(m_cache, m_coefficients, doMatrix ? axb.getOperator() : nullptr,
          doVector ? axb.getVector() : nullptr);
      }

      typename SystemType::CacheType m_cache;
      typename SystemType::CoefficientCache m_coefficients;
  };
}

#endif
