/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_ASSEMBLY_KOKKOS_H
#define RODIN_KOKKOS_ASSEMBLY_KOKKOS_H

/**
 * @file
 * @brief Device assembly of bilinear forms, linear forms, problems and
 * Dirichlet conditions on resident Kokkos meshes.
 */

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <KokkosSparse_CooMatrix.hpp>

#if defined(KOKKOS_ENABLE_SYCL)
#include <sycl/ext/intel/math.hpp>
#endif

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Kokkos/Utility.h"
#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Assembly/ConstraintMap.h"
#include "Rodin/Assembly/Sequential.h"
#include "Rodin/Variational/BilinearForm.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/PeriodicBC.h"
#include "Rodin/Variational/Problem.h"
#include "Rodin/Variational/Mult.h"
#include "Rodin/Variational/UnaryMinus.h"

#include "Rodin/Kokkos/ForwardDecls.h"
#include "Rodin/Kokkos/Math/LinearSystem.h"
#include "Rodin/Kokkos/Variational/Integral.h"

#include "SystemAssembly.h"

namespace Rodin::Assembly
{
  /**
   * @ingroup KokkosAssemblySpecializations
   * @brief Device assembly of a sparse bilinear form on a Kokkos mesh.
   *
   * Every local integrator of the form, possibly negated, must implement
   * Rodin::Kokkos::Variational::BilinearFormIntegrator. Each one
   * integrates its cells on the device, and the matrix is summed and
   * compressed on the device before being copied to the host.
   *
   * Global integrators are integrated after the local ones, as on the host.
   *
   * @tparam Solution Solution variable type.
   * @tparam TrialFES Trial space on a resident Kokkos mesh.
   * @tparam TestFES Test space on the same resident Kokkos mesh.
   */
  template <class Matrix, class Solution, class TrialFES, class TestFES>
    requires(std::is_same_v<Matrix, ::Rodin::Math::SparseMatrix<Real>> ||
      requires { typename Matrix::BackendType; })
  class Kokkos<Matrix, Variational::BilinearForm<Solution, TrialFES, TestFES, Matrix>>
    final : public AssemblyBase<Matrix,
              Variational::BilinearForm<Solution, TrialFES, TestFES, Matrix>>
  {
    public:
      /// @brief Operator type.
      using OperatorType = Matrix;

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

      static_assert(
        !Context::IsMPI<
          typename FormLanguage::Traits<TrialFES>::MeshType::ContextType>::value,
        "Native Kokkos algebra is not distributed. Use PETSc algebra with "
        "Rodin/PETSc/Kokkos.h for a distributed assembled system.");

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Integrates the form on the device and builds its matrix.
       * @param res Storage for the assembled operator.
       * @param input Input data.
       */
      void execute(OperatorType& res, const InputType& input) override
      {
        Internal::SystemAssembly<BackendType> system(
          input.getTrialFES().getMesh().getContext().getExecutionSpace(),
          input.getTestFES().getSize(), input.getTrialFES().getSize(), &m_symbolic);
        for (auto& lbfi : input.getLocalBFIs())
          system.add(lbfi, 0, 0);
        for (auto& gbfi : input.getGlobalBFIs())
          system.add(gbfi, 0, 0);
        system.finish(&res, static_cast<::Rodin::Math::Vector<Real>*>(nullptr));
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
      typename Internal::SystemAssembly<BackendType>::Symbolic m_symbolic;
  };

  /**
   * @ingroup KokkosAssemblySpecializations
   * @brief Device assembly of a linear form on a Kokkos mesh.
   *
   * Every integrator of the form, possibly negated, must implement
   * Rodin::Kokkos::Variational::LinearFormIntegrator. Each one
   * integrates its cells on the device, and the vector is summed on the
   * device before being copied to the host.
   *
   * @tparam FES Test space on a resident Kokkos mesh.
   */
  template <class Vector, class FES>
    requires(std::is_same_v<Vector, ::Rodin::Math::Vector<Real>> ||
      requires { typename Vector::MemorySpace; })
  class Kokkos<Vector, Variational::LinearForm<FES, Vector>> final
    : public AssemblyBase<Vector, Variational::LinearForm<FES, Vector>>
  {
    public:
      /// @brief Vector type.
      using VectorType = Vector;

      /// @brief Linear form type.
      using LinearFormType = Variational::LinearForm<FES, VectorType>;

      /// @brief Parent class.
      using Parent = AssemblyBase<VectorType, LinearFormType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the test mesh context.
      using BackendType = typename FormLanguage::Traits<
        typename FormLanguage::Traits<FES>::MeshType>::ContextType::BackendType;

      static_assert(
        !Context::IsMPI<typename FormLanguage::Traits<FES>::MeshType::ContextType>::value,
        "Native Kokkos algebra is not distributed. Use PETSc algebra with "
        "Rodin/PETSc/Kokkos.h for a distributed assembled system.");

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Assigns the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other) = default;

      /**
       * @brief Moves the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Integrates the form on the device and builds its vector.
       * @param res Storage for the assembled vector.
       * @param input Input data.
       */
      void execute(VectorType& res, const InputType& input) override
      {
        Internal::SystemAssembly<BackendType> system(
          input.getFES().getMesh().getContext().getExecutionSpace(),
          input.getFES().getSize(), 0, &m_symbolic);
        for (auto& lfi : input.getLFIs())
          system.add(lfi, 0, Real(1));
        system.finish(nullptr, &res);
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
      typename Internal::SystemAssembly<BackendType>::Symbolic m_symbolic;
  };

  namespace Internal
  {
    /// @brief Whether a matrix is a Rodin::Kokkos::Math::SparseMatrix of reals.
    template <class Matrix>
    struct IsKokkosMatrix : std::false_type
    {};

    /// @brief A Rodin::Kokkos::Math::SparseMatrix of reals.
    template <class Backend>
    struct IsKokkosMatrix<::Rodin::Kokkos::Math::SparseMatrix<Real, Backend>>
      : std::true_type
    {};

    /// @brief Matrix of a system assembled from a Problem on a Kokkos mesh.
    template <class Matrix>
    concept KokkosProblemMatrix =
      std::is_same_v<Matrix, ::Rodin::Math::SparseMatrix<Real>> ||
      IsKokkosMatrix<Matrix>::value;

    /**
     * @brief Builds a system on the host or with its matrix on the device.
     * @param system Contributions of the problem.
     * @param axb System receiving the matrix and the right-hand side.
     */
    template <class Backend, class Matrix, class Vector>
    void finish(
      SystemAssembly<Backend>& system, ::Rodin::Math::LinearSystem<Matrix, Vector>& axb)
    {
      if constexpr (IsKokkosMatrix<Matrix>::value)
      {
        static_assert(
          std::is_same_v<Matrix, typename SystemAssembly<Backend>::BackendMatrix>,
          "The matrix of the system must have the backend of the mesh.");
        system.finish(axb.getOperator(), axb.getVector());
      }
      else
      {
        system.finish(&axb.getOperator(), &axb.getVector());
      }
    }

  }

  /**
   * @ingroup KokkosAssemblySpecializations
   * @brief Device assembly of a single field problem on a Kokkos mesh.
   *
   * The local integrators and linear integrators of the problem, possibly
   * negated, are integrated on the device. Their entries, the preassembled
   * forms and the Dirichlet and periodic conditions are combined on the
   * device, eliminating the constrained degrees of freedom as the host
   * assemblies do. The system is then copied to the host, or only its
   * right-hand side when the matrix stays on the device.
   *
   * @tparam Matrix ::Rodin::Math::SparseMatrix<Real>, copied to the host, or
   * Rodin::Kokkos::Math::SparseMatrix<Real, Backend>, left on the device.
   * @tparam TrialFunction Trial function on a resident Kokkos mesh.
   * @tparam TestFunction Test function on the same mesh.
   */
  template <class Matrix, class Vector, class TrialFunction, class TestFunction>
    requires Internal::KokkosProblemMatrix<Matrix>
  class Kokkos<::Rodin::Math::LinearSystem<Matrix, Vector>,
    Variational::Problem<::Rodin::Math::LinearSystem<Matrix, Vector>, TrialFunction,
      TestFunction>>
    final : public AssemblyBase<::Rodin::Math::LinearSystem<Matrix, Vector>,
              Variational::Problem<::Rodin::Math::LinearSystem<Matrix, Vector>,
                TrialFunction, TestFunction>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType = ::Rodin::Math::LinearSystem<Matrix, Vector>;

      /// @brief Problem type.
      using ProblemType =
        Variational::Problem<LinearSystemType, TrialFunction, TestFunction>;

      /// @brief Parent class.
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the trial mesh context.
      using BackendType = typename FormLanguage::Traits<typename FormLanguage::Traits<
        typename FormLanguage::Traits<TrialFunction>::FESType>::MeshType>::ContextType::
        BackendType;

      static_assert(
        !Context::IsMPI<typename FormLanguage::Traits<typename FormLanguage::Traits<
          TrialFunction>::FESType>::MeshType::ContextType>::value,
        "Native Kokkos algebra is not distributed. Use PETSc algebra with "
        "Rodin/PETSc/Kokkos.h for a distributed assembled system.");

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Assigns the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other) = default;

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
        auto& pb = input.getProblemBody();
        const auto& u = input.getTrialFunction();
        const auto& v = input.getTestFunction();
        Internal::SystemAssembly<BackendType> system(
          u.getFiniteElementSpace().getMesh().getContext().getExecutionSpace(),
          v.getFiniteElementSpace().getSize(), u.getFiniteElementSpace().getSize(),
          &m_symbolic);
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
        for (auto& lbfi : pb.getLocalBFIs())
          system.add(lbfi, 0, 0);
        for (auto& gbfi : pb.getGlobalBFIs())
          system.add(gbfi, 0, 0);
        {
          for (const auto& bf : pb.getBFs())
            system.add(bf.getOperator(), 0, 0);
        }
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        for (auto& lfi : pb.getLFIs())
          system.add(lfi, 0, Real(-1));
        for (const auto& lf : pb.getLFs())
          system.add(lf.getVector(), 0);
        Internal::finish(system, axb);
      }

      /**
       * @brief Assembles the full system and keeps only one side.
       * @param axb Assembled linear system.
       * @param input Problem input.
       * @param target Side of the system to assemble.
       */
      void execute(
        LinearSystemType& axb, const InputType& input, Variational::AssemblyTarget target)
      {
        LinearSystemType scratch;
        execute(scratch, input);
        if (target == Variational::AssemblyTarget::LHS)
          axb.getOperator() = std::move(scratch.getOperator());
        else
          axb.getVector() = std::move(scratch.getVector());
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
      typename Internal::SystemAssembly<BackendType>::Symbolic m_symbolic;
  };

  /**
   * @ingroup KokkosAssemblySpecializations
   * @brief Device assembly of a mixed problem on a Kokkos mesh.
   *
   * Each field occupies a block of the system, at the offsets given by the
   * problem input. Integrators are integrated on the device and their entries
   * are shifted to the blocks of their trial and test functions, then the
   * constraints are eliminated as in the single field assembly. The masters
   * of an identification Dirichlet condition may belong to another field.
   *
   * @tparam Matrix As for the single field assembly.
   */
  template <class Matrix, class Vector, class U1, class U2, class U3, class... Us>
    requires Internal::KokkosProblemMatrix<Matrix>
  class Kokkos<::Rodin::Math::LinearSystem<Matrix, Vector>,
    Variational::Problem<::Rodin::Math::LinearSystem<Matrix, Vector>, U1, U2, U3, Us...>>
    final : public AssemblyBase<::Rodin::Math::LinearSystem<Matrix, Vector>,
              Variational::Problem<::Rodin::Math::LinearSystem<Matrix, Vector>, U1, U2,
                U3, Us...>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType = ::Rodin::Math::LinearSystem<Matrix, Vector>;

      /// @brief Problem type.
      using ProblemType = Variational::Problem<LinearSystemType, U1, U2, U3, Us...>;

      /// @brief Parent class.
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the mesh context of the first field.
      using BackendType = typename FormLanguage::Traits<typename FormLanguage::Traits<
        typename FormLanguage::Traits<U1>::FESType>::MeshType>::ContextType::BackendType;

      static_assert(
        !Context::IsMPI<typename FormLanguage::Traits<
          typename FormLanguage::Traits<U1>::FESType>::MeshType::ContextType>::value,
        "Native Kokkos algebra is not distributed. Use PETSc algebra with "
        "Rodin/PETSc/Kokkos.h for a distributed assembled system.");

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Assigns the assembly policy.
       * @param other Source assembler.
       * @returns Reference to this assembler.
       */
      Kokkos& operator=(const Kokkos& other) = default;

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

        const auto& u1 = input.getTrialFunctions().template get<0>().get();
        Internal::SystemAssembly<BackendType> system(
          u1.getFiniteElementSpace().getMesh().getContext().getExecutionSpace(),
          input.getTotalTestSize(), input.getTotalTrialSize(), &m_symbolic);
        for (auto& dbc : pb.getDBCs())
        {
          const Index offset = trialOffset(dbc.getOperand().getUUID());
          const auto master = dbc.getValueUUID();
          system.fix(dbc, offset, master ? trialOffset(*master) : offset);
        }
        for (auto& pbc : pb.getPBCs())
          system.identify(pbc, trialOffset(pbc.getOperand().getUUID()));
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
        // The problem reads a(u, v) + l(v) = 0, so b = -l. The host adds the
        // linear integrators before the preassembled forms, which matters
        // for the values they lift.
        for (auto& lfi : pb.getLFIs())
          system.add(lfi, testOffset(lfi.getTestFunction().getUUID()), Real(-1));
        {
          for (const auto& bf : pb.getBFs())
          {
            system.add(bf.getOperator(), testOffset(bf.getTestFunction().getUUID()),
              trialOffset(bf.getTrialFunction().getUUID()));
          }
        }
        for (const auto& lf : pb.getLFs())
          system.add(lf.getVector(), testOffset(lf.getTestFunction().getUUID()));
        Internal::finish(system, axb);
      }

      /**
       * @brief Assembles the full system and keeps only one side.
       * @param axb Assembled linear system.
       * @param input Problem input.
       * @param target Side of the system to assemble.
       */
      void execute(
        LinearSystemType& axb, const InputType& input, Variational::AssemblyTarget target)
      {
        LinearSystemType scratch;
        execute(scratch, input);
        if (target == Variational::AssemblyTarget::LHS)
          axb.getOperator() = std::move(scratch.getOperator());
        else
          axb.getVector() = std::move(scratch.getVector());
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
      typename Internal::SystemAssembly<BackendType>::Symbolic m_symbolic;
  };

  /**
   * @ingroup KokkosAssemblySpecializations
   * @brief Device assembly of a value Dirichlet condition on a Lagrange
   * Kokkos space.
   *
   * The constrained faces are the faces of codimension one carrying one of
   * the essential attributes or, when no attribute is given, the boundary
   * faces. Each constrained degree of freedom is assigned to its first face,
   * where the value is evaluated at its node. Grid function values are
   * therefore interpolated from the degrees of freedom of that face, which
   * matches the host for continuous data.
   *
   * The device arrays of an execution are kept for the next one, which
   * reuses them when the space and the constrained degrees of freedom have
   * the same sizes. A copy starts without arrays.
   *
   * @tparam Solution Solution variable type.
   * @tparam FES Device space, whose degrees of freedom provide their nodes.
   * @tparam ValueDerived Value function type.
   */
  template <class Solution, class FES, class ValueDerived>
    requires ::Rodin::Kokkos::Variational::Space<FES>
  class Kokkos<IndexMap<Real>,
    Variational::DirichletBC<Variational::TrialFunction<Solution, FES>,
      Variational::FunctionBase<ValueDerived>>>
    final : public AssemblyBase<IndexMap<Real>,
              Variational::DirichletBC<Variational::TrialFunction<Solution, FES>,
                Variational::FunctionBase<ValueDerived>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Device degrees of freedom.
      using DOFsType = typename FESType::DOFsType;

      /// @brief Dirichlet condition type.
      using DirichletBCType =
        Variational::DirichletBC<Variational::TrialFunction<Solution, FESType>,
          Variational::FunctionBase<ValueDerived>>;

      /// @brief Parent class.
      using Parent = AssemblyBase<IndexMap<Real>, DirichletBCType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Resident mesh type.
      using MeshType = typename FormLanguage::Traits<FESType>::MeshType;

      /// @brief Execution space of the mesh context.
      using ExecutionSpace = typename MeshType::ExecutionSpace;

      /// @brief Memory space of the mesh context.
      using MemorySpace = typename MeshType::MemorySpace;

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Evaluates the condition on the device.
       * @param res Constrained degrees of freedom and their values.
       * @param input Boundary condition input.
       */
      void execute(IndexMap<Real>& res, const InputType& input) override
      {
        const auto& fes = input.getOperand().getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto& essBdr = input.getEssentialBoundary();
        res.clear();
        const Index d = mesh.getDimension();
        if (d == 0)
          return;
        if (essBdr.empty() && !mesh.getConnectivity().hasIncidence(d - 1, d))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Compute face-to-cell incidence before imposing a boundary condition."
            << Alert::Raise;
        }
        const auto& execution = mesh.getContext().getExecutionSpace();
        const DOFsType dofs = fes.getDOFsView();
        const Index stride = dofs.getStride();
        const Index vdim = fes.getVectorDimension();
        // The arrays are kept between the executions: a time loop assembles
        // the same condition at every step, and an allocation or a release
        // waits for the device.
        if (!m_boundary || *m_boundary != essBdr)
        {
          m_attributes = ::Rodin::Kokkos::Geometry::uploadAttributes(execution, essBdr);
          m_boundary = essBdr;
        }
        const auto attributes = m_attributes;
        if (m_owners.extent(0) != static_cast<size_t>(stride))
        {
          m_owners = ::Kokkos::View<Index*, MemorySpace>();
          m_owners = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin boundary owners")),
            stride);
        }
        // Only the values of the constrained degrees of freedom, written by
        // the evaluation, are read.
        if (m_values.extent(0) != static_cast<size_t>(stride * vdim))
        {
          m_values = ::Kokkos::View<Real*, MemorySpace>();
          m_values = ::Kokkos::View<Real*, MemorySpace>(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin boundary values")),
            stride * vdim);
        }
        const auto owners = m_owners;
        const auto values = m_values;
        ::Kokkos::deep_copy(execution, owners, std::numeric_limits<Index>::max());
        const ::Rodin::Kokkos::Variational::View<Variational::FunctionBase<ValueDerived>>
          value(input.getValue());
        evaluate(execution, mesh.getView(), attributes, owners, values, value, dofs,
          mesh.getFaceCount(), vdim);

        // Only the constrained degrees of freedom, a fraction of the space,
        // are brought back to the host.
        const auto range =
          ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>(
            execution, 0, stride);
        Index count = 0;
        ::Kokkos::parallel_reduce(
          "Rodin boundary count", range, BoundaryCount{owners}, count);
        if (m_indices.extent(0) != static_cast<size_t>(count) ||
          m_compact.extent(0) != static_cast<size_t>(count * vdim))
        {
          m_indices = ::Kokkos::View<Index*, MemorySpace>();
          m_compact = ::Kokkos::View<Real*, MemorySpace>();
          m_indices = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin boundary indices")),
            count);
          m_compact = ::Kokkos::View<Real*, MemorySpace>(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin boundary compact values")),
            count * vdim);
          m_hostIndices =
            ::Kokkos::create_mirror_view(::Kokkos::WithoutInitializing, m_indices);
          m_hostValues =
            ::Kokkos::create_mirror_view(::Kokkos::WithoutInitializing, m_compact);
        }
        const auto indices = m_indices;
        const auto compact = m_compact;
        ::Kokkos::parallel_scan("Rodin boundary compaction", range,
          BoundaryCompaction{owners, values, indices, compact, stride, vdim});
        const auto hostIndices = m_hostIndices;
        const auto hostValues = m_hostValues;
        ::Kokkos::deep_copy(execution, hostIndices, indices);
        ::Kokkos::deep_copy(execution, hostValues, compact);
        execution.fence("Rodin boundary values");
        for (Index r = 0; r < vdim; ++r)
        {
          for (Index k = 0; k < count; ++k)
          {
            res.emplace_hint(
              res.end(), hostIndices(k) + r * stride, hostValues(k * vdim + r));
          }
        }
      }

      /// @brief Counts the constrained scalar degrees of freedom.
      struct BoundaryCount
      {
        /// @brief First constrained face of every scalar degree of freedom.
        ::Kokkos::View<Index*, MemorySpace> owners;

        /**
         * @brief Counts a degree of freedom when it is constrained.
         * @param i Scalar degree of freedom.
         * @param count Running count.
         */
        KOKKOS_INLINE_FUNCTION void operator()(Index i, Index& count) const
        {
          if (owners(i) != ::Kokkos::Experimental::finite_max_v<Index>)
            ++count;
        }
      };

      /**
       * @brief Gathers the constrained scalar degrees of freedom and their
       * values, by increasing degree of freedom.
       */
      struct BoundaryCompaction
      {
        /// @brief First constrained face of every scalar degree of freedom.
        ::Kokkos::View<Index*, MemorySpace> owners;

        /// @brief Component r of the value of the degree of freedom i at i + r * stride.
        ::Kokkos::View<Real*, MemorySpace> values;

        /// @brief Constrained scalar degrees of freedom.
        ::Kokkos::View<Index*, MemorySpace> indices;

        /// @brief Component r of the value of the constrained k at k * vdim + r.
        ::Kokkos::View<Real*, MemorySpace> compact;

        /// @brief Offset between components.
        Index stride;

        /// @brief Vector dimension of the space.
        Index vdim;

        /**
         * @brief Stores a constrained degree of freedom at its rank.
         * @param i Scalar degree of freedom.
         * @param position Rank of @p i among the constrained ones.
         * @param final Whether the ranks are final.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          Index i, Index& position, bool final) const
        {
          if (owners(i) == ::Kokkos::Experimental::finite_max_v<Index>)
            return;
          if (final)
          {
            indices(position) = i;
            for (Index r = 0; r < vdim; ++r)
              compact(position * vdim + r) = values(i + r * stride);
          }
          ++position;
        }
      };

      /**
       * @brief Assigns the degrees of freedom of a constrained face to their
       * first constrained face.
       */
      template <class AttributeView>
      struct BoundaryOwners
      {
        /// @brief Device mesh.
        ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase> view;

        /// @brief Essential attributes, empty to select the boundary.
        AttributeView attributes;

        /// @brief Device degrees of freedom.
        DOFsType dofs;

        /// @brief Dimension of the faces.
        Index faceDim;

        /**
         * @brief First constrained face of every scalar degree of freedom,
         * lowered atomically.
         */
        ::Kokkos::View<Index*, MemorySpace> owners;

        /**
         * @brief Claims the degrees of freedom of a face when it is
         * constrained.
         * @param f Face index.
         */
        KOKKOS_INLINE_FUNCTION void operator()(Index f) const
        {
          if (!isConstrained(view, attributes, f))
            return;
          const auto indices = dofs.bind(view.getPolytope(faceDim, f));
          for (Index q = 0; q < indices.size(); ++q)
            ::Kokkos::atomic_min(&owners(indices[q]), f);
        }
      };

      /**
       * @brief Evaluates the value at the nodes of the degrees of freedom
       * owned by a constrained face.
       */
      template <class AttributeView, class ValueView>
      struct BoundaryValues
      {
        /// @brief Device mesh.
        ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase> view;

        /// @brief Essential attributes, empty to select the boundary.
        AttributeView attributes;

        /// @brief Dimension of the faces.
        Index faceDim;

        /// @brief Device degrees of freedom.
        DOFsType dofs;

        /// @brief First constrained face of every scalar degree of freedom.
        ::Kokkos::View<Index*, MemorySpace> owners;

        /**
         * @brief Receives component r of the value of the degree of
         * freedom i at i + r * stride.
         */
        ::Kokkos::View<Real*, MemorySpace> values;

        /// @brief Offset between components.
        Index stride;

        /// @brief Vector dimension of the space.
        Index vdim;

        /// @brief Device view of the value function.
        ValueView value;

        /**
         * @brief Evaluates the value at the degrees of freedom owned by a
         * face when it is constrained.
         * @param f Face index.
         */
        KOKKOS_INLINE_FUNCTION void operator()(Index f) const
        {
          if (!isConstrained(view, attributes, f))
            return;
          const auto face = view.getPolytope(faceDim, f);
          const auto indices = dofs.bind(face);
          for (Index q = 0; q < indices.size(); ++q)
          {
            const Index dof = indices[q];
            if (owners(dof) != f)
              continue;
            const auto p = face.getPoint(dofs.getNode(face, q));
            store(values, dof, stride, vdim, value.getValue(p));
          }
        }
      };

      /**
       * @brief Selects the constrained faces and evaluates the value.
       * @param execution Execution space instance.
       * @param view Device mesh.
       * @param attributes Essential attributes; empty selects the boundary.
       * @param owners First constrained face of every scalar degree of freedom.
       * @param values Value of component r of the degree of freedom i at
       * i + r * stride.
       * @param value Device view of the value function.
       * @param dofs Device degrees of freedom.
       * @param faceCount Number of faces.
       * @param vdim Vector dimension of the space.
       *
       * The counts are given by the host, since the topology of the device
       * mesh is not readable on the host.
       */
      template <class AttributeView, class ValueView>
      static void evaluate(const ExecutionSpace& execution,
        const ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase>& view,
        const AttributeView& attributes,
        const ::Kokkos::View<Index*, MemorySpace>& owners,
        const ::Kokkos::View<Real*, MemorySpace>& values, const ValueView& value,
        const DOFsType& dofs, Index faceCount, Index vdim)
      {
        const Index stride = dofs.getStride();
        const auto policy =
          ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>(
            execution, 0, faceCount);
        const Index faceDim = view.getDimension() - 1;
        ::Kokkos::parallel_for("Rodin boundary owners", policy,
          BoundaryOwners<AttributeView>{view, attributes, dofs, faceDim, owners});
        ::Kokkos::parallel_for("Rodin boundary values", policy,
          BoundaryValues<AttributeView, ValueView>{
            view, attributes, faceDim, dofs, owners, values, stride, vdim, value});
      }

      /**
       * @brief Stores a scalar value of a degree of freedom.
       * @param values Boundary values.
       * @param dof Scalar degree of freedom.
       * @param stride Offset between components.
       * @param vdim Vector dimension of the space.
       * @param v Value at the node.
       */
      KOKKOS_INLINE_FUNCTION static void store(
        const ::Kokkos::View<Real*, MemorySpace>& values, Index dof,
        [[maybe_unused]] Index stride, [[maybe_unused]] Index vdim, Real v)
      {
        values(dof) = v;
      }

      /**
       * @brief Stores the components of a vector value of a degree of freedom.
       * @param values Boundary values, component r of the degree of freedom i
       * at i + r * stride.
       * @param dof Scalar degree of freedom.
       * @param stride Offset between components.
       * @param vdim Vector dimension of the space.
       * @param v Value at the node.
       */
      KOKKOS_INLINE_FUNCTION static void store(
        const ::Kokkos::View<Real*, MemorySpace>& values, Index dof, Index stride,
        Index vdim, const ::Rodin::Kokkos::Geometry::Vector& v)
      {
        for (Index r = 0; r < vdim; ++r)
          values(dof + r * stride) = v(r);
      }

      /**
       * @brief Tests whether a face is constrained.
       * @param view Device mesh.
       * @param attributes Essential attributes; empty selects the boundary.
       * @param f Face index.
       * @return Whether the face carries the condition.
       *
       * On a distributed mesh, a ghost face incident to one local cell may
       * be interior, so the boundary is read from the owned faces.
       */
      template <class AttributeView>
      KOKKOS_INLINE_FUNCTION static bool isConstrained(
        const ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase>& view,
        const AttributeView& attributes, Index f)
      {
        const Index d = view.getDimension();
        if (attributes.extent(0) == 0)
        {
          return view.isOwned(d - 1, f) &&
            view.getConnectivity().getIncidence(d - 1, d, f).size() == 1;
        }
        const auto face = view.getPolytope(d - 1, f);
        return ::Rodin::Kokkos::Geometry::isSelected(face, attributes);
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
      /// @brief Essential attributes of the resident ones, when uploaded.
      Optional<FlatSet<::Rodin::Geometry::Attribute>> m_boundary;

      /// @brief Resident essential attributes.
      ::Kokkos::View<::Rodin::Geometry::Attribute*, MemorySpace> m_attributes;

      /// @brief First constrained face of every scalar degree of freedom.
      ::Kokkos::View<Index*, MemorySpace> m_owners;

      /// @brief Component r of the value of the degree of freedom i at i + r * stride.
      ::Kokkos::View<Real*, MemorySpace> m_values;

      /// @brief Constrained scalar degrees of freedom.
      ::Kokkos::View<Index*, MemorySpace> m_indices;

      /// @brief Component r of the value of the constrained k at k * vdim + r.
      ::Kokkos::View<Real*, MemorySpace> m_compact;

      /// @brief Host copy of the constrained scalar degrees of freedom.
      decltype(::Kokkos::create_mirror_view(::Kokkos::WithoutInitializing,
        std::declval<const ::Kokkos::View<Index*, MemorySpace>&>())) m_hostIndices;

      /// @brief Host copy of the values of the constrained degrees of freedom.
      decltype(::Kokkos::create_mirror_view(::Kokkos::WithoutInitializing,
        std::declval<const ::Kokkos::View<Real*, MemorySpace>&>())) m_hostValues;
  };

  /**
   * @ingroup KokkosAssemblySpecializations
   * @brief Identification Dirichlet condition @f$ u = A(v) @f$ on a Kokkos
   * space.
   *
   * The coefficients are the slave degrees of freedom functionals applied to
   * the master basis, evaluated on the host from the host finite elements of
   * the device space. This is the host algorithm, so the coefficients are
   * identical.
   *
   * @tparam Solution Solution variable type.
   * @tparam FES1 Device space of the slave function.
   * @tparam Derived2 Type of @f$ A(v) @f$.
   * @tparam FES2 Space of the master function.
   * @tparam Sp Shape function space type of @f$ A(v) @f$.
   */
  template <class Solution, class FES1, class Derived2, class FES2,
    Variational::ShapeFunctionSpaceType Sp>
    requires ::Rodin::Kokkos::Variational::Space<FES1>
  class Kokkos<IndexMap<std::pair<IndexArray, ::Rodin::Math::Vector<Real>>>,
    Variational::DirichletBC<Variational::TrialFunction<Solution, FES1>,
      Variational::ShapeFunctionBase<Derived2, FES2, Sp>>>
    final
    : public AssemblyBase<IndexMap<std::pair<IndexArray, ::Rodin::Math::Vector<Real>>>,
        Variational::DirichletBC<Variational::TrialFunction<Solution, FES1>,
          Variational::ShapeFunctionBase<Derived2, FES2, Sp>>>
  {
    public:
      /// @brief Slave to master coefficients.
      using OutputType = IndexMap<std::pair<IndexArray, ::Rodin::Math::Vector<Real>>>;

      /// @brief Dirichlet condition type.
      using DirichletBCType =
        Variational::DirichletBC<Variational::TrialFunction<Solution, FES1>,
          Variational::ShapeFunctionBase<Derived2, FES2, Sp>>;

      /// @brief Parent class.
      using Parent = AssemblyBase<OutputType, DirichletBCType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

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
        : Parent(std::move(other))
      {}

      /**
       * @brief Computes the identification coefficients.
       * @param res Slave degrees of freedom and their master coefficients.
       * @param input Boundary condition input.
       */
      void execute(OutputType& res, const InputType& input) override
      {
        Sequential<OutputType, DirichletBCType>().execute(res, input);
      }

      /**
       * @brief Copies the assembly.
       * @return Copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }
  };

  /**
   * @brief Default assembly for a test space on a Kokkos mesh.
   * @tparam Host Host context of the mesh, local or distributed.
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Host, class Backend>
  class Default<Context::Kokkos<Host, Backend>>
  {
    public:
      /**
       * @brief Device assembly of the object.
       * @tparam LinearAlgebraType Target linear algebra type.
       * @tparam Object Variational form type to assemble.
       */
      template <class LinearAlgebraType, class Object>
      using Type = Kokkos<LinearAlgebraType, Object>;
  };

  /**
   * @brief Default assembly for trial and test spaces on Kokkos meshes.
   * @tparam Host Host context of both meshes, local or distributed.
   * @tparam Backend Kokkos backend of both mesh contexts.
   */
  template <class Host, class Backend>
  class Default<Context::Kokkos<Host, Backend>, Context::Kokkos<Host, Backend>>
  {
    public:
      /**
       * @brief Device assembly of the object.
       * @tparam LinearAlgebraType Target linear algebra type.
       * @tparam Object Variational form type to assemble.
       */
      template <class LinearAlgebraType, class Object>
      using Type = Kokkos<LinearAlgebraType, Object>;
  };
}

#endif
