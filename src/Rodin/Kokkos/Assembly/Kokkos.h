/*
 *          Copyright Carlos BRITO PACHECO 2026.
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
#include <limits>
#include <utility>
#include <variant>
#include <vector>

#include <Kokkos_Core.hpp>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Assembly/ConstraintMap.h"
#include "Rodin/Variational/BilinearForm.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/Problem.h"
#include "Rodin/Variational/UnaryMinus.h"

#include "Rodin/Kokkos/ForwardDecls.h"
#include "Rodin/Kokkos/Variational/Integral.h"

namespace Rodin::Assembly::Internal
{
  /**
   * @brief Finds the device implementation of a local integrator.
   * @tparam Backend Kokkos backend of the mesh context.
   * @param lbfi Local integrator, possibly negated.
   * @param sign Sign of the contributions.
   * @return Device integrator and the sign of its contributions.
   */
  template <class Backend>
  std::pair<const ::Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>*,
    Real>
  getDeviceIntegrator(
    const Variational::LocalBilinearFormIntegratorBase<Real>& lbfi, Real sign = Real(1))
  {
    using Negated =
      Variational::UnaryMinus<Variational::LocalBilinearFormIntegratorBase<Real>>;
    if (const auto* negated = dynamic_cast<const Negated*>(&lbfi))
      return getDeviceIntegrator<Backend>(negated->getOperand(), -sign);
    const auto* integrator = dynamic_cast<
      const ::Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>*>(&lbfi);
    if (!integrator)
      Alert::Exception() << "Local integrator has no device implementation."
                         << Alert::Raise;
    return {integrator, sign};
  }

  /**
   * @brief Finds the device implementation of a linear integrator.
   * @tparam Backend Kokkos backend of the mesh context.
   * @param lfi Integrator, possibly negated.
   * @param sign Sign of the contributions.
   * @return Device integrator and the sign of its contributions.
   */
  template <class Backend>
  std::pair<const ::Rodin::Kokkos::Variational::DeviceLinearFormIntegrator<Backend>*,
    Real>
  getDeviceIntegrator(
    const Variational::LinearFormIntegratorBase<Real>& lfi, Real sign = Real(1))
  {
    using Negated = Variational::UnaryMinus<Variational::LinearFormIntegratorBase<Real>>;
    if (const auto* negated = dynamic_cast<const Negated*>(&lfi))
      return getDeviceIntegrator<Backend>(negated->getOperand(), -sign);
    const auto* integrator = dynamic_cast<
      const ::Rodin::Kokkos::Variational::DeviceLinearFormIntegrator<Backend>*>(&lfi);
    if (!integrator)
      Alert::Exception() << "Integrator has no device implementation." << Alert::Raise;
    return {integrator, sign};
  }

  /**
   * @brief Integrates a local integrator on the device and visits its
   * nonzero entries on the host.
   * @tparam Backend Kokkos backend of the mesh context.
   * @param lbfi Local integrator, possibly negated.
   * @param sign Sign of the contributions.
   * @param f Called as f(row, col, value), with indices local to the test
   * and trial spaces.
   */
  template <class Backend, class F>
  void forEachDeviceEntry(
    const Variational::LocalBilinearFormIntegratorBase<Real>& lbfi, Real sign, F&& f)
  {
    const auto [integrator, s] = getDeviceIntegrator<Backend>(lbfi, sign);
    const auto device = integrator->assemble();
    const auto rows =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), device.rows);
    const auto cols =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), device.cols);
    const auto values =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), device.values);
    for (size_t e = 0; e < values.extent(0); ++e)
    {
      if (values(e) != Real(0))
        f(rows(e), cols(e), s * values(e));
    }
  }

  /**
   * @brief Integrates a linear integrator on the device and visits its
   * nonzero entries on the host.
   * @tparam Backend Kokkos backend of the mesh context.
   * @param lfi Integrator, possibly negated.
   * @param sign Sign of the contributions.
   * @param f Called as f(index, value), with indices local to the test space.
   */
  template <class Backend, class F>
  void forEachDeviceEntry(
    const Variational::LinearFormIntegratorBase<Real>& lfi, Real sign, F&& f)
  {
    const auto [integrator, s] = getDeviceIntegrator<Backend>(lfi, sign);
    const auto device = integrator->assemble();
    const auto indices =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), device.indices);
    const auto values =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), device.values);
    for (size_t e = 0; e < values.extent(0); ++e)
    {
      if (values(e) != Real(0))
        f(indices(e), s * values(e));
    }
  }

  /**
   * @brief Host linear system receiving device contributions, with value
   * constraints eliminated.
   *
   * Rows of fixed degrees of freedom become identity rows holding the fixed
   * value, and the columns of fixed degrees of freedom move to the
   * right-hand side, exactly as the host problem assemblies eliminate value
   * Dirichlet conditions. Identification constraints are not supported.
   */
  class ConstrainedSystem
  {
    public:
      /**
       * @brief Creates an empty system.
       * @param rows Number of rows, the size of the test spaces.
       * @param cols Number of columns, the size of the trial spaces.
       */
      ConstrainedSystem(size_t rows, size_t cols)
        : m_rows(rows),
          m_cols(cols),
          m_constraints(std::max(rows, cols)),
          m_vector(Math::Vector<Real>::Zero(rows))
      {}

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
          m_constraints.setFixed(static_cast<Index>(offset + local), value);
      }

      /**
       * @brief Adds a matrix contribution.
       * @param row Row index. @param col Column index. @param value Contribution.
       */
      void addMatrix(Index row, Index col, Real value)
      {
        m_triplets.emplace_back(row, col, value);
      }

      /**
       * @brief Adds a right-hand side contribution.
       * @param row Row index. @param value Contribution.
       */
      void addVector(Index row, Real value)
      {
        m_vector.coeffRef(row) += value;
      }

      /**
       * @brief Eliminates the constraints and builds the system.
       * @param A Assembled operator.
       * @param b Assembled right-hand side.
       */
      void finish(Math::SparseMatrix<Real>& A, Math::Vector<Real>& b)
      {
        size_t write = 0;
        for (size_t read = 0; read < m_triplets.size(); ++read)
        {
          const auto& t = m_triplets[read];
          const Index row = static_cast<Index>(t.row());
          const Index col = static_cast<Index>(t.col());
          if (m_constraints.isFixed(row))
            continue;
          if (m_constraints.isFixed(col) && row != col)
          {
            m_vector.coeffRef(row) -= t.value() * m_constraints.getFixedValue(col);
            continue;
          }
          if (write != read)
            m_triplets[write] = t;
          ++write;
        }
        m_triplets.resize(write);
        for (Index i = 0; i < static_cast<Index>(m_rows); ++i)
        {
          if (m_constraints.isFixed(i))
          {
            m_triplets.emplace_back(i, i, Real(1));
            m_vector.coeffRef(i) = m_constraints.getFixedValue(i);
          }
        }
        A.resize(m_rows, m_cols);
        A.setFromTriplets(m_triplets.begin(), m_triplets.end());
        b = std::move(m_vector);
      }

    private:
      size_t m_rows;
      size_t m_cols;
      ConstraintMap<Real> m_constraints;
      Math::Vector<Real> m_vector;
      std::vector<Math::SparseTriplet<Real>> m_triplets;
  };

  /**
   * @brief Rejects the problem terms without a device implementation.
   * @param pb Problem body.
   */
  template <class ProblemBody>
  void requireDeviceProblem(const ProblemBody& pb)
  {
    if (pb.getGlobalBFIs().size() > 0)
    {
      Alert::Exception() << "Global integrators are not supported on Kokkos meshes."
                         << Alert::Raise;
    }
    if (pb.getPBCs().size() > 0)
    {
      Alert::Exception() << "Periodic conditions are not supported on Kokkos meshes."
                         << Alert::Raise;
    }
  }
}

namespace Rodin::Assembly
{
  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a sparse bilinear form on a Kokkos mesh.
   *
   * Every local integrator of the form, possibly negated, must implement
   * Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator. Each one
   * integrates its cells on the device and returns its local matrices as
   * triplets, which are copied to the host and summed into the sparse
   * matrix.
   *
   * Global integrators are not supported.
   *
   * @tparam Solution Solution variable type.
   * @tparam TrialFES Trial space on a resident Kokkos mesh.
   * @tparam TestFES Test space on the same resident Kokkos mesh.
   */
  template <class Solution, class TrialFES, class TestFES>
  class Kokkos<Math::SparseMatrix<Real>,
    Variational::BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Real>>>
    final
    : public AssemblyBase<Math::SparseMatrix<Real>,
        Variational::BilinearForm<Solution, TrialFES, TestFES, Math::SparseMatrix<Real>>>
  {
    public:
      /// @brief Operator type.
      using OperatorType = Math::SparseMatrix<Real>;

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
        if (input.getGlobalBFIs().size() > 0)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Global integrators are not supported on Kokkos meshes." << Alert::Raise;
        }
        std::vector<Math::SparseTriplet<Real>> triplets;
        for (const auto& lbfi : input.getLocalBFIs())
        {
          Internal::forEachDeviceEntry<BackendType>(
            lbfi, Real(1), [&](Index row, Index col, Real value) {
              triplets.emplace_back(row, col, value);
            });
        }
        res.resize(input.getTestFES().getSize(), input.getTrialFES().getSize());
        res.setFromTriplets(triplets.begin(), triplets.end());
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a linear form on a Kokkos mesh.
   *
   * Every integrator of the form, possibly negated, must implement
   * Rodin::Kokkos::Variational::DeviceLinearFormIntegrator. Each one
   * integrates its cells on the device and returns its local vectors as
   * entries, which are copied to the host and summed into the vector.
   *
   * @tparam FES Test space on a resident Kokkos mesh.
   */
  template <class FES>
  class Kokkos<Math::Vector<Real>, Variational::LinearForm<FES, Math::Vector<Real>>> final
    : public AssemblyBase<Math::Vector<Real>,
        Variational::LinearForm<FES, Math::Vector<Real>>>
  {
    public:
      /// @brief Vector type.
      using VectorType = Math::Vector<Real>;

      /// @brief Linear form type.
      using LinearFormType = Variational::LinearForm<FES, VectorType>;

      /// @brief Parent class.
      using Parent = AssemblyBase<VectorType, LinearFormType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the test mesh context.
      using BackendType = typename FormLanguage::Traits<
        typename FormLanguage::Traits<FES>::MeshType>::ContextType::BackendType;

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
        res.resize(input.getFES().getSize());
        res.setZero();
        for (const auto& lfi : input.getLFIs())
        {
          Internal::forEachDeviceEntry<BackendType>(
            lfi, Real(1), [&](Index i, Real value) { res.coeffRef(i) += value; });
        }
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a single field problem on a Kokkos mesh.
   *
   * The local integrators and linear integrators of the problem, possibly
   * negated, are integrated on the device. Their entries, the preassembled
   * forms and the value Dirichlet conditions are combined on the host,
   * eliminating the constrained degrees of freedom as the host assemblies
   * do.
   *
   * Global integrators, periodic conditions and identification conditions
   * are not supported.
   *
   * @tparam TrialFunction Trial function on a resident Kokkos mesh.
   * @tparam TestFunction Test function on the same mesh.
   */
  template <class TrialFunction, class TestFunction>
  class Kokkos<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
    Variational::Problem<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
      TrialFunction, TestFunction>>
    final : public AssemblyBase<
              Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
              Variational::Problem<
                Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
                TrialFunction, TestFunction>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType =
        Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

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
        Internal::requireDeviceProblem(pb);
        Internal::ConstrainedSystem system(
          v.getFiniteElementSpace().getSize(), u.getFiniteElementSpace().getSize());
        for (auto& dbc : pb.getDBCs())
        {
          if (dbc.getOperand().getUUID() == u.getUUID())
            system.fix(dbc, 0);
        }
        for (const auto& lbfi : pb.getLocalBFIs())
        {
          Internal::forEachDeviceEntry<BackendType>(lbfi, Real(1),
            [&](Index row, Index col, Real value) { system.addMatrix(row, col, value); });
        }
        for (const auto& bf : pb.getBFs())
        {
          const auto& op = bf.getOperator();
          for (int k = 0; k < op.outerSize(); ++k)
          {
            for (typename Math::SparseMatrix<Real>::InnerIterator it(op, k); it; ++it)
              system.addMatrix(it.row(), it.col(), it.value());
          }
        }
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        for (const auto& lfi : pb.getLFIs())
        {
          Internal::forEachDeviceEntry<BackendType>(
            lfi, Real(-1), [&](Index i, Real value) { system.addVector(i, value); });
        }
        for (const auto& lf : pb.getLFs())
        {
          const auto& vec = lf.getVector();
          for (Eigen::Index i = 0; i < vec.size(); ++i)
            system.addVector(i, vec.coeff(i));
        }
        system.finish(axb.getOperator(), axb.getVector());
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
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a mixed problem on a Kokkos mesh.
   *
   * Each field occupies a block of the system, at the offsets given by the
   * problem input. Integrators are integrated on the device and their entries
   * are shifted to the blocks of their trial and test functions, then the
   * constraints are eliminated as in the single field assembly.
   *
   * Global integrators, periodic conditions and identification conditions
   * are not supported.
   */
  template <class U1, class U2, class U3, class... Us>
  class Kokkos<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
    Variational::Problem<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
      U1, U2, U3, Us...>>
    final : public AssemblyBase<
              Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>,
              Variational::Problem<
                Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>, U1, U2,
                U3, Us...>>
  {
    public:
      /// @brief Linear system type.
      using LinearSystemType =
        Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

      /// @brief Problem type.
      using ProblemType = Variational::Problem<LinearSystemType, U1, U2, U3, Us...>;

      /// @brief Parent class.
      using Parent = AssemblyBase<LinearSystemType, ProblemType>;

      /// @brief Input type.
      using InputType = typename Parent::InputType;

      /// @brief Kokkos backend of the mesh context of the first field.
      using BackendType = typename FormLanguage::Traits<typename FormLanguage::Traits<
        typename FormLanguage::Traits<U1>::FESType>::MeshType>::ContextType::BackendType;

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

        Internal::requireDeviceProblem(pb);
        Internal::ConstrainedSystem system(
          input.getTotalTestSize(), input.getTotalTrialSize());
        for (auto& dbc : pb.getDBCs())
          system.fix(dbc, trialOffset(dbc.getOperand().getUUID()));
        for (const auto& lbfi : pb.getLocalBFIs())
        {
          const Index row = testOffset(lbfi.getTestFunction().getUUID());
          const Index col = trialOffset(lbfi.getTrialFunction().getUUID());
          Internal::forEachDeviceEntry<BackendType>(
            lbfi, Real(1), [&](Index i, Index j, Real value) {
              system.addMatrix(row + i, col + j, value);
            });
        }
        for (const auto& bf : pb.getBFs())
        {
          const Index row = testOffset(bf.getTestFunction().getUUID());
          const Index col = trialOffset(bf.getTrialFunction().getUUID());
          const auto& op = bf.getOperator();
          for (int k = 0; k < op.outerSize(); ++k)
          {
            for (typename Math::SparseMatrix<Real>::InnerIterator it(op, k); it; ++it)
              system.addMatrix(row + it.row(), col + it.col(), it.value());
          }
        }
        // The problem reads a(u, v) + l(v) = 0, so b = -l.
        for (const auto& lfi : pb.getLFIs())
        {
          const Index row = testOffset(lfi.getTestFunction().getUUID());
          Internal::forEachDeviceEntry<BackendType>(lfi, Real(-1),
            [&](Index i, Real value) { system.addVector(row + i, value); });
        }
        for (const auto& lf : pb.getLFs())
        {
          const Index row = testOffset(lf.getTestFunction().getUUID());
          const auto& vec = lf.getVector();
          for (Eigen::Index i = 0; i < vec.size(); ++i)
            system.addVector(row + i, vec.coeff(i));
        }
        system.finish(axb.getOperator(), axb.getVector());
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
  };

  /**
   * @ingroup RodinAssembly
   * @brief Device assembly of a value Dirichlet condition on a P1 Kokkos space.
   *
   * The constrained faces are the faces of codimension one carrying one of
   * the essential attributes or, when no attribute is given, the boundary
   * faces. Each constrained vertex is assigned to its first face, where the
   * value is evaluated at the vertex. Grid function values are therefore
   * interpolated from the vertices of that face, which matches the host for
   * continuous P1 data.
   *
   * @tparam Solution Solution variable type.
   * @tparam Range Range type of the space, Real or a spatial vector.
   * @tparam Backend Kokkos backend of the mesh context.
   * @tparam ValueDerived Value function type.
   */
  template <class Solution, class Range, class Backend, class ValueDerived>
  class Kokkos<IndexMap<Real>,
    Variational::DirichletBC<
      Variational::TrialFunction<Solution,
        ::Rodin::Kokkos::Variational::DeviceP1Space<Range, Backend>>,
      Variational::FunctionBase<ValueDerived>>>
    final : public AssemblyBase<IndexMap<Real>,
              Variational::DirichletBC<
                Variational::TrialFunction<Solution,
                  ::Rodin::Kokkos::Variational::DeviceP1Space<Range, Backend>>,
                Variational::FunctionBase<ValueDerived>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = ::Rodin::Kokkos::Variational::DeviceP1Space<Range, Backend>;

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
        const Index vertexCount = mesh.getVertexCount();
        const Index vdim = fes.getVectorDimension();
        const auto attributes =
          ::Rodin::Kokkos::Variational::Internal::uploadAttributes(execution, essBdr);
        ::Kokkos::View<Index*, MemorySpace> owners(
          ::Kokkos::view_alloc(execution, std::string("Rodin boundary owners")),
          vertexCount);
        ::Kokkos::deep_copy(execution, owners, std::numeric_limits<Index>::max());
        ::Kokkos::View<Real*, MemorySpace> values(
          ::Kokkos::view_alloc(execution, std::string("Rodin boundary values")),
          vertexCount * vdim);
        const ::Rodin::Kokkos::Variational::DeviceView<
          Variational::FunctionBase<ValueDerived>>
          value(input.getValue());
        evaluate(execution, mesh.getDeviceView(), attributes, owners, values, value,
          mesh.getFaceCount(), vertexCount, vdim);

        const auto host = ::Kokkos::create_mirror_view(owners);
        const auto hostValues = ::Kokkos::create_mirror_view(values);
        ::Kokkos::deep_copy(execution, host, owners);
        ::Kokkos::deep_copy(execution, hostValues, values);
        execution.fence("Rodin boundary values");
        for (Index r = 0; r < vdim; ++r)
        {
          for (Index i = 0; i < vertexCount; ++i)
          {
            if (host(i) != std::numeric_limits<Index>::max())
              res.insert(
                res.end(), {i + r * vertexCount, hostValues(i + r * vertexCount)});
          }
        }
      }

      /**
       * @brief Selects the constrained faces and evaluates the value.
       * @param execution Execution space instance.
       * @param view Device mesh.
       * @param attributes Essential attributes; empty selects the boundary.
       * @param owners First constrained face of every vertex.
       * @param values Value of component r of vertex i at i + r * vertexCount.
       * @param value Device view of the value function.
       * @param faceCount Number of faces.
       * @param vertexCount Number of vertices.
       * @param vdim Vector dimension of the space.
       *
       * The counts are given by the host, since the topology of the device
       * mesh is not readable on the host. Public as the enclosing function of
       * device lambdas.
       */
      template <class AttributeView, class ValueView>
      static void evaluate(const ExecutionSpace& execution,
        const ::Rodin::Kokkos::Geometry::MeshView& view, const AttributeView& attributes,
        const ::Kokkos::View<Index*, MemorySpace>& owners,
        const ::Kokkos::View<Real*, MemorySpace>& values, const ValueView& value,
        Index faceCount, Index vertexCount, Index vdim)
      {
        const auto policy =
          ::Kokkos::RangePolicy<ExecutionSpace, ::Kokkos::IndexType<Index>>(
            execution, 0, faceCount);
        const Index faceDim = view.getDimension() - 1;
        ::Kokkos::parallel_for(
          "Rodin boundary owners", policy, KOKKOS_LAMBDA(Index f) {
            if (!isConstrained(view, attributes, f))
              return;
            const auto vertices = view.getPolytope(faceDim, f).getVertices();
            for (Index a = 0; a < vertices.size(); ++a)
              ::Kokkos::atomic_min(&owners(vertices[a]), f);
          });
        ::Kokkos::parallel_for(
          "Rodin boundary values", policy, KOKKOS_LAMBDA(Index f) {
            if (!isConstrained(view, attributes, f))
              return;
            const auto face = view.getPolytope(faceDim, f);
            const auto vertices = face.getVertices();
            for (Index a = 0; a < vertices.size(); ++a)
            {
              const Index vertex = vertices[a];
              if (owners(vertex) != f)
                continue;
              const auto p = face.getPoint(getReferenceVertex(face, a));
              store(values, vertex, vertexCount, vdim, value.getValue(p));
            }
          });
      }

      /**
       * @brief Stores a scalar value of a vertex.
       * @param values Boundary values.
       * @param vertex Vertex index.
       * @param vertexCount Number of vertices.
       * @param vdim Vector dimension of the space.
       * @param v Value at the vertex.
       */
      KOKKOS_INLINE_FUNCTION static void store(
        const ::Kokkos::View<Real*, MemorySpace>& values, Index vertex,
        [[maybe_unused]] Index vertexCount, [[maybe_unused]] Index vdim, Real v)
      {
        values(vertex) = v;
      }

      /**
       * @brief Stores the components of a vector value of a vertex.
       * @param values Boundary values, component r of vertex i at
       * i + r * vertexCount.
       * @param vertex Vertex index.
       * @param vertexCount Number of vertices.
       * @param vdim Vector dimension of the space.
       * @param v Value at the vertex.
       */
      KOKKOS_INLINE_FUNCTION static void store(
        const ::Kokkos::View<Real*, MemorySpace>& values, Index vertex, Index vertexCount,
        Index vdim, const ::Rodin::Kokkos::Geometry::Vector& v)
      {
        for (Index r = 0; r < vdim; ++r)
          values(vertex + r * vertexCount) = v(r);
      }

      /**
       * @brief Tests whether a face is constrained.
       * @param view Device mesh.
       * @param attributes Essential attributes; empty selects the boundary.
       * @param f Face index.
       * @return Whether the face carries the condition.
       */
      template <class AttributeView>
      KOKKOS_INLINE_FUNCTION static bool isConstrained(
        const ::Rodin::Kokkos::Geometry::MeshView& view, const AttributeView& attributes,
        Index f)
      {
        const Index d = view.getDimension();
        if (attributes.extent(0) == 0)
          return view.getConnectivity().getIncidence(d - 1, d, f).size() == 1;
        const auto face = view.getPolytope(d - 1, f);
        return ::Rodin::Kokkos::Variational::isSelected(face, attributes);
      }

      /**
       * @brief Reference coordinates of a vertex of a polytope.
       * @param polytope Device polytope.
       * @param a Local vertex index.
       * @return Reference coordinates of the P1 node of local vertex a.
       */
      KOKKOS_INLINE_FUNCTION static ::Rodin::Kokkos::Geometry::Vector getReferenceVertex(
        const ::Rodin::Kokkos::Geometry::PolytopeView& polytope, Index a)
      {
        const Index rdim = polytope.getDimension();
        ::Rodin::Kokkos::Geometry::Vector result(rdim);
        if (polytope.getGeometry() == Geometry::Polytope::Type::Quadrilateral)
        {
          result(0) = (a == 1 || a == 2) ? Real(1) : Real(0);
          result(1) = (a == 2 || a == 3) ? Real(1) : Real(0);
        }
        else if (a > 0)
        {
          result(a - 1) = Real(1);
        }
        return result;
      }

      /**
       * @brief Creates a polymorphic copy.
       * @return Pointer to a new copy.
       */
      Kokkos* copy() const noexcept override
      {
        return new Kokkos(*this);
      }
  };

  /**
   * @brief Default assembly for a test space on a Kokkos mesh.
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class Default<Context::Kokkos<Context::Local, Backend>>
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
   * @tparam Backend Kokkos backend of both mesh contexts.
   */
  template <class Backend>
  class Default<Context::Kokkos<Context::Local, Backend>,
    Context::Kokkos<Context::Local, Backend>>
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
