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
 * @brief Device assembly of bilinear and linear forms on resident Kokkos
 * meshes.
 */

#include <vector>

#include <Kokkos_Core.hpp>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Assembly/AssemblyBase.h"
#include "Rodin/Variational/BilinearForm.h"
#include "Rodin/Variational/LinearForm.h"
#include "Rodin/Variational/UnaryMinus.h"

#include "Rodin/Kokkos/ForwardDecls.h"
#include "Rodin/Kokkos/Variational/Integral.h"

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

      /// @brief Device integrator interface.
      using DeviceIntegratorType =
        ::Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<BackendType>;

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
          add(triplets, lbfi, Real(1));
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

    private:
      /**
       * @brief Integrates one integrator on the device and appends its entries.
       * @param triplets Host triplets.
       * @param lbfi Local integrator, possibly negated.
       * @param sign Sign of the contributions.
       */
      void add(std::vector<Math::SparseTriplet<Real>>& triplets,
        const Variational::LocalBilinearFormIntegratorBase<Real>& lbfi, Real sign) const
      {
        using Negated =
          Variational::UnaryMinus<Variational::LocalBilinearFormIntegratorBase<Real>>;
        if (const auto* negated = dynamic_cast<const Negated*>(&lbfi))
        {
          add(triplets, negated->getOperand(), -sign);
          return;
        }
        const auto* integrator = dynamic_cast<const DeviceIntegratorType*>(&lbfi);
        if (!integrator)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Local integrator has no device implementation." << Alert::Raise;
        }
        const auto device = integrator->assemble();
        const auto rows =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), device.rows);
        const auto cols =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), device.cols);
        const auto values =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), device.values);
        triplets.reserve(triplets.size() + values.extent(0));
        for (size_t e = 0; e < values.extent(0); ++e)
        {
          if (values(e) != Real(0))
            triplets.emplace_back(rows(e), cols(e), sign * values(e));
        }
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

      /// @brief Device integrator interface.
      using DeviceIntegratorType =
        ::Rodin::Kokkos::Variational::DeviceLinearFormIntegrator<BackendType>;

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
          add(res, lfi, Real(1));
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
      /**
       * @brief Integrates one integrator on the device and sums its entries.
       * @param res Assembled vector.
       * @param lfi Integrator, possibly negated.
       * @param sign Sign of the contributions.
       */
      void add(VectorType& res, const Variational::LinearFormIntegratorBase<Real>& lfi,
        Real sign) const
      {
        using Negated =
          Variational::UnaryMinus<Variational::LinearFormIntegratorBase<Real>>;
        if (const auto* negated = dynamic_cast<const Negated*>(&lfi))
        {
          add(res, negated->getOperand(), -sign);
          return;
        }
        const auto* integrator = dynamic_cast<const DeviceIntegratorType*>(&lfi);
        if (!integrator)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Integrator has no device implementation." << Alert::Raise;
        }
        const auto device = integrator->assemble();
        const auto indices =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), device.indices);
        const auto values =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), device.values);
        for (size_t e = 0; e < values.extent(0); ++e)
          res.coeffRef(indices(e)) += sign * values(e);
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
