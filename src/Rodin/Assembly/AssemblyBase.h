/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_ASSEMBLYBASE_H
#define RODIN_ASSEMBLY_ASSEMBLYBASE_H

#include "Rodin/FormLanguage/List.h"

#include "Rodin/Math.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "ForwardDecls.h"

namespace Rodin::Assembly
{
  template <class TrialFES, class TestFES, class OperatorType>
  class AssemblyBase<OperatorType, Variational::BilinearForm<TrialFES, TestFES, OperatorType>>
    : public FormLanguage::Base
  {
    public:
      struct Input
      {
        const Geometry::MeshBase& mesh;
        const TrialFES& trialFES;
        const TestFES& testFES;
        FormLanguage::List<Variational::LocalBilinearFormIntegratorBase>& bfis;
      };

      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase&) = default;

      AssemblyBase(AssemblyBase&&) = default;

      virtual OperatorType execute(const Input& data) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };

  template <class FES, class VectorType>
  class AssemblyBase<VectorType, Variational::LinearForm<FES, VectorType>>
    : public FormLanguage::Base
  {
    public:
      struct Input
      {
        const Geometry::MeshBase& mesh;
        const Variational::FiniteElementSpaceBase& fes;
        FormLanguage::List<Variational::LinearFormIntegratorBase>& lfis;
      };

      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase&) = default;

      AssemblyBase(AssemblyBase&&) = default;

      virtual VectorType execute(const Input& data) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };

  // template <class OperatorType>
  // class AssemblyBase<
  //   OperatorType,
  //   Integral<
  //     Dot<Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
  //         ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>>
  // {
  //   public:
  //     struct Input
  //     {
  //       const Geometry::MeshBase& mesh;
  //       const Variational::FiniteElementSpaceBase& fes;
  //     };

  //     AssemblyBase() = default;

  //     AssemblyBase(const AssemblyBase&) = default;

  //     AssemblyBase(AssemblyBase&&) = default;

  //     virtual  execute(const Input& data) const = 0;

  //     virtual AssemblyBase* copy() const noexcept = 0;
  // };
}

#endif