/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_ASSEMBLYBASE_H
#define RODIN_ASSEMBLY_ASSEMBLYBASE_H

#include "Rodin/Math.h"
#include "Rodin/Tuple.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/FormLanguage/List.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "ForwardDecls.h"

namespace Rodin::Assembly
{
  template <class OperatorType, class TrialFES, class TestFES>
  class AssemblyBase<OperatorType, Variational::BilinearForm<TrialFES, TestFES, OperatorType>>
    : public FormLanguage::Base
  {
    public:
      struct Input
      {
        const Geometry::MeshBase& mesh;
        const TrialFES& trialFES;
        const TestFES& testFES;
        FormLanguage::List<Variational::LocalBilinearFormIntegratorBase>& lbfis;
        FormLanguage::List<Variational::GlobalBilinearFormIntegratorBase>& gbfis;
      };

      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase&) = default;

      AssemblyBase(AssemblyBase&&) = default;

      virtual OperatorType execute(const Input& data) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };

  template <class OperatorType, class ... TrialFES, class ... TestFES>
  class AssemblyBase<OperatorType, Tuple<Variational::BilinearForm<TrialFES, TestFES, OperatorType>...>>
  {
    public:
      static_assert(sizeof...(TrialFES) == sizeof...(TestFES));
      using Input =
        Tuple<typename AssemblyBase<OperatorType, Variational::BilinearForm<TrialFES, TestFES, OperatorType>>::Input...>;
  };

  template <class VectorType, class FES>
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
}

#endif
