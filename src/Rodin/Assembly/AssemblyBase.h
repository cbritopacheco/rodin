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
  struct BilinearAssemblyInput
  {
    const Geometry::MeshBase& mesh;
    const Variational::FiniteElementSpaceBase& trialFES;
    const Variational::FiniteElementSpaceBase& testFES;
    FormLanguage::List<Variational::BilinearFormIntegratorBase>& bfis;
  };

  template <class OperatorType>
  class AssemblyBase<Variational::BilinearFormBase<OperatorType>>
    : public FormLanguage::Base
  {
    public:
      AssemblyBase() = default;

      AssemblyBase(const AssemblyBase&) = default;

      AssemblyBase(AssemblyBase&&) = default;

      virtual OperatorType execute(const BilinearAssemblyInput& data) const = 0;

      virtual AssemblyBase* copy() const noexcept = 0;
  };

  template <class VectorType>
  class AssemblyBase<Variational::LinearFormBase<VectorType>>
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
