/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ASSEMBLY_H
#define RODIN_VARIATIONAL_ASSEMBLY_H

#include <variant>
#include <ostream>

#include <mfem.hpp>

#include "Rodin/Geometry/Element.h"
#include "Rodin/FormLanguage/List.h"

#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

namespace Rodin::Variational::Assembly
{
   template <class OperatorType>
   class AssemblyBase<BilinearFormBase<OperatorType>>
      : public FormLanguage::Base
   {
      public:
         struct Input
         {
            const Geometry::MeshBase& mesh;
            const FiniteElementSpaceBase& trialFES;
            const FiniteElementSpaceBase& testFES;
            const FormLanguage::List<BilinearFormIntegratorBase>& bfis;
         };

         AssemblyBase() = default;

         AssemblyBase(const AssemblyBase&) = default;

         AssemblyBase(AssemblyBase&&) = default;

         virtual OperatorType execute(const Input& data) const = 0;

         virtual AssemblyBase* copy() const noexcept = 0;
   };
}

#endif
