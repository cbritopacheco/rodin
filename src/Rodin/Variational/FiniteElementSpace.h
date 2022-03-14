/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENTSPACE_H
#define RODIN_VARIATIONAL_FINITEELEMENTSPACE_H

#include <mfem.hpp>

#include "Rodin/Mesh.h"
#include "Rodin/Utility.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class FiniteElementSpaceBase
   {
      public:
         virtual void update() = 0;

         virtual Mesh& getMesh() = 0;
         virtual const Mesh& getMesh() const = 0;

         virtual int getNumberOfDofs() const = 0;
         virtual int getVectorDimension() const = 0;

         virtual mfem::FiniteElementSpace& getFES() = 0;
         virtual const mfem::FiniteElementSpace& getFES() const = 0;

         virtual mfem::FiniteElementCollection& getFEC() = 0;
         virtual const mfem::FiniteElementCollection& getFEC() const = 0;

         mfem::Array<int> getEssentialTrueDOFs(const std::set<int>& bdrAttr)
         {
            mfem::Array<int> essTrueDofList;
            int maxBdrAttr = *getMesh().getBoundaryAttributes().rbegin();
            getFES().GetEssentialTrueDofs(Utility::set2marker(bdrAttr, maxBdrAttr), essTrueDofList);
            return essTrueDofList;
         }

         mfem::Array<int> getEssentialTrueDOFs(const std::set<int>& bdrAttr, int component)
         {
            mfem::Array<int> essTrueDofList;
            int maxBdrAttr = *getMesh().getBoundaryAttributes().rbegin();
            getFES().GetEssentialTrueDofs(Utility::set2marker(bdrAttr, maxBdrAttr), essTrueDofList, component);
            return essTrueDofList;
         }
   };
}

#endif
