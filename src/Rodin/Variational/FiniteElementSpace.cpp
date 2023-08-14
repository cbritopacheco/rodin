/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
  size_t FiniteElementSpaceBase::getVectorDimension() const
  {
   return getHandle().GetVDim();
  }

  mfem::Array<int> FiniteElementSpaceBase::getEssentialTrueDOFs(
    const std::set<Geometry::Attribute>& bdrAttr) const
  {
   mfem::Array<int> essTrueDofList;
   int maxBdrAttr = *getMesh().getBoundaryAttributes().rbegin();
   getHandle().GetEssentialTrueDofs(Utility::set2marker(bdrAttr, maxBdrAttr), essTrueDofList);
   return essTrueDofList;
  }

  mfem::Array<int> FiniteElementSpaceBase::getEssentialTrueDOFs(
    const std::set<Geometry::Attribute>& bdrAttr, size_t component) const
  {
   mfem::Array<int> essTrueDofList;
   int maxBdrAttr = *getMesh().getBoundaryAttributes().rbegin();
   getHandle().GetEssentialTrueDofs(Utility::set2marker(bdrAttr, maxBdrAttr), essTrueDofList, component);
   return essTrueDofList;
  }
}
