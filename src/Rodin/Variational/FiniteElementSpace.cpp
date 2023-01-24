/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
  size_t FiniteElementSpaceBase::getNumberOfDofs() const
  {
    return getHandle().GetNDofs();
  }

  int FiniteElementSpaceBase::getVectorDimension() const
  {
    return getHandle().GetVDim();
  }

  void FiniteElementSpaceBase::update()
  {
    getHandle().Update();
  }

  mfem::Array<int> FiniteElementSpaceBase::getEssentialTrueDOFs(
      const std::set<int>& bdrAttr) const
  {
    mfem::Array<int> essTrueDofList;
    int maxBdrAttr = *getMesh().getBoundaryAttributes().rbegin();
    const_cast<mfem::FiniteElementSpace&>(getHandle()
        ).GetEssentialTrueDofs(
          Utility::set2marker(bdrAttr, maxBdrAttr), essTrueDofList);
    return essTrueDofList;
  }

  mfem::Array<int> FiniteElementSpaceBase::getEssentialTrueDOFs(
      const std::set<int>& bdrAttr, int component) const
  {
    mfem::Array<int> essTrueDofList;
    int maxBdrAttr = *getMesh().getBoundaryAttributes().rbegin();
    const_cast<mfem::FiniteElementSpace&>(getHandle()
        ).GetEssentialTrueDofs(
          Utility::set2marker(bdrAttr, maxBdrAttr), essTrueDofList, component);
    return essTrueDofList;
  }
}
