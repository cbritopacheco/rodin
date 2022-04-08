#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
   int FiniteElementSpaceBase::getNumberOfDofs() const
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
