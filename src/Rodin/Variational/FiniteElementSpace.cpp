#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
   int FiniteElementSpaceBase::getNumberOfDofs() const
   {
      return getFES().GetNDofs();
   }

   int FiniteElementSpaceBase::getVectorDimension() const
   {
      return getFES().GetVDim();
   }

   int FiniteElementSpaceBase::getOrder() const
   {
      return getFEC().GetOrder();
   }

   void FiniteElementSpaceBase::update()
   {
      getFES().Update();
   }

   mfem::Array<int> FiniteElementSpaceBase::getEssentialTrueDOFs(
         const std::set<int>& bdrAttr) const
   {
      mfem::Array<int> essTrueDofList;
      int maxBdrAttr = *getMesh().getBoundaryAttributes().rbegin();
      const_cast<mfem::FiniteElementSpace&>(getFES()
            ).GetEssentialTrueDofs(
               Utility::set2marker(bdrAttr, maxBdrAttr), essTrueDofList);
      return essTrueDofList;
   }

   mfem::Array<int> FiniteElementSpaceBase::getEssentialTrueDOFs(
         const std::set<int>& bdrAttr, int component) const
   {
      mfem::Array<int> essTrueDofList;
      int maxBdrAttr = *getMesh().getBoundaryAttributes().rbegin();
      const_cast<mfem::FiniteElementSpace&>(getFES()
            ).GetEssentialTrueDofs(
               Utility::set2marker(bdrAttr, maxBdrAttr), essTrueDofList, component);
      return essTrueDofList;
   }
}
