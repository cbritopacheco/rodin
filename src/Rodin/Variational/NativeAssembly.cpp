#include "NativeAssembly.h"

namespace Rodin::Variational
{
   mfem::SparseMatrix NativeAssembly<BilinearFormBase<mfem::SparseMatrix>>
   ::execute(const Input& input) const
   {
      OperatorType res(input.testFES.getSize(), input.trialFES.getSize());
      res = 0.0;

      for (const auto& bfi : input.bfis)
      {
         switch (bfi.getRegion())
         {
            case Geometry::Region::Domain:
            {
               for (int i = 0; i < input.mesh.template count<Geometry::Element>(); i++)
               {
                  const auto& element = input.mesh.template get<Geometry::Element>(i);
                  if (bfi.getAttributes().size() == 0
                        || bfi.getAttributes().count(element.getAttribute()))
                  {
                     res.AddSubMatrix(
                           input.testFES.getDOFs(element),
                           input.trialFES.getDOFs(element),
                           bfi.getElementMatrix(element));
                  }
               }
               break;
            }
            case Geometry::Region::Boundary:
            {
               assert(false);
               // mat.AddSubMatrix(
               //       testFes.getDOFs(element),
               //       testFes.getDOFs(element),
               //       bfi.getElementMatrix(element), true);
               break;
            }
            case Geometry::Region::Interface:
            {
               assert(false);
               break;
            }
         }
      }

      return res;
   }
}
