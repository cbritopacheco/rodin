#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "Native.h"

namespace Rodin::Variational::Assembly
{
   mfem::SparseMatrix
   Native<BilinearFormBase<mfem::SparseMatrix>>
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
               for (int i = 0; i < input.mesh.count<Geometry::Element>(); i++)
               {
                  const auto& element = input.mesh.get<Geometry::Element>(i);
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
               for (int i = 0; i < input.mesh.count<Geometry::Boundary>(); i++)
               {
                  const auto& element = input.mesh.get<Geometry::Boundary>(i);
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


