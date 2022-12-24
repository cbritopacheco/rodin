#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
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
               for (auto it = input.mesh.getElement(); !it.end(); ++it)
               {
                  const auto& element = *it;
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
               for (auto it = input.mesh.getBoundary(); !it.end(); ++it)
               {
                  const auto& element = *it;
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
               for (auto it = input.mesh.getInterface(); !it.end(); ++it)
               {
                  const auto& element = *it;
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
         }
      }
      return res;
   }

   mfem::Vector
   Native<LinearFormBase<mfem::Vector>>
   ::execute(const Input& input) const
   {
      VectorType res(input.fes.getSize());
      res = 0.0;

      for (const auto& lfi : input.lfis)
      {
         switch (lfi.getRegion())
         {
            case Geometry::Region::Domain:
            {
               for (auto it = input.mesh.getElement(); !it.end(); ++it)
               {
                  const auto& element = *it;
                  if (lfi.getAttributes().size() == 0
                        || lfi.getAttributes().count(element.getAttribute()))
                  {
                     res.AddElementVector(input.fes.getDOFs(element), lfi.getElementVector(element));
                  }
               }
               break;
            }
            case Geometry::Region::Boundary:
            {
               for (auto it = input.mesh.getBoundary(); !it.end(); ++it)
               {
                  const auto& element = *it;
                  std::cout << "miaow: " << element.getAttribute() << std::endl;
                  if (lfi.getAttributes().size() == 0
                        || lfi.getAttributes().count(element.getAttribute()))
                  {
                     res.AddElementVector(input.fes.getDOFs(element), lfi.getElementVector(element));
                  }
               }
               break;
            }
            case Geometry::Region::Interface:
            {
               for (auto it = input.mesh.getInterface(); !it.end(); ++it)
               {
                  const auto& element = *it;
                  if (lfi.getAttributes().size() == 0
                        || lfi.getAttributes().count(element.getAttribute()))
                  {
                     res.AddElementVector(input.fes.getDOFs(element), lfi.getElementVector(element));
                  }
               }
               break;
            }
         }
      }

      return res;
   }
}


