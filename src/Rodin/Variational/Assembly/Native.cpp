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

      FormLanguage::List<BilinearFormIntegratorBase> domainBFIs;
      FormLanguage::List<BilinearFormIntegratorBase> boundaryBFIs;
      FormLanguage::List<BilinearFormIntegratorBase> interfaceBFIs;

      for (const auto& bfi : input.bfis)
      {
         switch (bfi.getRegion())
         {
            case Integrator::Region::Domain:
            {
               domainBFIs.add(bfi);
               break;
            }
            case Integrator::Region::Boundary:
            {
               boundaryBFIs.add(bfi);
               break;
            }
            case Integrator::Region::Interface:
            {
               interfaceBFIs.add(bfi);
               break;
            }
         }
      }

      if (domainBFIs.size() > 0)
      {
         for (auto it = input.mesh.getElement(); !it.end(); ++it)
         {
            std::optional<Geometry::Element> cached;
            for (const auto& bfi : domainBFIs)
            {
               const Geometry::Attribute attr = input.mesh.getElementAttribute(it.getIndex());
               if (bfi.getAttributes().size() == 0 || bfi.getAttributes().count(attr))
               {
                  if (!cached.has_value())
                     cached.emplace(std::move(*it));
                  assert(cached.has_value());
                  const auto& element = cached.value();
                  res.AddSubMatrix(
                        input.testFES.getDOFs(element), input.trialFES.getDOFs(element),
                        bfi.getMatrix(element));
               }
            }
         }
      }

      if (boundaryBFIs.size() > 0)
      {
         for (auto it = input.mesh.getBoundary(); !it.end(); ++it)
         {
            std::optional<Geometry::Boundary> cached;
            for (const auto& bfi : boundaryBFIs)
            {
               const Geometry::Attribute attr = input.mesh.getFaceAttribute(it.getIndex());
               if (bfi.getAttributes().size() == 0 || bfi.getAttributes().count(attr))
               {
                  if (!cached.has_value())
                     cached.emplace(std::move(*it));
                  assert(cached.has_value());
                  const auto& boundary = cached.value();
                  res.AddSubMatrix(
                        input.testFES.getDOFs(boundary), input.trialFES.getDOFs(boundary),
                        bfi.getMatrix(boundary));
               }
            }
         }
      }

      if (interfaceBFIs.size() > 0)
      {
         for (auto it = input.mesh.getInterface(); !it.end(); ++it)
         {
            std::optional<Geometry::Interface> cached;
            for (const auto& bfi : interfaceBFIs)
            {
               const Geometry::Attribute attr = input.mesh.getFaceAttribute(it.getIndex());
               if (bfi.getAttributes().size() == 0 || bfi.getAttributes().count(attr))
               {
                  if (!cached.has_value())
                     cached.emplace(std::move(*it));
                  assert(cached.has_value());
                  const auto& interface = cached.value();
                  res.AddSubMatrix(
                        input.testFES.getDOFs(interface),
                        input.trialFES.getDOFs(interface),
                        bfi.getMatrix(interface));
               }
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

      FormLanguage::List<LinearFormIntegratorBase> domainLFIs;
      FormLanguage::List<LinearFormIntegratorBase> boundaryLFIs;
      FormLanguage::List<LinearFormIntegratorBase> interfaceLFIs;

      for (const auto& lfi : input.lfis)
      {
         switch (lfi.getRegion())
         {
            case Integrator::Region::Domain:
            {
               domainLFIs.add(lfi);
               break;
            }
            case Integrator::Region::Boundary:
            {
               boundaryLFIs.add(lfi);
               break;
            }
            case Integrator::Region::Interface:
            {
               interfaceLFIs.add(lfi);
               break;
            }
         }
      }

      if (domainLFIs.size() > 0)
      {
         for (auto it = input.mesh.getElement(); !it.end(); ++it)
         {
            std::optional<Geometry::Element> cached;
            for (const auto& lfi : domainLFIs)
            {
               const Geometry::Attribute attr = input.mesh.getElementAttribute(it.getIndex());
               if (lfi.getAttributes().size() == 0 || lfi.getAttributes().count(attr))
               {
                  if (!cached.has_value())
                     cached.emplace(std::move(*it));
                  assert(cached.has_value());
                  const auto& element = cached.value();
                  res.AddElementVector(
                     input.fes.getDOFs(element), lfi.getVector(element));
               }
            }
         }
      }

      if (boundaryLFIs.size() > 0)
      {
         for (auto it = input.mesh.getBoundary(); !it.end(); ++it)
         {
            std::optional<Geometry::Boundary> cached;
            for (const auto& lfi : boundaryLFIs)
            {
               const Geometry::Attribute attr = input.mesh.getFaceAttribute(it.getIndex());
               if (lfi.getAttributes().size() == 0 || lfi.getAttributes().count(attr))
               {
                  if (!cached.has_value())
                     cached.emplace(std::move(*it));
                  assert(cached.has_value());
                  const auto& boundary = cached.value();
                  res.AddElementVector(
                     input.fes.getDOFs(boundary), lfi.getVector(boundary));
               }
            }
         }
      }

      if (interfaceLFIs.size() > 0)
      {
         for (auto it = input.mesh.getInterface(); !it.end(); ++it)
         {
            std::optional<Geometry::Interface> cached;
            for (const auto& lfi : interfaceLFIs)
            {
               const Geometry::Attribute attr = input.mesh.getFaceAttribute(it.getIndex());
               if (lfi.getAttributes().size() == 0 || lfi.getAttributes().count(attr))
               {
                  if (!cached.has_value())
                     cached.emplace(std::move(*it));
                  assert(cached.has_value());
                  const auto& boundary = cached.value();
                  res.AddElementVector(
                     input.fes.getDOFs(boundary), lfi.getVector(boundary));
               }
            }
         }
      }
      return res;
   }
}


