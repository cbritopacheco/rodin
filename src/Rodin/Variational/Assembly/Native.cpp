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
    FormLanguage::List<BilinearFormIntegratorBase> facesBFIs;
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
        case Integrator::Region::Faces:
        {
          facesBFIs.add(bfi);
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
        const auto& element = *it;
        const Geometry::Attribute attr = element.getAttribute();
        for (const auto& bfi : domainBFIs)
        {
          if (bfi.getAttributes().size() == 0 || bfi.getAttributes().count(attr))
          {
            res.AddSubMatrix(
                input.testFES.getDOFs(element), input.trialFES.getDOFs(element),
                bfi.getMatrix(element));
          }
        }
      }
    }

    if (facesBFIs.size() > 0 || boundaryBFIs.size() > 0 || interfaceBFIs.size() > 0)
    {
      for (auto it = input.mesh.getFace(); !it.end(); ++it)
      {
        const auto& face = *it;
        const Geometry::Attribute attr = face.getAttribute();
        for (const auto& bfi : facesBFIs)
        {
          if (bfi.getAttributes().size() == 0 || bfi.getAttributes().count(attr))
          {
            res.AddSubMatrix(
                input.testFES.getDOFs(face), input.trialFES.getDOFs(face),
                bfi.getMatrix(face));
          }
        }

        if (face.isBoundary())
        {
          for (const auto& bfi : boundaryBFIs)
          {
            const Geometry::Attribute attr = input.mesh.getFaceAttribute(it->getIndex());
            if (bfi.getAttributes().size() == 0 || bfi.getAttributes().count(attr))
            {
              res.AddSubMatrix(
                  input.testFES.getDOFs(face), input.trialFES.getDOFs(face),
                  bfi.getMatrix(face));
            }
          }
        }

        if (face.isInterface())
        {
          for (const auto& bfi : interfaceBFIs)
          {
            const Geometry::Attribute attr = input.mesh.getFaceAttribute(it->getIndex());
            if (bfi.getAttributes().size() == 0 || bfi.getAttributes().count(attr))
            {
              res.AddSubMatrix(
                  input.testFES.getDOFs(face), input.trialFES.getDOFs(face),
                  bfi.getMatrix(face));
            }
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
    FormLanguage::List<LinearFormIntegratorBase> facesLFIs;
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
        case Integrator::Region::Faces:
        {
          facesLFIs.add(lfi);
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
        const auto& element = *it;
        const Geometry::Attribute attr = element.getAttribute();
        for (const auto& lfi : domainLFIs)
        {
          if (lfi.getAttributes().size() == 0 || lfi.getAttributes().count(attr))
          {
            res.AddElementVector(
              input.fes.getDOFs(element), lfi.getVector(element));
          }
        }
      }
    }

    if (facesLFIs.size() > 0 || boundaryLFIs.size() > 0 || interfaceLFIs.size() > 0)
    {
      for (auto it = input.mesh.getFace(); !it.end(); ++it)
      {
        const auto& face = *it;
        const Geometry::Attribute attr = face.getAttribute();
        for (const auto& bfi : facesLFIs)
        {
          if (bfi.getAttributes().size() == 0 || bfi.getAttributes().count(attr))
          {
            res.AddElementVector(input.fes.getDOFs(face), bfi.getVector(face));
          }
        }

        if (face.isBoundary())
        {
          for (const auto& lfi : boundaryLFIs)
          {
            if (lfi.getAttributes().size() == 0 || lfi.getAttributes().count(attr))
            {
              res.AddElementVector(input.fes.getDOFs(face), lfi.getVector(face));
            }
          }
        }

        if (face.isInterface())
        {
          for (const auto& lfi : interfaceLFIs)
          {
            const Geometry::Attribute attr = input.mesh.getFaceAttribute(it->getIndex());
            if (lfi.getAttributes().size() == 0 || lfi.getAttributes().count(attr))
            {
              res.AddElementVector(input.fes.getDOFs(face), lfi.getVector(face));
            }
          }
        }
      }
    }
    return res;
  }
}

