/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "Serial.h"

namespace Rodin::Assembly
{
  void
  Serial<Variational::LinearFormBase<Math::Vector>>
  ::add(Math::Vector& out, const Math::Vector& in, const IndexArray& s)
  {
    assert(in.size() == s.size());
    size_t i = 0;
    for (const auto& global : s)
      out.coeffRef(global) += in.coeff(i++);
  }

  Math::Vector
  Serial<Variational::LinearFormBase<Math::Vector>>
  ::execute(const Input& input) const
  {
    VectorType res(input.fes.getSize());
    res.setZero();

    for (auto& lfi : input.lfis)
    {
      const auto& attrs = lfi.getAttributes();
      switch (lfi.getRegion())
      {
        case Variational::Integrator::Region::Domain:
        {
          for (auto it = input.mesh.getCell(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.fes.getDOFs(d, i);
              lfi.assemble(*it);
              add(res, lfi.getVector(), dofs);
            }
          }
          break;
        }
        case Variational::Integrator::Region::Faces:
        {
          for (auto it = input.mesh.getFace(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.fes.getDOFs(d, i);
              lfi.assemble(*it);
              add(res, lfi.getVector(), dofs);
            }
          }
          break;
        }
        case Variational::Integrator::Region::Boundary:
        {
          for (auto it = input.mesh.getBoundary(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.fes.getDOFs(d, i);
              lfi.assemble(*it);
              add(res, lfi.getVector(), dofs);
            }
          }
          break;
        }
        case Variational::Integrator::Region::Interface:
        {
          for (auto it = input.mesh.getInterface(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.fes.getDOFs(d, i);
              lfi.assemble(*it);
              add(res, lfi.getVector(), dofs);
            }
          }
          break;
        }
      }
    }

    return res;
  }
}


