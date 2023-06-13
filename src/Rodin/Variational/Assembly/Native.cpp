/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "Native.h"

namespace Rodin::Variational::Assembly
{
  void Native<BilinearFormBase<Math::SparseMatrix>>
  ::add(Math::SparseMatrix& out, const Math::Matrix& in,
      const IndexSet& rows, const IndexSet& cols)
  {
    assert(static_cast<size_t>(in.rows()) == rows.size());
    assert(static_cast<size_t>(in.cols()) == cols.size());
    size_t i = 0;
    for (const auto& r : rows)
    {
      size_t j = 0;
      for (const auto& c : cols)
        out.coeffRef(r, c) = in.coeff(i, j++);
      i++;
    }
  }

  Math::SparseMatrix Native<BilinearFormBase<Math::SparseMatrix>>
  ::execute(const Input& input) const
  {
    OperatorType res(input.testFES.getSize(), input.trialFES.getSize());
    res.setZero();

    for (const auto& bfi : input.bfis)
    {
      const auto& attrs = bfi.getAttributes();
      switch (bfi.getRegion())
      {
        case Integrator::Region::Domain:
        {
          for (auto it = input.mesh.getElement(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              const auto mat = bfi.getMatrix(*it);
              add(res, mat, testDOFs, trialDOFs);
            }
          }
          break;
        }
        case Integrator::Region::Faces:
        {
          for (auto it = input.mesh.getFace(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              const auto mat = bfi.getMatrix(*it);
              add(res, mat, testDOFs, trialDOFs);
            }
          }
          break;
        }
        case Integrator::Region::Boundary:
        {
          for (auto it = input.mesh.getBoundary(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              const auto mat = bfi.getMatrix(*it);
              add(res, mat, testDOFs, trialDOFs);
            }
          }
          break;
        }
        case Integrator::Region::Interface:
        {
          for (auto it = input.mesh.getInterface(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              const auto mat = bfi.getMatrix(*it);
              add(res, mat, testDOFs, trialDOFs);
            }
          }
          break;
        }
      }
    }
    return res;
  }

  void Native<LinearFormBase<Math::Vector>>
  ::add(Math::Vector& out, const Math::Vector& in, const IndexSet& s)
  {
    assert(static_cast<size_t>(in.size()) == s.size());
    size_t i = 0;
    for (const auto& global : s)
      out.coeffRef(global) += in.coeff(i++);
  }

  Math::Vector Native<LinearFormBase<Math::Vector>>
  ::execute(const Input& input) const
  {
    VectorType res(input.fes.getSize());
    res.setZero();

    for (const auto& lfi : input.lfis)
    {
      const auto& attrs = lfi.getAttributes();
      switch (lfi.getRegion())
      {
        case Integrator::Region::Domain:
        {
          for (auto it = input.mesh.getElement(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.fes.getDOFs(d, i);
              const auto vec = lfi.getVector(*it);
              add(res, vec, dofs);
            }
          }
          break;
        }
        case Integrator::Region::Faces:
        {
          for (auto it = input.mesh.getFace(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.fes.getDOFs(d, i);
              const auto vec = lfi.getVector(*it);
              add(res, vec, dofs);
            }
          }
          break;
        }
        case Integrator::Region::Boundary:
        {
          for (auto it = input.mesh.getBoundary(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.fes.getDOFs(d, i);
              const auto vec = lfi.getVector(*it);
              add(res, vec, dofs);
            }
          }
          break;
        }
        case Integrator::Region::Interface:
        {
          for (auto it = input.mesh.getInterface(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& dofs = input.fes.getDOFs(d, i);
              const auto vec = lfi.getVector(*it);
              add(res, vec, dofs);
            }
          }
          break;
        }
      }
    }

    return res;
  }
}


