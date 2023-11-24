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
  Math::SparseMatrix
  Serial<Variational::BilinearFormBase<Math::SparseMatrix>>
  ::execute(const BilinearAssemblyInput& input) const
  {
    Serial<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>> assembly;
    const auto triplets = assembly.execute(input);
    OperatorType res(input.testFES.getSize(), input.trialFES.getSize());
    res.setFromTriplets(triplets.begin(), triplets.end());
    return res;
  }

  void
  Serial<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
  ::add(std::vector<Eigen::Triplet<Scalar>>& out, const Math::Matrix& in,
      const IndexArray& rows, const IndexArray& cols)
  {
    assert(rows.size() >= 0);
    assert(cols.size() >= 0);
    assert(in.rows() == rows.size());
    assert(in.cols() == cols.size());
    for (size_t i = 0; i < static_cast<size_t>(rows.size()); i++)
    {
      for (size_t j = 0; j < static_cast<size_t>(cols.size()); j++)
      {
        const Scalar s = in(i, j);
        if (s != Scalar(0))
          out.emplace_back(rows(i), cols(j), s);
      }
    }
  }

  std::vector<Eigen::Triplet<Scalar>>
  Serial<Variational::BilinearFormBase<std::vector<Eigen::Triplet<Scalar>>>>
  ::execute(const BilinearAssemblyInput& input) const
  {
    std::vector<Eigen::Triplet<Scalar>> res;
    res.reserve(input.testFES.getSize() * std::log(input.trialFES.getSize()));
    for (auto& bfi : input.bfis)
    {
      const auto& attrs = bfi.getAttributes();
      switch (bfi.getRegion())
      {
        case Variational::Integrator::Region::Domain:
        {
          for (auto it = input.mesh.getCell(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              bfi.assemble(*it);
              add(res, bfi.getMatrix(), testDOFs, trialDOFs);
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
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              bfi.assemble(*it);
              add(res, bfi.getMatrix(), testDOFs, trialDOFs);
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
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              bfi.assemble(*it);
              add(res, bfi.getMatrix(), testDOFs, trialDOFs);
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
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              bfi.assemble(*it);
              add(res, bfi.getMatrix(), testDOFs, trialDOFs);
            }
          }
          break;
        }
      }
    }
    return res;
  }

  void
  Serial<Variational::BilinearFormBase<Math::Matrix>>
  ::add(Math::Matrix& out, const Math::Matrix& in,
      const IndexArray& rows, const IndexArray& cols)
  {
    assert(rows.size() >= 0);
    assert(cols.size() >= 0);
    assert(in.rows() == rows.size());
    assert(in.cols() == cols.size());
    for (size_t i = 0; i < static_cast<size_t>(rows.size()); i++)
    {
      for (size_t j = 0; j < static_cast<size_t>(cols.size()); j++)
        out(rows(i), cols(j)) += in(i, j);
    }
  }

  Math::Matrix
  Serial<Variational::BilinearFormBase<Math::Matrix>>
  ::execute(const BilinearAssemblyInput& input) const
  {
    Math::Matrix res(input.testFES.getSize(), input.trialFES.getSize());;
    res.setZero();
    for (auto& bfi : input.bfis)
    {
      const auto& attrs = bfi.getAttributes();
      switch (bfi.getRegion())
      {
        case Variational::Integrator::Region::Domain:
        {
          for (auto it = input.mesh.getCell(); !it.end(); ++it)
          {
            if (attrs.size() == 0 || attrs.count(it->getAttribute()))
            {
              const size_t d = it->getDimension();
              const size_t i = it->getIndex();
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              bfi.assemble(*it);
              add(res, bfi.getMatrix(), testDOFs, trialDOFs);
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
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              bfi.assemble(*it);
              add(res, bfi.getMatrix(), testDOFs, trialDOFs);
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
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              bfi.assemble(*it);
              add(res, bfi.getMatrix(), testDOFs, trialDOFs);
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
              const auto& trialDOFs = input.trialFES.getDOFs(d, i);
              const auto& testDOFs = input.testFES.getDOFs(d, i);
              bfi.assemble(*it);
              add(res, bfi.getMatrix(), testDOFs, trialDOFs);
            }
          }
          break;
        }
      }
    }
    return res;
  }

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


