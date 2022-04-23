/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "GridFunction.h"

#include "Sum.h"
#include "UnaryMinus.h"
#include "Mult.h"
#include "Division.h"

#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
   GridFunctionBase& GridFunctionBase::update()
   {
      getHandle().Update();
      return *this;
   }

   double GridFunctionBase::max() const
   {
      return getHandle().Max();
   }

   double GridFunctionBase::min() const
   {
      return getHandle().Min();
   }

   void GridFunctionBase::save(const boost::filesystem::path& filename)
   {
      getHandle().Save(filename.string().c_str());
   }

   GridFunctionBase& GridFunctionBase::load(const boost::filesystem::path& filename)
   {
      std::ifstream in(filename.string());
      getHandle().Load(in, getHandle().Size());
      return *this;
   }

   std::pair<const double*, int> GridFunctionBase::getData() const
   {
      return {static_cast<const double*>(getHandle().GetData()), getHandle().Size()};
   }

   GridFunctionBase& GridFunctionBase::setData(std::unique_ptr<double[]> data, int size)
   {
      getHandle().SetDataAndSize(data.release(), size);
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator+=(double t)
   {
      getHandle() += t;
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator-=(double t)
   {
      getHandle() -= t;
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator*=(double t)
   {
      getHandle() *= t;
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator/=(double t)
   {
      getHandle() /= t;
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator=(double v)
   {
      getHandle() = v;
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator+=(const ScalarFunctionBase& v)
   {
      return project(ScalarFunction(*this) + v);
   }

   GridFunctionBase& GridFunctionBase::operator-=(const ScalarFunctionBase& v)
   {
      return project(ScalarFunction(*this) - v);
   }

   GridFunctionBase& GridFunctionBase::operator*=(const ScalarFunctionBase& v)
   {
      if (getFiniteElementSpace().getVectorDimension() > 1)
         return project(VectorFunction(*this) * v);
      else
         return project(ScalarFunction(*this) * v);
   }

   GridFunctionBase& GridFunctionBase::operator/=(const ScalarFunctionBase& v)
   {
      if (getFiniteElementSpace().getVectorDimension() > 1)
         return project(VectorFunction(*this) / v);
      else
         return project(ScalarFunction(*this) / v);
   }

   GridFunctionBase& GridFunctionBase::operator+=(const VectorFunctionBase& v)
   {
      return project(VectorFunction(*this) + v);
   }

   GridFunctionBase& GridFunctionBase::operator-=(const VectorFunctionBase& v)
   {
      return project(VectorFunction(*this) - v);
   }

   GridFunctionBase& GridFunctionBase::project(
         const ScalarFunctionBase& s, const std::set<int>& attrs)
   {
      assert(getFiniteElementSpace().getVectorDimension() == 1);
      auto iv = s.build();

      if (attrs.size() == 0)
         getHandle().ProjectCoefficient(*iv);
      else
      {
         int maxAttr = getFiniteElementSpace()
                      .getMesh()
                      .getHandle().attributes.Max();
         mfem::Array<int> marker(maxAttr);
         marker = 0;
         for (const auto& attr : attrs)
         {
            assert(attr - 1 < maxAttr);
            marker[attr - 1] = 1;
         }
         getHandle().ProjectCoefficient(*iv, marker);
      }
      return *this;
   }

   GridFunctionBase& GridFunctionBase::project(
         const VectorFunctionBase& s, const std::set<int>& attrs)
   {
      assert(getFiniteElementSpace().getVectorDimension() == s.getDimension());
      auto iv = s.build();

      if (attrs.size() == 0)
         getHandle().ProjectCoefficient(*iv);
      else
      {
         int maxAttr = getFiniteElementSpace()
                      .getMesh()
                      .getHandle().attributes.Max();
         mfem::Array<int> marker(maxAttr);
         marker = 0;
         for (const auto& attr : attrs)
         {
            assert(attr - 1 < maxAttr);
            marker[attr - 1] = 1;
         }
         getHandle().ProjectCoefficient(*iv, marker);
      }
      return *this;
   }

   GridFunctionBase& GridFunctionBase::project(const Restriction<ScalarFunctionBase>& s)
   {
      assert(getFiniteElementSpace().getVectorDimension() == 1);
      auto iv = s.getScalarFunction().build();
      getHandle() = std::numeric_limits<double>::quiet_NaN();
      mfem::Array<int> vdofs;
      mfem::Vector vals;
      const auto& fes = getFiniteElementSpace().getHandle();
      const auto& attrs = s.getAttributes();
      for (int i = 0; i < fes.GetNE(); i++)
      {
         if (attrs.count(fes.GetAttribute(i)) > 0)
         {
            fes.GetElementVDofs(i, vdofs);
            vals.SetSize(vdofs.Size());
            fes.GetFE(i)->Project(
                  *iv, *fes.GetElementTransformation(i), vals);
            getHandle().SetSubVector(vdofs, vals);
         }
      }
      return *this;
   }

   void GridFunctionBase::transfer(GridFunctionBase& other)
   {
      assert(getFiniteElementSpace().getVectorDimension() ==
            other.getFiniteElementSpace().getVectorDimension());
      if (getFiniteElementSpace().getMesh().isSubMesh())
      {
         // If we are here the this means that we are in a submesh of the
         // underlying target finite element space. Hence we should seek
         // out to copy the grid function at the corresponding nodes
         // given by the vertex map given in the Submesh object.
         auto& submesh = static_cast<const SubMesh<Traits::Serial>&>(
               getFiniteElementSpace().getMesh());
         if (&submesh.getParent() == &other.getFiniteElementSpace().getMesh())
         {
            int vdim = getFiniteElementSpace().getVectorDimension();
            const auto& s2pv = submesh.getVertexMap();
            if (vdim == 1)
            {
               int size = getHandle().Size();
               for (int i = 0; i < size; i++)
                  other.getHandle()[i] = getHandle()[s2pv.at(i)];
            }
            else
            {
               int nv = getFiniteElementSpace().getHandle().GetNV();
               int pnv = other.getFiniteElementSpace().getHandle().GetNV();

               assert(getFiniteElementSpace().getHandle().GetOrdering() ==
                        getFiniteElementSpace().getHandle().GetOrdering());
               switch(getFiniteElementSpace().getHandle().GetOrdering())
               {
                  case mfem::Ordering::byNODES:
                  {
                     for (int i = 0; i < vdim; i++)
                        for (int j = 0; j < nv; j++)
                           other.getHandle()[s2pv.at(j) + i * pnv] = getHandle()[j + i * nv];
                     return;
                  }
                  case mfem::Ordering::byVDIM:
                  {
                     for (int i = 0; i < nv; i++)
                        for (int j = 0; j < vdim; j++)
                           other.getHandle()[s2pv.at(i) * vdim + j] = getHandle()[i * vdim + j];
                     return;
                  }
               }
            }
         }
      }
      else
      {
         // If the meshes are equal or where obtained from refinements
         // one could use the mfem functionality to make a GridTransfer.
         // Alternatively, if the mesh is equal but the finite element
         // spaces are not, mfem also contains the TransferOperator class
         // which can come in useful.
         Alert::Exception("Unimplemented. Sorry.").raise();
      }
   }
}
