/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "GridFunction.h"

#include "Rodin/IO/GridFunctionLoader.h"

#include "Sum.h"
#include "Mult.h"
#include "Minus.h"
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

   void GridFunctionBase::transfer(GridFunctionBase& dst)
   {
      assert(getFiniteElementSpace().getVectorDimension() ==
            dst.getFiniteElementSpace().getVectorDimension());
      if (getFiniteElementSpace().getMesh().isSubMesh() && (
            &dst.getFiniteElementSpace().getMesh()) ==
            &static_cast<const SubMesh<Traits::Serial>&>(
               getFiniteElementSpace().getMesh()).getParent())
      {
         // If we are here the this means that we are in a submesh of the
         // underlying target finite element space. Hence we should seek
         // out to copy the grid function at the corresponding nodes
         // given by the vertex map given in the Submesh object.
         auto& submesh = static_cast<const SubMesh<Traits::Serial>&>(
               getFiniteElementSpace().getMesh());
         if (&submesh.getParent() == &dst.getFiniteElementSpace().getMesh())
         {
            int vdim = getFiniteElementSpace().getVectorDimension();
            const auto& s2pv = submesh.getVertexMap();
            int nv = getFiniteElementSpace().getHandle().GetNV();
            int pnv = dst.getFiniteElementSpace().getHandle().GetNV();

            assert(getFiniteElementSpace().getHandle().GetOrdering() ==
                     getFiniteElementSpace().getHandle().GetOrdering());
            switch(getFiniteElementSpace().getHandle().GetOrdering())
            {
               case mfem::Ordering::byNODES:
               {
                  for (int i = 0; i < vdim; i++)
                     for (int j = 0; j < nv; j++)
                        dst.getHandle()[s2pv.left.at(j) + i * pnv] = getHandle()[j + i * nv];
                  return;
               }
               case mfem::Ordering::byVDIM:
               {
                  for (int i = 0; i < nv; i++)
                     for (int j = 0; j < vdim; j++)
                        dst.getHandle()[s2pv.left.at(i) * vdim + j] = getHandle()[i * vdim + j];
                  return;
               }
            }
         }
      }
      else if (dst.getFiniteElementSpace().getMesh().isSubMesh() && (
            &getFiniteElementSpace().getMesh() ==
            &static_cast<const SubMesh<Traits::Serial>&>(
               dst.getFiniteElementSpace().getMesh()).getParent()))
      {
         auto& submesh = static_cast<const SubMesh<Traits::Serial>&>(
               dst.getFiniteElementSpace().getMesh());
         int vdim = getFiniteElementSpace().getVectorDimension();
         const auto& s2pv = submesh.getVertexMap();
         int nv = getFiniteElementSpace().getHandle().GetNV();
         int pnv = dst.getFiniteElementSpace().getHandle().GetNV();

         assert(getFiniteElementSpace().getHandle().GetOrdering() ==
                  getFiniteElementSpace().getHandle().GetOrdering());
         switch(getFiniteElementSpace().getHandle().GetOrdering())
         {
            case mfem::Ordering::byNODES:
            {
               for (int i = 0; i < vdim; i++)
                  for (int j = 0; j < nv; j++)
                     if (s2pv.right.count(j))
                        dst.getHandle()[s2pv.right.at(j) + i * pnv] = getHandle()[j + i * nv];
               return;
            }
            case mfem::Ordering::byVDIM:
            {
               for (int i = 0; i < nv; i++)
                  for (int j = 0; j < vdim; j++)
                     if (s2pv.right.count(i))
                        dst.getHandle()[s2pv.right.at(i) * vdim + j] = getHandle()[i * vdim + j];
               return;
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
