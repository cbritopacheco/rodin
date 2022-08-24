/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "GridFunction.h"

#include "Rodin/IO/GridFunctionLoader.h"

#include "Exceptions.h"

#include "Sum.h"
#include "Mult.h"
#include "Minus.h"
#include "Division.h"
#include "Component.h"
#include "BooleanFunction.h"
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

   GridFunctionBase& GridFunctionBase::operator+=(const GridFunctionBase& rhs)
   {
      if (this == &rhs)
      {
         operator*=(2.0);
      }
      else
      {
         assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
         getHandle() += rhs.getHandle();
      }
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator-=(double t)
   {
      getHandle() -= t;
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator-=(const GridFunctionBase& rhs)
   {
      if (this == &rhs)
      {
         operator=(0.0);
      }
      else
      {
         assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
         getHandle() -= rhs.getHandle();
      }
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator*=(double t)
   {
      getHandle() *= t;
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator*=(const GridFunctionBase& rhs)
   {
      if (this == &rhs)
      {
         for (auto& v : getHandle())
            v *= v;
      }
      else
      {
         assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
         getHandle() *= rhs.getHandle();
      }
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator/=(double t)
   {
      getHandle() /= t;
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator/=(const GridFunctionBase& rhs)
   {
      if (this == &rhs)
      {
         operator=(1.0);
      }
      else
      {
         assert(&getFiniteElementSpace() == &rhs.getFiniteElementSpace());
         getHandle() /= rhs.getHandle();
      }
      return *this;
   }

   GridFunctionBase& GridFunctionBase::operator=(double v)
   {
      getHandle() = v;
      return *this;
   }

   GridFunctionBase& GridFunctionBase::project(
         const FunctionBase& s, const std::set<int>& attrs)
   {
      auto va = s.build();
      std::visit(Utility::Overloaded{
         [&](Internal::ScalarProxyFunction& iv)
         {
            if (attrs.size() == 0)
               getHandle().ProjectCoefficient(iv);
            else
            {
               mfem::Array<int> vdofs;
               const auto& fes = getFiniteElementSpace().getHandle();
               for (int i = 0; i < fes.GetNE(); i++)
               {
                  if (attrs.count(fes.GetAttribute(i)) > 0)
                  {
                     fes.GetElementVDofs(i, vdofs);
                     getHandle().ProjectCoefficient(iv, vdofs);
                  }
               }
            }
         },
         [&](Internal::VectorProxyFunction& iv)
         {
            if (attrs.size() == 0)
               getHandle().ProjectCoefficient(iv);
            else
            {
               mfem::Array<int> vdofs;
               const auto& fes = getFiniteElementSpace().getHandle();
               for (int i = 0; i < fes.GetNE(); i++)
               {
                  if (attrs.count(fes.GetAttribute(i)) > 0)
                  {
                     fes.GetElementVDofs(i, vdofs);
                     getHandle().ProjectCoefficient(iv, vdofs);
                  }
               }
            }
         },
         [&](Internal::MatrixProxyFunction& iv)
         {
            UnexpectedRangeTypeException(
                  {RangeType::Scalar, RangeType::Vector}, RangeType::Matrix).raise();
         }}, va);
      return *this;
   }

   GridFunctionBase& GridFunctionBase::project(const Restriction<FunctionBase>& s)
   {
      assert(false); // Not implemented
      // assert(getFiniteElementSpace().getVectorDimension() == 1);
      // auto iv = s.getScalarFunction().build();
      // getHandle() = std::numeric_limits<double>::quiet_NaN();
      // mfem::Array<int> vdofs;
      // mfem::Vector vals;
      // const auto& fes = getFiniteElementSpace().getHandle();
      // const auto& attrs = s.getAttributes();
      // for (int i = 0; i < fes.GetNE(); i++)
      // {
      //    if (attrs.count(fes.GetAttribute(i)) > 0)
      //    {
      //       fes.GetElementVDofs(i, vdofs);
      //       vals.SetSize(vdofs.Size());
      //       fes.GetFE(i)->Project(
      //             *iv, *fes.GetElementTransformation(i), vals);
      //       getHandle().SetSubVector(vdofs, vals);
      //    }
      // }
      return *this;
   }

   void GridFunctionBase::transfer(GridFunctionBase& dst)
   {
      assert(getFiniteElementSpace().getVectorDimension() ==
            dst.getFiniteElementSpace().getVectorDimension());
      if (getFiniteElementSpace().getMesh().isSubMesh() && (
            &dst.getFiniteElementSpace().getMesh()) ==
            &static_cast<const SubMesh<Context::Serial>&>(
               getFiniteElementSpace().getMesh()).getParent())
      {
         // If we are here the this means that we are in a submesh of the
         // underlying target finite element space. Hence we should seek
         // out to copy the grid function at the corresponding nodes
         // given by the vertex map given in the Submesh object.
         auto& submesh = static_cast<const SubMesh<Context::Serial>&>(
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
            &static_cast<const SubMesh<Context::Serial>&>(
               dst.getFiniteElementSpace().getMesh()).getParent()))
      {
         auto& submesh = static_cast<const SubMesh<Context::Serial>&>(
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

   void GridFunctionBase::getValue(
         mfem::Vector& value,
         mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip) const
   {
      switch (getRangeType())
      {
         case RangeType::Scalar:
         {
            value.SetSize(1);
            value(0) = getHandle().GetValue(trans, ip);
            break;
         }
         case RangeType::Vector:
         {
            getHandle().GetVectorValue(trans, ip, value);
            break;
         }
         case RangeType::Matrix:
         {
            assert(false);
         }
      }
   }

   RangeShape GridFunctionBase::getRangeShape() const
   {
      return {getFiniteElementSpace().getVectorDimension(), 1};
   }

   Component<FunctionBase> GridFunctionBase::x() const
   {
      assert(getFiniteElementSpace().getVectorDimension() >= 1);
      return Component(*this, 0);
   }

   Component<FunctionBase> GridFunctionBase::y() const
   {
      assert(getFiniteElementSpace().getVectorDimension() >= 2);
      return Component(*this, 1);
   }

   Component<FunctionBase> GridFunctionBase::z() const
   {
      assert(getFiniteElementSpace().getVectorDimension() >= 3);
      return Component(*this, 2);
   }

   template <>
   GridFunction<H1<Context::Serial>>&
   GridFunction<H1<Context::Serial>>::load(
         const boost::filesystem::path& filename, IO::FileFormat fmt)
   {
      mfem::named_ifgzstream input(filename.c_str());
      if (!input)
      {
         Alert::Exception()
            << "Failed to open " << filename << " for reading."
            << Alert::Raise;
      }

      switch (fmt)
      {
         case IO::FileFormat::MFEM:
         {
            IO::GridFunctionLoader<IO::FileFormat::MFEM, H1<Context::Serial>> loader(*this);
            loader.load(input);
            break;
         }
         case IO::FileFormat::MEDIT:
         {
            IO::GridFunctionLoader<IO::FileFormat::MEDIT, H1<Context::Serial>> loader(*this);
            loader.load(input);
            break;
         }
         default:
         {
            Alert::Exception()
               << "Loading from \"" << fmt << "\" format unssuported."
               << Alert::Raise;
         }
      }
      return *this;
   }

   template <>
   void GridFunction<H1<Context::Serial>>
   ::save(const boost::filesystem::path& filename, IO::FileFormat fmt, int precision)
   const
   {
      std::ofstream output(filename.c_str());
      if (!output)
      {
         Alert::Exception()
            << "Failed to open " << filename << " for writing."
            << Alert::Raise;
      }

      output.precision(precision);
      switch (fmt)
      {
         case IO::FileFormat::MFEM:
         {
            IO::GridFunctionPrinter<IO::FileFormat::MFEM, H1<Context::Serial>> printer(*this);
            printer.print(output);
            break;
         }
         case IO::FileFormat::MEDIT:
         {
            IO::GridFunctionPrinter<IO::FileFormat::MEDIT, H1<Context::Serial>> printer(*this);
            printer.print(output);
            break;
         }
         default:
         {
            Alert::Exception()
               << "Saving to \"" << fmt << "\" format unssuported."
               << Alert::Raise;
         }
      }
   }

   std::set<Vertex> GridFunctionBase::where(
         const BooleanFunctionBase& p,
         const std::set<int>& attrs,
         std::function<int(mfem::ElementTransformation&)> order) const
   {
      std::set<Vertex> result;
      const auto& fes = getFiniteElementSpace();
      const auto& mesh = fes.getMesh();
      for (int i = 0; i < mesh.count<Element>(); i++)
      {
         if (attrs.size() == 0 || attrs.count(mesh.get<Element>(i).getAttribute()))
         {
            mfem::ElementTransformation* trans =
               fes.getHandle().GetElementTransformation(i);

            const mfem::IntegrationRule* ir =
               &mfem::IntRules.Get(trans->GetGeometryType(), order(*trans));

            for (int j = 0; j < ir->GetNPoints(); j++)
            {
               const mfem::IntegrationPoint& ip = ir->IntPoint(j);
               trans->SetIntPoint(&ip);
               if (p.getValue(*trans, ip))
               {
                  mfem::Vector v;
                  trans->Transform(ip, v);
                  Vertex vx(std::move(v));
                  vx.setElementTransformation(trans);
                  vx.setIntegrationPoint(&ip);
                  result.insert(std::move(vx));
               }
            }
         }
      }
      return result;
   }

   int GridFunctionBase::getDimension() const
   {
      return getFiniteElementSpace().getVectorDimension();
   }


   RangeType GridFunctionBase::getRangeType() const
   {
      auto shape = getRangeShape();
      if (shape.height() == 1 && shape.width() == 1)
      {
         return RangeType::Scalar;
      }
      else if (shape.height() > 1 && shape.width() == 1)
      {
         return RangeType::Vector;
      }
      else
      {
         assert(false);
         return RangeType::Matrix;
      }
   }

   VectorFunctionBase* GridFunctionBase::copy() const noexcept
   {
      return new Internal::GridFunctionEvaluator(*this);
   }
}

namespace Rodin::Variational::Internal
{
   GridFunctionEvaluator::GridFunctionEvaluator(const GridFunctionBase& gf)
      : m_gf(gf)
   {}

   GridFunctionEvaluator::GridFunctionEvaluator(const GridFunctionEvaluator& other)
      : VectorFunctionBase(other),
        m_gf(other.m_gf)
   {}

   GridFunctionEvaluator::GridFunctionEvaluator(GridFunctionEvaluator&& other)
      : VectorFunctionBase(std::move(other)),
        m_gf(other.m_gf)
   {}

   RangeShape GridFunctionEvaluator::getRangeShape() const
   {
      return m_gf.getRangeShape();
   }

   void GridFunctionEvaluator::getValue(
         mfem::Vector& value,
         mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip) const
   {
      m_gf.getValue(value, trans, ip);
   }
}
