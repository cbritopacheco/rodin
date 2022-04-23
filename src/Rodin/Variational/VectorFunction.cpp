/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "VectorFunction.h"

#include "Component.h"
#include "GridFunction.h"

namespace Rodin::Variational
{
   Component<VectorFunctionBase> VectorFunctionBase::operator()(int i) const
   {
      assert(0 <= i);
      assert(i < getDimension());
      return Component<VectorFunctionBase>(*this, i);
   }

   Component<VectorFunctionBase> VectorFunctionBase::x() const
   {
      assert(getDimension() >= 1);
      return Component<VectorFunctionBase>(*this, 0);
   }

   Component<VectorFunctionBase> VectorFunctionBase::y() const
   {
      assert(getDimension() >= 2);
      return Component<VectorFunctionBase>(*this, 1);
   }

   Component<VectorFunctionBase> VectorFunctionBase::z() const
   {
      assert(getDimension() >= 3);
      return Component<VectorFunctionBase>(*this, 2);
   }

   std::unique_ptr<mfem::VectorCoefficient> VectorFunctionBase::build() const
   {
      return std::unique_ptr<mfem::VectorCoefficient>(new Internal::ProxyVectorFunction(*this));
   }

   // ---- VectorFunction<GridFunctionBase> ----------------------------------
   VectorFunction<GridFunctionBase>::VectorFunction(const GridFunctionBase& u)
      :  m_dimension(u.getFiniteElementSpace().getVectorDimension()),
         m_u(u)
   {}

   void VectorFunction<GridFunctionBase>::getValue(
         mfem::Vector& value,
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const
   {
      mfem::Mesh* gfMesh = m_u.getHandle().FESpace()->GetMesh();
      if (trans.mesh == gfMesh)
      {
         m_u.getHandle().GetVectorValue(trans, ip, value);
      }
      else
      {
         mfem::IntegrationPoint coarseIp;
         mfem::ElementTransformation* coarseTrans = refinedToCoarse(
               *gfMesh, trans, ip, coarseIp);
         m_u.getHandle().GetVectorValue(*coarseTrans, coarseIp, value);
      }
   }
}

