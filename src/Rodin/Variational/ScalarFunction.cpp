/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "ScalarFunction.h"

#include "Rodin/Mesh/SubMesh.h"

#include "Utility.h"
#include "Restriction.h"
#include "GridFunction.h"

namespace Rodin::Variational
{
   std::unique_ptr<mfem::Coefficient> ScalarFunctionBase::build() const
   {
      return std::unique_ptr<mfem::Coefficient>(new Internal::ProxyScalarFunction(*this));
   }

   Restriction<ScalarFunctionBase> ScalarFunctionBase::restrictTo(int attr)
   {
      return restrictTo(std::set<int>{attr});
   }

   Restriction<ScalarFunctionBase> ScalarFunctionBase::restrictTo(
         const std::set<int>& attrs)
   {
      return Restriction<ScalarFunctionBase>(*this).to(attrs);
   }

   // ---- ScalarFunction<GridFunctionBase> ----------------------------------
   ScalarFunction<GridFunctionBase>::ScalarFunction(const GridFunctionBase& u)
      : m_u(u)
   {
      assert(u.getFiniteElementSpace().getVectorDimension() == 1);
   }

   double ScalarFunction<GridFunctionBase>::getValue(
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const
   {
      mfem::Mesh* gfMesh = m_u.getHandle().FESpace()->GetMesh();
      if (trans.mesh == gfMesh)
      {
         return m_u.getHandle().GetValue(trans, ip);
      }
      else
      {
         mfem::IntegrationPoint coarseIp;
         mfem::ElementTransformation* coarseTrans = refinedToCoarse(
               *gfMesh, trans, ip, coarseIp);
         return m_u.getHandle().GetValue(*coarseTrans, coarseIp);
      }
   }

}
