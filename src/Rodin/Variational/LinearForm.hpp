/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORM_HPP
#define RODIN_VARIATIONAL_LINEARFORM_HPP

#include "Rodin/Alert.h"
#include "FiniteElementSpace.h"

#include "LinearForm.h"

namespace Rodin::Variational
{
   template <class FES>
   constexpr
   LinearForm<FES, Context::Serial, mfem::Vector>
   ::LinearForm(TestFunction<FES>& v)
      :  m_v(v)
   {}

   template <class FES>
   constexpr
   double
   LinearForm<FES, Context::Serial, mfem::Vector>
   ::operator()(const GridFunction<FES>& u) const
   {
      assert(m_lf);
      return *m_lf * u.getHandle();
   }

   template <class FES>
   void
   LinearForm<FES, Context::Serial, mfem::Vector>::assemble()
   {
      m_lf.reset(new mfem::LinearForm(&m_v.getFiniteElementSpace().getHandle()));

      std::vector<mfem::Array<int>> attrs;
      attrs.reserve(getIntegrators().size());
      for (const auto& lfi : getIntegrators())
      {
         switch (lfi.getRegion())
         {
            case Geometry::Region::Boundary:
            {
               if (lfi.getAttributes().size() == 0)
               {
                  m_lf->AddBoundaryIntegrator(lfi.build().release());
               }
               else
               {
                  const int size =
                     m_v.getFiniteElementSpace().getMesh().getHandle().bdr_attributes.Max();
                  m_lf->AddBoundaryIntegrator(
                        lfi.build().release(),
                        attrs.emplace_back(Utility::set2marker(lfi.getAttributes(), size)));
               }
               break;
            }
            case Geometry::Region::Domain:
            {
               if (lfi.getAttributes().size() == 0)
               {
                  m_lf->AddDomainIntegrator(lfi.build().release());
               }
               else
               {
                  const int size =
                     m_v.getFiniteElementSpace().getMesh().getHandle().attributes.Max();
                  m_lf->AddDomainIntegrator(
                        lfi.build().release(),
                        attrs.emplace_back(Utility::set2marker(lfi.getAttributes(), size)));
               }
               break;
            }
         }
      }

      m_lf->Assemble();
   }
}

#endif
