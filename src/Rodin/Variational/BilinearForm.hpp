/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_HPP
#define RODIN_VARIATIONAL_BILINEARFORM_HPP

#include <cassert>

#include "Rodin/Alert.h"

#include "FiniteElementSpace.h"
#include "BilinearFormIntegrator.h"

#include "BilinearForm.h"

namespace Rodin::Variational
{
   template <class TrialFES, class TestFES>
   constexpr
   double
   BilinearForm<TrialFES, TestFES, Context::Serial, mfem::Operator>
   ::operator()(const GridFunction<TrialFES>& u, const GridFunction<TestFES>& v) const
   {
      assert(m_bf);
      return m_bf->InnerProduct(u.getHandle(), v.getHandle());
   }

   template <class TrialFES, class TestFES>
   BilinearForm<TrialFES, TestFES, Context::Serial, mfem::Operator>&
   BilinearForm<TrialFES, TestFES, Context::Serial, mfem::Operator>::update()
   {
      assert(m_bf);
      m_bf->Update();
      return *this;
   }

   template <class TrialFES, class TestFES>
   void
   BilinearForm<TrialFES, TestFES, Context::Serial, mfem::Operator>::assemble()
   {
      m_bf.reset(new mfem::BilinearForm(&m_v.getFiniteElementSpace().getHandle()));

      std::vector<mfem::Array<int>> attrs;
      attrs.reserve(getIntegrators().size());
      for (const auto& bfi : getIntegrators())
      {
         switch (bfi.getIntegratorRegion())
         {
            case IntegratorRegion::Boundary:
            {
               if (bfi.getAttributes().size() == 0)
               {
                  m_bf->AddBoundaryIntegrator(bfi.build().release());
               }
               else
               {
                  const int size =
                     m_v.getFiniteElementSpace().getMesh().getHandle().bdr_attributes.Max();
                  m_bf->AddBoundaryIntegrator(
                        bfi.build().release(),
                        attrs.emplace_back(Utility::set2marker(bfi.getAttributes(), size)));
               }
               break;
            }
            case IntegratorRegion::Domain:
            {
               if (bfi.getAttributes().size() == 0)
               {
                  m_bf->AddDomainIntegrator(bfi.build().release());
               }
               else
               {
                  const int size =
                     m_v.getFiniteElementSpace().getMesh().getHandle().attributes.Max();
                  m_bf->AddDomainIntegrator(
                        bfi.build().release(),
                        attrs.emplace_back(Utility::set2marker(bfi.getAttributes(), size)));
               }
               break;
            }
         }
      }

      m_bf->Assemble();
   }
}

#endif
