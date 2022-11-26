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
   double BilinearForm<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix>
   ::operator()(const GridFunction<TrialFES>& u, const GridFunction<TestFES>& v) const
   {
      assert(m_operator);
      return m_operator->InnerProduct(u.getHandle(), v.getHandle());
   }

   template <class TrialFES, class TestFES>
   void
   BilinearForm<TrialFES, TestFES, Context::Serial, mfem::SparseMatrix>::assemble()
   {
      const auto& trialFes = getTrialFunction().getFiniteElementSpace();
      const auto& testFes = getTestFunction().getFiniteElementSpace();

      assert(&getTrialFunction().getFiniteElementSpace().getMesh() ==
            &getTestFunction().getFiniteElementSpace().getMesh());

      const auto& mesh = getTrialFunction().getFiniteElementSpace().getMesh();

      m_operator.reset(new OperatorType(testFes.getSize(), trialFes.getSize()));
      *m_operator = 0.0;

      for (const auto& bfi : getIntegrators())
      {
         switch (bfi.getRegion())
         {
            case Geometry::Region::Domain:
            {
               for (int i = 0; i < mesh.template count<Geometry::Element>(); i++)
               {
                  const auto& element = mesh.template get<Geometry::Element>(i);
                  if (bfi.getAttributes().size() == 0
                        || bfi.getAttributes().count(element.getAttribute()))
                  {
                     m_operator->AddSubMatrix(
                           testFes.getDOFs(element),
                           trialFes.getDOFs(element),
                           bfi.getElementMatrix(element));
                  }
               }
               break;
            }
            case Geometry::Region::Boundary:
            {
               assert(false);
               // mat.AddSubMatrix(
               //       testFes.getDOFs(element),
               //       testFes.getDOFs(element),
               //       bfi.getElementMatrix(element), true);
               break;
            }
            case Geometry::Region::Interface:
            {
               assert(false);
               break;
            }
         }
      }
   }
}

#endif
