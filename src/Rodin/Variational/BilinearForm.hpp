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

#include "Assembly/AssemblyBase.h"

#include "BilinearForm.h"
#include "BilinearFormIntegrator.h"
#include "FiniteElementSpace.h"


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
      assert(&getTrialFunction().getFiniteElementSpace().getMesh() ==
            &getTestFunction().getFiniteElementSpace().getMesh());
      const auto& trialFes = getTrialFunction().getFiniteElementSpace();
      const auto& testFes = getTestFunction().getFiniteElementSpace();
      const auto& mesh = getTrialFunction().getFiniteElementSpace().getMesh();
      m_operator.reset(
            new OperatorType(
               getAssembly().execute({mesh, trialFes, testFes, getIntegrators()})));
   }
}

#endif
