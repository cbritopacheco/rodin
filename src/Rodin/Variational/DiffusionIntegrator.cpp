/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Problem.h"

#include "DiffusionIntegrator.h"

namespace Rodin::Variational
{
   DiffusionIntegrator::DiffusionIntegrator(const DiffusionIntegrator& other)
      : m_lambda(other.m_lambda->copy()), m_bf(other.m_bf)
   {}

   DiffusionIntegrator&
   DiffusionIntegrator::setBilinearForm(BilinearFormBase& bf)
   {
      m_bf.emplace(bf);
      return *this;
   }

   void DiffusionIntegrator::eval()
   {
      assert(m_bf);
      m_lambda->buildMFEMCoefficient();
      m_bf->get()
           .getHandle()
           .AddDomainIntegrator(
                 new mfem::DiffusionIntegrator(m_lambda->getMFEMCoefficient()));
   }

   DiffusionIntegrator& DiffusionIntegrator::toggleSign()
   {
      m_lambda.reset(new ScalarCoefficient(-(*m_lambda)));
      return *this;
   }

   DiffusionIntegrator*
   DiffusionIntegrator::copy() const noexcept
   {
      return new DiffusionIntegrator(*this);
   }
}

