/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "LinearFormIntegratorUnaryMinus.h"

#include "LinearFormIntegratorSum.h"

namespace Rodin::Variational::FormLanguage
{
   LinearFormIntegratorSum::LinearFormIntegratorSum(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorBase& rhs)
      : m_lhs(lhs.copy()), m_rhs(rhs.copy())
   {}

   LinearFormIntegratorSum::LinearFormIntegratorSum(const LinearFormIntegratorSum& other)
      :  m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
   {}

   LinearFormIntegratorBase& LinearFormIntegratorSum::getLHS()
   {
      return *m_lhs;
   }

   LinearFormIntegratorBase& LinearFormIntegratorSum::getRHS()
   {
      return *m_rhs;
   }

   void LinearFormIntegratorSum::buildMFEMLinearFormIntegrator()
   {
      m_lhs->buildMFEMLinearFormIntegrator();
      m_rhs->buildMFEMLinearFormIntegrator();
      m_mfemLFI = std::make_unique<Internal::LinearFormIntegratorSum>(
            m_lhs->getMFEMLinearFormIntegrator(), m_rhs->getMFEMLinearFormIntegrator());
   }

   mfem::LinearFormIntegrator*
   LinearFormIntegratorSum::releaseMFEMLinearFormIntegrator()
   {
      assert(m_mfemLFI);
      return m_mfemLFI.release();
   }

   mfem::LinearFormIntegrator&
   LinearFormIntegratorSum::getMFEMLinearFormIntegrator()
   {
      assert(m_mfemLFI);
      return *m_mfemLFI;
   }

   LinearFormIntegratorSum operator+(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorBase& rhs)
   {
      return LinearFormIntegratorSum(lhs, rhs);
   }

   LinearFormIntegratorSum operator-(
         const LinearFormIntegratorBase& lhs, const LinearFormIntegratorBase& rhs)
   {
      return LinearFormIntegratorSum(lhs, LinearFormIntegratorUnaryMinus(rhs));
   }
}

