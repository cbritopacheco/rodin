/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "FormLanguage/ScalarUnaryMinus.h"

#include "Problem.h"

#include "DiffusionIntegrator.h"

namespace Rodin::Variational
{
   DiffusionIntegrator::DiffusionIntegrator()
      : DiffusionIntegrator(ScalarCoefficient(1.0))
   {}

   DiffusionIntegrator::DiffusionIntegrator(const ScalarCoefficientBase& lambda)
      : m_lambda(lambda.copy())
   {}

   DiffusionIntegrator::DiffusionIntegrator(const DiffusionIntegrator& other)
      :  m_attr(other.m_attr),
         m_lambda(other.m_lambda->copy())
   {}

   void DiffusionIntegrator::buildMFEMBilinearFormIntegrator()
   {
      m_lambda->buildMFEMCoefficient();
      m_bfi = std::make_unique<mfem::DiffusionIntegrator>(
            m_lambda->getMFEMCoefficient());
   }

   mfem::BilinearFormIntegrator&
   DiffusionIntegrator::getMFEMBilinearFormIntegrator()
   {
      assert(m_bfi);
      return *m_bfi;
   }

   mfem::BilinearFormIntegrator*
   DiffusionIntegrator::releaseMFEMBilinearFormIntegrator()
   {
      return m_bfi.release();
   }
}

