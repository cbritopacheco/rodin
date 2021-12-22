/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "FormLanguage/ScalarUnaryMinus.h"

#include "Problem.h"

#include "VectorDiffusionIntegrator.h"

namespace Rodin::Variational
{
   VectorDiffusionIntegrator::VectorDiffusionIntegrator()
      : VectorDiffusionIntegrator(ScalarCoefficient(1.0))
   {}

   VectorDiffusionIntegrator::VectorDiffusionIntegrator(const ScalarCoefficientBase& lambda)
      : m_lambda(lambda.copy())
   {}

   VectorDiffusionIntegrator::VectorDiffusionIntegrator(const VectorDiffusionIntegrator& other)
      : m_lambda(other.m_lambda->copy())
   {}

   void VectorDiffusionIntegrator::buildMFEMBilinearFormIntegrator()
   {
      m_lambda->buildMFEMCoefficient();
      m_bfi = std::make_unique<mfem::VectorDiffusionIntegrator>(
            m_lambda->getMFEMCoefficient());
   }

   mfem::BilinearFormIntegrator&
   VectorDiffusionIntegrator::getMFEMBilinearFormIntegrator()
   {
      assert(m_bfi);
      return *m_bfi;
   }

   mfem::BilinearFormIntegrator*
   VectorDiffusionIntegrator::releaseMFEMBilinearFormIntegrator()
   {
      return m_bfi.release();
   }
}


