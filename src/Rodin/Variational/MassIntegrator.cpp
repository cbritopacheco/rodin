#include "ScalarCoefficient.h"

#include "MassIntegrator.h"

namespace Rodin::Variational
{
   MassIntegrator::MassIntegrator()
      : MassIntegrator(ScalarCoefficient(1.0))
   {}

   MassIntegrator::MassIntegrator(const ScalarCoefficientBase& lambda)
      : m_lambda(lambda.copy())
   {}

   MassIntegrator::MassIntegrator(const MassIntegrator& other)
      :  m_attr(other.m_attr),
         m_lambda(other.m_lambda->copy())
   {}

   void MassIntegrator::buildMFEMBilinearFormIntegrator()
   {
      m_lambda->buildMFEMCoefficient();
      m_mfemBFI
         = std::make_unique<mfem::MassIntegrator>(
               m_lambda->getMFEMCoefficient());
   }

   mfem::BilinearFormIntegrator&
   MassIntegrator::getMFEMBilinearFormIntegrator()
   {
      assert(m_mfemBFI);
      return *m_mfemBFI;
   }

   mfem::BilinearFormIntegrator* MassIntegrator::releaseMFEMBilinearFormIntegrator()
   {
      return m_mfemBFI.release();
   }
}

