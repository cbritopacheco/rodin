#include "ScalarCoefficient.h"

#include "VectorMassIntegrator.h"

namespace Rodin::Variational
{
   VectorMassIntegrator::VectorMassIntegrator()
      : VectorMassIntegrator(ScalarCoefficient(1.0))
   {}

   VectorMassIntegrator::VectorMassIntegrator(const ScalarCoefficientBase& lambda)
      : m_lambda(lambda.copy())
   {}

   VectorMassIntegrator::VectorMassIntegrator(const VectorMassIntegrator& other)
      :  m_attr(other.m_attr),
         m_lambda(other.m_lambda->copy())
   {}

   void VectorMassIntegrator::buildMFEMBilinearFormIntegrator()
   {
      m_lambda->buildMFEMCoefficient();
      m_mfemBFI
         = std::make_unique<mfem::VectorMassIntegrator>(
               m_lambda->getMFEMCoefficient());
   }

   mfem::BilinearFormIntegrator&
   VectorMassIntegrator::getMFEMBilinearFormIntegrator()
   {
      assert(m_mfemBFI);
      return *m_mfemBFI;
   }

   mfem::BilinearFormIntegrator* VectorMassIntegrator::releaseMFEMBilinearFormIntegrator()
   {
      return m_mfemBFI.release();
   }
}


