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

   void MassIntegrator::build()
   {
      m_lambda->build();
      m_mfemBFI
         = std::make_unique<mfem::MassIntegrator>(
               m_lambda->get());
   }

   mfem::BilinearFormIntegrator&
   MassIntegrator::get()
   {
      assert(m_mfemBFI);
      return *m_mfemBFI;
   }

   mfem::BilinearFormIntegrator* MassIntegrator::release()
   {
      return m_mfemBFI.release();
   }
}

