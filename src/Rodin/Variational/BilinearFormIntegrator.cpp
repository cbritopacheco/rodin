#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      BilinearFormIntegrator::BilinearFormIntegrator(const BilinearFormIntegratorBase& bfi)
         : m_bfi(bfi.copy())
      {}

      void BilinearFormIntegrator::AssembleElementMatrix(
            const mfem::FiniteElement& fe,
            mfem::ElementTransformation& trans, mfem::DenseMatrix& mat)
      {
         m_bfi->getElementMatrix(fe, trans, mat);
      }

      void BilinearFormIntegrator::AssembleElementMatrix2(
            const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
            mfem::ElementTransformation& trans, mfem::DenseMatrix& mat)
      {
         m_bfi->getElementMatrix(trial, test, trans, mat);
      }
   }
}
