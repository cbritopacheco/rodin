#include "Rodin/Alert.h"

#include "Transpose.h"

#include "MatrixCoefficient.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      MatrixCoefficient::MatrixCoefficient(const MatrixCoefficientBase& mat)
         :  mfem::MatrixCoefficient(mat.getRows(), mat.getColumns()),
            m_mat(mat.copy())
      {}

      void MatrixCoefficient::Eval(
            mfem::DenseMatrix& value,
            mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
      {
         if (trans.ElementType == mfem::ElementTransformation::BDR_ELEMENT
               && trans.mesh->FaceIsInterior((trans.mesh->GetBdrFace(trans.ElementNo))))
         {
            m_mat->getValueOnInteriorBoundary(value, trans, ip);
         }
         else
         {
            m_mat->getValue(value, trans, ip);
         }
      }
   }

   void MatrixCoefficientBase::getValueOnInteriorBoundary(
         mfem::DenseMatrix& value, mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   {
      std::optional<int> traceDomain = getTraceDomain();
      if (!traceDomain)
      {
         // Attempt to extend the coefficient to the whole domain. This will
         // probably work for constant coefficients and piecewise coefficients
         // which are defined everywhere.
         getValue(value, trans, ip);
      }
      else
      {
         // Default behaviour is to extend the values on the trace domain up to
         // the interior boundary.
         mfem::FaceElementTransformations* ft =
            trans.mesh->GetFaceElementTransformations(trans.mesh->GetBdrFace(trans.ElementNo));
         ft->SetAllIntPoints(&ip);
         if (ft->GetElement1Transformation().Attribute == *m_traceDomain)
            return getValue(value, ft->GetElement1Transformation(), ip);
         else if (ft->GetElement2Transformation().Attribute == *m_traceDomain)
            return getValue(value, ft->GetElement2Transformation(), ip);
         else
         {
            // The boundary over which we are evaluating must be the interface
            // between the trace domain and some other domain, i.e. it is not
            // the boundary that was specified!
            Alert::Exception()
               << "Invalid boundary for trace domain " << *traceDomain
               << Alert::Raise;
         }
      }
      value = NAN;
      return;
   }

   Transpose MatrixCoefficientBase::T() const
   {
      return Transpose(*this);
   }
}

