#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      VectorCoefficient::VectorCoefficient(const VectorCoefficientBase& v)
         :  mfem::VectorCoefficient(v.getDimension()),
            m_v(v.copy())
      {}

      void VectorCoefficient::Eval(
            mfem::Vector& value,
            mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
      {
         if (trans.ElementType == mfem::ElementTransformation::BDR_ELEMENT
               && trans.mesh->FaceIsInterior((trans.mesh->GetBdrFace(trans.ElementNo))))
         {
            m_v->getValueOnInteriorBoundary(value, trans, ip);
         }
         else
         {
            m_v->getValue(value, trans, ip);
         }
      }
   }

   void VectorCoefficientBase::getValueOnInteriorBoundary(
         mfem::Vector& value, mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
      const
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
}

