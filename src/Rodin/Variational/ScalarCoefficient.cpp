/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Restriction.h"

#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
  //  double ScalarCoefficientBase::getValueOnInteriorBoundary(
  //        mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
  //  const
  //  {
  //     std::optional<int> traceDomain = getTraceDomain();
  //     if (!traceDomain)
  //     {
  //        // Attempt to extend the coefficient to the whole domain. This will
  //        // probably work for constant coefficients and piecewise coefficients
  //        // which are defined everywhere.
  //        return getValue(trans, ip);
  //     }
  //     else
  //     {
  //        // Default behaviour is to extend the values on the trace domain up to
  //        // the interior boundary.
  //        mfem::FaceElementTransformations* ft =
  //           trans.mesh->GetFaceElementTransformations(trans.mesh->GetBdrFace(trans.ElementNo));
  //        ft->SetAllIntPoints(&ip);
  //        if (ft->GetElement1Transformation().Attribute == *m_traceDomain)
  //           return getValue(ft->GetElement1Transformation(), ip);
  //        else if (ft->GetElement2Transformation().Attribute == *m_traceDomain)
  //           return getValue(ft->GetElement2Transformation(), ip);
  //        else
  //        {
  //           // The boundary over which we are evaluating must be the interface
  //           // between the trace domain and some other domain, i.e. it is not
  //           // the boundary that was specified!
  //           Alert::Exception()
  //              << "Invalid boundary for trace domain " << *traceDomain
  //              << Alert::Raise;
  //        }
  //     }
  //     return NAN;
  //  }

   std::unique_ptr<mfem::Coefficient> ScalarCoefficientBase::build() const
   {
      return std::unique_ptr<mfem::Coefficient>(new Internal::ProxyScalarCoefficient(*this));
   }

   Restriction<ScalarCoefficientBase> ScalarCoefficientBase::restrictTo(int attr)
   {
      return restrictTo(std::set<int>{attr});
   }

   Restriction<ScalarCoefficientBase> ScalarCoefficientBase::restrictTo(
         const std::set<int>& attrs)
   {
      return Restriction<ScalarCoefficientBase>(*this).to(attrs);
   }
}
