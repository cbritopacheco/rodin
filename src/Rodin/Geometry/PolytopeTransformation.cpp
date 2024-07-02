#include <Eigen/Dense>

#include "Rodin/Variational/FiniteElement.h"

#include "Mesh.h"
#include "Polytope.h"

#include "PolytopeTransformation.h"

namespace Rodin::Geometry
{
  void PolytopeTransformation::inverse(const Math::SpatialVector<Real>& pc, Math::SpatialVector<Real>& rc) const
  {
    const size_t pdim = getPhysicalDimension();
    const size_t rdim = getReferenceDimension();
    assert(pc.size() >= 0);
    assert(static_cast<size_t>(pc.size()) == pdim);
    // All elements have 0 as reference coordinate
    Math::SpatialVector<Real> rc0 = Math::SpatialVector<Real>::Zero(rdim);
    Math::SpatialVector<Real> pc0 = transform(rc0);
    Math::SpatialMatrix<Real> jac = jacobian(rc0);
    if (rdim == pdim)
      rc = rc0 + jac.partialPivLu().solve(pc - pc0);
    else
      rc = rc0 + (jac.transpose() * jac).partialPivLu().solve(jac.transpose() * (pc - pc0));
  }
}
