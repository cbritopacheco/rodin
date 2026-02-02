#include <Eigen/Dense>

#include "Rodin/Variational/FiniteElement.h"

#include "Mesh.h"
#include "Polytope.h"

#include "PolytopeTransformation.h"

namespace Rodin::Geometry
{
  void PolytopeTransformation::inverse(Math::SpatialVector<Real>& rc, const Math::SpatialVector<Real>& pc) const
  {
    const size_t pdim = this->getPhysicalDimension();
    const size_t rdim = this->getReferenceDimension();

    assert(pc.size() >= 0);
    assert(static_cast<size_t>(pc.size()) == pdim);

    Math::SpatialPoint zero(rdim);
    zero.setZero();

    Math::SpatialPoint pc0;
    Math::SpatialMatrix<Real> jac;

    this->transform(pc0, zero);
    this->jacobian(jac, zero);

    if (rdim == pdim)
      rc = jac.solve(pc - pc0);
    else
      rc = (jac.transpose() * jac).solve(jac.transpose() * (pc - pc0));
  }
}
