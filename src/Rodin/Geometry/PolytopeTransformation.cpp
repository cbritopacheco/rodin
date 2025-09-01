#include <Eigen/Dense>

#include "Rodin/Variational/FiniteElement.h"

#include "Mesh.h"
#include "Polytope.h"

#include "PolytopeTransformation.h"

namespace Rodin::Geometry
{
  void PolytopeTransformation::inverse(Math::SpatialVector<Real>& rc, const Math::SpatialVector<Real>& pc) const
  {
    static thread_local Math::SpatialVector<Real> s_zero;
    static thread_local Math::SpatialVector<Real> s_pc0;
    static thread_local Math::SpatialMatrix<Real> s_jac;

    const size_t pdim = getPhysicalDimension();
    const size_t rdim = getReferenceDimension();

    assert(pc.size() >= 0);
    assert(static_cast<size_t>(pc.size()) == pdim);

    s_zero.resize(rdim);
    s_zero.setZero();

    this->transform(s_pc0, s_zero);
    this->jacobian(s_jac, s_zero);

    if (rdim == pdim)
      rc = s_jac.partialPivLu().solve(pc - s_pc0);
    else
      rc = (s_jac.transpose() * s_jac).partialPivLu().solve(s_jac.transpose() * (pc - s_pc0));
  }
}
