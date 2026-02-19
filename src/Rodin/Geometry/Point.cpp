#include <Eigen/Cholesky>

#include "Rodin/Configure.h"

#include "Rodin/Variational/QuadratureRule.h"

#include "Mesh.h"
#include "Polytope.h"
#include "PolytopeTransformation.h"

#include "Point.h"

namespace Rodin::Geometry
{
  // ---- Point --------------------------------------------------------------
  PointBase::PointBase(const Polytope& polytope)
    : m_polytope(polytope)
  {}

  PointBase::PointBase(const Polytope& polytope, const Math::SpatialPoint& pc)
    : m_polytope(polytope), m_pc(pc)
  {}

  PointBase::PointBase(const Polytope& polytope, Math::SpatialPoint&& pc)
    : m_polytope(polytope), m_pc(std::move(pc))
  {}

  bool PointBase::operator<(const PointBase& p) const
  {
    assert(this->getDimension() == p.getDimension());
    const auto& lhs = getCoordinates(Coordinates::Physical);
    const auto& rhs = p.getCoordinates(Coordinates::Physical);
    for (int i = 0; i < lhs.size() - 1; i++)
    {
      if (lhs(i) < rhs(i))
        return true;
      if (rhs(i) > lhs(i))
        return false;
    }
    return (lhs(lhs.size() - 1) < rhs(rhs.size() - 1));
  }

  PointBase::PointBase(const PointBase& other)
    : m_polytope(other.m_polytope),
      m_pc(other.m_pc),
      m_jacobian(other.m_jacobian),
      m_jacobianInverse(other.m_jacobianInverse),
      m_jacobianDeterminant(other.m_jacobianDeterminant),
      m_distortion(other.m_distortion)
  {}

  PointBase::PointBase(PointBase&& other)
    : m_polytope(std::move(other.m_polytope)),
      m_pc(std::move(other.m_pc)),
      m_jacobian(std::move(other.m_jacobian)),
      m_jacobianInverse(std::move(other.m_jacobianInverse)),
      m_jacobianDeterminant(std::move(other.m_jacobianDeterminant)),
      m_distortion(std::move(other.m_distortion))
  {}

  const Polytope& PointBase::getPolytope() const
  {
    return m_polytope;
  }

  PointBase& PointBase::setPolytope(const Polytope& polytope)
  {
    m_polytope = polytope;
    m_pc.reset();
    m_jacobian.reset();
    m_jacobianInverse.reset();
    m_jacobianDeterminant.reset();
    m_distortion.reset();
    return *this;
  }

  const Math::SpatialVector<Real>& PointBase::getCoordinates(Coordinates coords) const
  {
    if (coords == Coordinates::Physical)
    {
      return getPhysicalCoordinates();
    }
    else
    {
      assert(coords == Coordinates::Reference);
      return getReferenceCoordinates();
    }
  }

  const Math::SpatialVector<Real>& PointBase::getPhysicalCoordinates() const
  {
    if (!m_pc)
    {
      auto& pc = m_pc.emplace();
      this->getPolytope().getTransformation().transform(pc, this->getReferenceCoordinates());
    }
    return *m_pc;
  }

  const Math::SpatialMatrix<Real>& PointBase::getJacobian() const
  {
    if (!m_jacobian)
    {
      auto& jac = m_jacobian.emplace();
      this->getPolytope().getTransformation().jacobian(jac, this->getReferenceCoordinates());
      return jac;
    }
    assert(m_jacobian);
    return *m_jacobian;
  }

  const Math::SpatialMatrix<Real>& PointBase::getJacobianInverse() const
  {
    if (!m_jacobianInverse)
    {
      const auto& polytope = this->getPolytope();
      const size_t rdim = Polytope::Traits(polytope.getGeometry()).getDimension();
      const size_t sdim = polytope.getMesh().getSpaceDimension();
      assert(rdim <= sdim);
      if (rdim == sdim)
      {
        switch (rdim)
        {
          case 1:
          {
            const auto& jac = this->getJacobian();
            const auto& det = m_jacobianDeterminant.emplace(jac(0, 0));
            auto& inv = m_jacobianInverse.emplace(1, 1);
            inv(0, 0) = 1 / det;
            return inv;
          }
          case 2:
          {
            const auto& jac = this->getJacobian();
            const Real a = jac(0, 0);
            const Real b = jac(0, 1);
            const Real c = jac(1, 0);
            const Real d = jac(1, 1);
            const auto& det = m_jacobianDeterminant.emplace(a * d - b * c);
            assert(det != 0);

            auto& inv = m_jacobianInverse.emplace(2, 2);
            inv(0, 0) = d / det;
            inv(0, 1) = -b / det;
            inv(1, 0) = -c / det;
            inv(1, 1) = a / det;
            return inv;
          }
          case 3:
          {
            const auto& jac = this->getJacobian();
            const Real a = jac(0, 0);
            const Real b = jac(0, 1);
            const Real c = jac(0, 2);
            const Real d = jac(1, 0);
            const Real e = jac(1, 1);
            const Real f = jac(1, 2);
            const Real g = jac(2, 0);
            const Real h = jac(2, 1);
            const Real i = jac(2, 2);
            const Real A = e * i - f * h;
            const Real B = -(d * i - f * g);
            const Real C = d * h - e * g;
            const Real D = -(b * i - c * h);
            const Real E = a * i - c * g;
            const Real F = -(a * h - b * g);
            const Real G = b * f - c * e;
            const Real H = - (a * f  - c * d);
            const Real I = a * e - b * d;
            const auto& det = m_jacobianDeterminant.emplace(a * A + b * B + c * C);
            assert(det != 0);

            auto& inv = m_jacobianInverse.emplace(3, 3);
            inv(0, 0) = A / det;
            inv(0, 1) = D / det;
            inv(0, 2) = G / det;
            inv(1, 0) = B / det;
            inv(1, 1) = E / det;
            inv(1, 2) = H / det;
            inv(2, 0) = C / det;
            inv(2, 1) = F / det;
            inv(2, 2) = I / det;
            return inv;
          }
          default:
          {
            return m_jacobianInverse.emplace(this->getJacobian().inverse());
          }
        }
      }
      else
      {
        return m_jacobianInverse.emplace(this->getJacobian().pseudoInverse());
      }
    }
    assert(m_jacobianInverse);
    return *m_jacobianInverse;
  }

  Real PointBase::getJacobianDeterminant() const
  {
    if (!m_jacobianDeterminant)
    {
      const auto& jac = this->getJacobian();
      const auto rows = jac.rows();
      const auto cols = jac.cols();
      if (rows == cols)
      {
        switch (rows)
        {
          case 1:
          {
            return m_jacobianDeterminant.emplace(jac(0, 0));
          }
          case 2:
          {
            const Real a = jac(0, 0);
            const Real b = jac(0, 1);
            const Real c = jac(1, 0);
            const Real d = jac(1, 1);
            return m_jacobianDeterminant.emplace(a * d - b * c);
          }
          case 3:
          {
            const Real a = jac(0, 0);
            const Real b = jac(0, 1);
            const Real c = jac(0, 2);
            const Real d = jac(1, 0);
            const Real e = jac(1, 1);
            const Real f = jac(1, 2);
            const Real g = jac(2, 0);
            const Real h = jac(2, 1);
            const Real i = jac(2, 2);
            const Real A = e * i - f * h;
            const Real B = -(d * i - f * g);
            const Real C = d * h - e * g;
            return m_jacobianDeterminant.emplace(a * A + b * B + c * C);
          }
          default:
          {
            return m_jacobianDeterminant.emplace(jac.determinant());
          }
        }
      }
      else
      {
        return m_jacobianDeterminant.emplace(
            Math::sqrt((jac.transpose() * jac).determinant()));
      }
    }
    assert(m_jacobianDeterminant);
    return *m_jacobianDeterminant;
  }

  Real PointBase::getDistortion() const
  {
    if (m_distortion)
      return *m_distortion;

    const auto& J = this->getJacobian();
    const int m = J.rows(); // physical dim
    const int n = J.cols(); // reference dim

    // Expected: 0 <= m,n <= 3 (as per your constraints)
    Real dist = 0;

    // 0D measure convention: a point has measure 1
    if (n == 0)
    {
      dist = 1;
      return m_distortion.emplace(dist);
    }

    // Prefer explicit low-dim formulas where they are cheap/stable.
    if (m == n)
    {
      // Square: length/area/volume element = |det J|
      // (Eigen defines det of 1x1,2x2,3x3; 0x0 handled above)
      dist = Math::abs(this->getJacobianDeterminant());
      return m_distortion.emplace(dist);
    }

    // Embedded measures with explicit formulas (more stable than det(J^T J))
    if (n == 1)
    {
      // Curve in mD: ||dX/dξ||
      dist = (m > 0) ? J.col(0).norm() : Real(0);
      return m_distortion.emplace(dist);
    }

    if (n == 2)
    {
      if (m == 3)
      {
        // Surface in 3D: ||a x b||
        const auto a = J.col(0);
        const auto b = J.col(1);
        dist = a.cross(b).norm();
        return m_distortion.emplace(dist);
      }
      // For m==0,1 the generic path below yields 0 anyway (rank deficiency),
      // for m==2 it would have been caught by m==n above.
    }

    // Fully generic k-dimensional measure: sqrt(det(J^T J))
    // Works for any 0<=m,n<=3, including m<n (det will go to 0 if rank-deficient).
    const Real detG = (J.transpose() * J).determinant();

    // Numerical guard: detG should be >= 0; clamp tiny negatives due to roundoff
    dist = Math::sqrt(detG >= 0 ? detG : 0);

    return m_distortion.emplace(dist);
  }

  size_t PointBase::getDimension(Coordinates coords) const
  {
    const auto& polytope = getPolytope();
    switch (coords)
    {
      case Coordinates::Physical:
        return polytope.getMesh().getSpaceDimension();
      case Coordinates::Reference:
        return polytope.getMesh().getDimension();
      default:
      {
        assert(false);
        return 0;
      }
    }
  }

  Point::Point(const Polytope& polytope, const Math::SpatialPoint& rc)
    : PointBase(polytope),
      m_rc(std::cref(rc))
  {}

  Point::Point(const Polytope& polytope, Math::SpatialPoint&& rc)
    : PointBase(polytope),
      m_rc(std::move(rc))
  {}

  Point::Point(const Polytope& polytope, const Math::SpatialPoint& rc, const Math::SpatialPoint& pc)
    : PointBase(polytope, pc),
      m_rc(std::cref(rc))
  {}

  Point::Point(const Polytope& polytope, const Math::SpatialPoint& rc, Math::SpatialPoint&& pc)
    : PointBase(polytope, std::move(pc)),
      m_rc(std::cref(rc))
  {}

  Point::Point(const Polytope& polytope, Math::SpatialPoint&& rc, const Math::SpatialPoint& pc)
    : PointBase(polytope, pc),
      m_rc(std::move(rc))
  {}

  Point::Point(const Polytope& polytope, Math::SpatialPoint&& rc, Math::SpatialPoint&& pc)
    : PointBase(polytope, std::move(pc)),
      m_rc(std::move(rc))
  {}

  Point::Point(const Point& other)
    : PointBase(other),
      m_rc(other.m_rc)
  {}

  Point::Point(Point&& other)
    : PointBase(std::move(other)),
      m_rc(std::move(other.m_rc))
  {}

  const Math::SpatialVector<Real>& Point::getReferenceCoordinates() const
  {
    if (std::holds_alternative<const Math::SpatialPoint>(m_rc))
    {
      return std::get<const Math::SpatialPoint>(m_rc);
    }
    else
    {
      return std::get<std::reference_wrapper<Math::SpatialPoint>>(m_rc).get();
    }
  }
}
