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
  PointBase::PointBase(std::reference_wrapper<const Polytope> polytope, std::reference_wrapper<const PolytopeTransformation> trans,
      const Math::SpatialVector<Real>& pc)
    : m_polytopeStorage(PolytopeStorage::Reference),
      m_polytope(polytope), m_trans(trans), m_pc(pc)
  {}

  PointBase::PointBase(std::reference_wrapper<const Polytope> polytope, std::reference_wrapper<const PolytopeTransformation> trans)
    : m_polytopeStorage(PolytopeStorage::Reference),
      m_polytope(polytope), m_trans(trans)
  {}

  PointBase::PointBase(Polytope&& polytope, std::reference_wrapper<const PolytopeTransformation> trans,
      const Math::SpatialVector<Real>& pc)
    : m_polytopeStorage(PolytopeStorage::Value),
      m_polytope(std::move(polytope)), m_trans(trans), m_pc(pc)
  {}

  PointBase::PointBase(Polytope&& polytope, std::reference_wrapper<const PolytopeTransformation> trans)
    : m_polytopeStorage(PolytopeStorage::Value),
      m_polytope(std::move(polytope)), m_trans(trans)
  {}

  bool PointBase::operator<(const PointBase& p) const
  {
    assert(getDimension() == p.getDimension());
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
    : m_polytopeStorage(other.m_polytopeStorage),
      m_polytope(other.m_polytope),
      m_trans(other.m_trans),
      m_pc(other.m_pc),
      m_jacobian(other.m_jacobian),
      m_jacobianInverse(other.m_jacobianInverse),
      m_jacobianDeterminant(other.m_jacobianDeterminant),
      m_distortion(other.m_distortion)
  {}

  PointBase::PointBase(PointBase&& other)
    : m_polytopeStorage(std::move(other.m_polytopeStorage)),
      m_polytope(std::move(other.m_polytope)),
      m_trans(std::move(other.m_trans)),
      m_pc(std::move(other.m_pc)),
      m_jacobian(std::move(other.m_jacobian)),
      m_jacobianInverse(std::move(other.m_jacobianInverse)),
      m_jacobianDeterminant(std::move(other.m_jacobianDeterminant)),
      m_distortion(std::move(other.m_distortion))
  {}

  const Polytope& PointBase::getPolytope() const
  {
    if (m_polytopeStorage == PolytopeStorage::Value)
    {
      return std::get<const Polytope>(m_polytope);
    }
    else
    {
      assert(m_polytopeStorage == PolytopeStorage::Reference);
      return std::get<std::reference_wrapper<const Polytope>>(m_polytope);
    }
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
    if (!m_pc.read().has_value())
    {
      m_pc.write(
          [&](auto& obj)
          {
            obj.emplace(m_trans.get().transform(getReferenceCoordinates()));
          });
    }
    assert(m_pc.read().has_value());
    return m_pc.read().value();
  }

  const Math::SpatialMatrix<Real>& PointBase::getJacobian() const
  {
    if (!m_jacobian.read().has_value())
    {
      m_jacobian.write(
          [&](auto& obj)
          {
            obj.emplace(m_trans.get().jacobian(getReferenceCoordinates()));
          });
    }
    assert(m_jacobian.read().has_value());
    return m_jacobian.read().value();
  }

  const Math::SpatialMatrix<Real>& PointBase::getJacobianInverse() const
  {
    if (!m_jacobianInverse.read().has_value())
    {
      const auto& polytope = getPolytope();
      const size_t rdim = Polytope::getGeometryDimension(polytope.getGeometry());
      const size_t sdim = polytope.getMesh().getSpaceDimension();
      assert(rdim <= sdim);
      if (rdim == sdim)
      {
        switch (rdim)
        {
          case 1:
          {
            Math::SpatialMatrix<Real> inv(1, 1);
            inv.coeffRef(0, 0) = 1 / getJacobian().coeff(0, 0);
            m_jacobianDeterminant.write(
                [&](auto& obj)
                {
                  obj.emplace(getJacobian().coeff(0, 0));
                });
            assert(m_jacobianDeterminant.read().has_value());
            m_jacobianInverse.write(
                [&](auto& obj)
                {
                  obj.emplace(std::move(inv));
                });
            assert(m_jacobianInverse.read().has_value());
            return m_jacobianInverse.read().value();
          }
          case 2:
          {
            const auto& jac = getJacobian();
            const Real a = jac.coeff(0, 0);
            const Real b = jac.coeff(0, 1);
            const Real c = jac.coeff(1, 0);
            const Real d = jac.coeff(1, 1);
            const Real det = a * d - b * c;
            assert(det != 0);
            m_jacobianDeterminant.write(
                [&](auto& obj)
                {
                  obj.emplace(det);
                });
            assert(m_jacobianDeterminant.read().has_value());
            Math::SpatialMatrix<Real> inv(2, 2);
            inv.coeffRef(0, 0) = d / det;
            inv.coeffRef(0, 1) = -b / det;
            inv.coeffRef(1, 0) = -c / det;
            inv.coeffRef(1, 1) = a / det;
            m_jacobianInverse.write(
                [&](auto& obj)
                {
                  obj.emplace(std::move(inv));
                });
            assert(m_jacobianInverse.read().has_value());
            return m_jacobianInverse.read().value();
          }
          case 3:
          {
            const auto& jac = getJacobian();
            const Real a = jac.coeff(0, 0);
            const Real b = jac.coeff(0, 1);
            const Real c = jac.coeff(0, 2);
            const Real d = jac.coeff(1, 0);
            const Real e = jac.coeff(1, 1);
            const Real f = jac.coeff(1, 2);
            const Real g = jac.coeff(2, 0);
            const Real h = jac.coeff(2, 1);
            const Real i = jac.coeff(2, 2);
            const Real A = e * i - f * h;
            const Real B = -(d * i - f * g);
            const Real C = d * h - e * g;
            const Real D = -(b * i - c * h);
            const Real E = a * i - c * g;
            const Real F = -(a * h - b * g);
            const Real G = b * f - c * e;
            const Real H = - (a * f  - c * d);
            const Real I = a * e - b * d;

            const Real det = a * A + b * B + c * C;
            m_jacobianDeterminant.write(
                [&](auto& obj)
                {
                  obj.emplace(det);
                });
            assert(m_jacobianDeterminant.read().has_value());
            assert(det != 0);
            Math::SpatialMatrix<Real> inv(3, 3);
            inv.coeffRef(0, 0) = A / det;
            inv.coeffRef(0, 1) = D / det;
            inv.coeffRef(0, 2) = G / det;
            inv.coeffRef(1, 0) = B / det;
            inv.coeffRef(1, 1) = E / det;
            inv.coeffRef(1, 2) = H / det;
            inv.coeffRef(2, 0) = C / det;
            inv.coeffRef(2, 1) = F / det;
            inv.coeffRef(2, 2) = I / det;
            m_jacobianInverse.write(
                [&](auto& obj)
                {
                  obj.emplace(std::move(inv));
                });
            assert(m_jacobianInverse.read().has_value());
            return m_jacobianInverse.read().value();
          }
          default:
          {
            m_jacobianInverse.write(
                [&](auto& obj)
                {
                  obj.emplace(getJacobian().inverse());
                });
            assert(m_jacobianInverse.read().has_value());
            return m_jacobianInverse.read().value();
          }
        }
      }
      else
      {
        m_jacobianInverse.write(
            [&](auto& obj)
            {
              obj.emplace(getJacobian().completeOrthogonalDecomposition().pseudoInverse());
            });
        assert(m_jacobianInverse.read().has_value());
        return m_jacobianInverse.read().value();
      }
    }
    assert(m_jacobianInverse.read().has_value());
    return m_jacobianInverse.read().value();
  }

  Real PointBase::getJacobianDeterminant() const
  {
    if (!m_jacobianDeterminant.read().has_value())
    {
      const auto& jac = getJacobian();
      const auto rows = jac.rows();
      const auto cols = jac.cols();
      if (rows == cols)
      {
        switch (rows)
        {
          case 1:
          {
            m_jacobianDeterminant.write(
                [&](auto& obj)
                {
                  obj.emplace(jac.coeff(0, 0));
                });
            assert(m_jacobianDeterminant.read().has_value());
            return m_jacobianDeterminant.read().value();
          }
          case 2:
          {
            const Real a = jac.coeff(0, 0);
            const Real b = jac.coeff(0, 1);
            const Real c = jac.coeff(1, 0);
            const Real d = jac.coeff(1, 1);
            m_jacobianDeterminant.write(
                [&](auto& obj)
                {
                  obj.emplace(a * d - b * c);
                });
            assert(m_jacobianDeterminant.read().has_value());
            return m_jacobianDeterminant.read().value();
          }
          case 3:
          {
            const Real a = jac.coeff(0, 0);
            const Real b = jac.coeff(0, 1);
            const Real c = jac.coeff(0, 2);
            const Real d = jac.coeff(1, 0);
            const Real e = jac.coeff(1, 1);
            const Real f = jac.coeff(1, 2);
            const Real g = jac.coeff(2, 0);
            const Real h = jac.coeff(2, 1);
            const Real i = jac.coeff(2, 2);
            const Real A = e * i - f * h;
            const Real B = -(d * i - f * g);
            const Real C = d * h - e * g;
            m_jacobianDeterminant.write(
                [&](auto& obj)
                {
                  obj.emplace(a * A + b * B + c * C);
                });
            assert(m_jacobianDeterminant.read().has_value());
            return m_jacobianDeterminant.read().value();
          }
          default:
          {
            m_jacobianDeterminant.write(
                [&](auto& obj)
                {
                  obj.emplace(jac.determinant());
                });
            assert(m_jacobianDeterminant.read().has_value());
            return m_jacobianDeterminant.read().value();
          }
        }
      }
      else
      {
        m_jacobianDeterminant.write(
            [&](auto& obj)
            {
              obj.emplace(Math::sqrt((jac.transpose() * jac).determinant()));
            });
        assert(m_jacobianDeterminant.read().has_value());
        return m_jacobianDeterminant.read().value();
      }
    }
    assert(m_jacobianDeterminant.read().has_value());
    return m_jacobianDeterminant.read().value();
  }

  Real PointBase::getDistortion() const
  {
    if (!m_distortion.read().has_value())
    {
      const auto& jac = getJacobian();
      const auto rows = jac.rows();
      const auto cols = jac.cols();
      if (rows == cols)
      {
        m_distortion.write(
            [&](auto& obj)
            {
              obj.emplace(getJacobianDeterminant());
            });
        assert(m_distortion.read().has_value());
        return m_distortion.read().value();
      }
      else
      {
        if (jac.rows() == 2 && jac.cols() == 1)
        {
          m_distortion.write(
              [&](auto& obj)
              {
                obj.emplace(Math::sqrt(jac.coeff(0) * jac.coeff(0) + jac.coeff(1) * jac.coeff(1)));
              });
          assert(m_distortion.read().has_value());
          return m_distortion.read().value();
        }
        else
        {
          m_distortion.write(
              [&](auto& obj)
              {
                obj.emplace(Math::sqrt(Math::abs((jac.transpose() * jac).determinant())));
              });
          assert(m_distortion.read().has_value());
          return m_distortion.read().value();
        }
      }
    }
    assert(m_distortion.read().has_value());
    return m_distortion.read().value();
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

  Point::Point(
      std::reference_wrapper<const Polytope> polytope,
      std::reference_wrapper<const PolytopeTransformation> trans,
      std::reference_wrapper<const Math::SpatialVector<Real>> rc,
      const Math::SpatialVector<Real>& pc)
    : PointBase(polytope, trans, pc), m_rcStorage(RCStorage::Reference), m_rc(rc)
  {}

  Point::Point(
      std::reference_wrapper<const Polytope> polytope,
      std::reference_wrapper<const PolytopeTransformation> trans,
      Math::SpatialVector<Real>&& rc,
      const Math::SpatialVector<Real>& pc)
    : PointBase(polytope, trans, pc), m_rcStorage(RCStorage::Value), m_rc(std::move(rc))
  {}

  Point::Point(
      std::reference_wrapper<const Polytope> polytope,
      std::reference_wrapper<const PolytopeTransformation> trans,
      std::reference_wrapper<const Math::SpatialVector<Real>> rc)
    : PointBase(polytope, trans), m_rcStorage(RCStorage::Reference), m_rc(rc)
  {}

  Point::Point(
      std::reference_wrapper<const Polytope> polytope,
      std::reference_wrapper<const PolytopeTransformation> trans,
      Math::SpatialVector<Real>&& rc)
    : PointBase(polytope, trans), m_rcStorage(RCStorage::Value), m_rc(std::move(rc))
  {}

  Point::Point(
      Polytope&& polytope,
      std::reference_wrapper<const PolytopeTransformation> trans,
      std::reference_wrapper<const Math::SpatialVector<Real>> rc,
      const Math::SpatialVector<Real>& pc)
    : PointBase(std::move(polytope), trans, pc), m_rcStorage(RCStorage::Reference), m_rc(rc)
  {}

  Point::Point(
      Polytope&& polytope,
      std::reference_wrapper<const PolytopeTransformation> trans,
      Math::SpatialVector<Real>&& rc,
      const Math::SpatialVector<Real>& pc)
    : PointBase(std::move(polytope), trans, pc), m_rcStorage(RCStorage::Value), m_rc(std::move(rc))
  {}

  Point::Point(
      Polytope&& polytope,
      std::reference_wrapper<const PolytopeTransformation> trans,
      std::reference_wrapper<const Math::SpatialVector<Real>> rc)
    : PointBase(std::move(polytope), trans), m_rcStorage(RCStorage::Reference), m_rc(rc)
  {}

  Point::Point(
      Polytope&& polytope,
      std::reference_wrapper<const PolytopeTransformation> trans,
      Math::SpatialVector<Real>&& rc)
    : PointBase(std::move(polytope), trans), m_rcStorage(RCStorage::Value), m_rc(std::move(rc))
  {}

  Point::Point(const Point& other)
    : PointBase(other),
      m_rcStorage(other.m_rcStorage),
      m_rc(other.m_rc)
  {}

  Point::Point(Point&& other)
    : PointBase(std::move(other)),
      m_rcStorage(std::move(other.m_rcStorage)),
      m_rc(std::move(other.m_rc))
  {}

  const Math::SpatialVector<Real>& Point::getReferenceCoordinates() const
  {
    if (m_rcStorage == RCStorage::Value)
    {
      return std::get<const Math::SpatialVector<Real>>(m_rc);
    }
    else
    {
      assert(m_rcStorage == RCStorage::Reference);
      return std::get<std::reference_wrapper<const Math::SpatialVector<Real>>>(m_rc);
    }
  }
}
