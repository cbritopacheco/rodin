#include <Eigen/Cholesky>

#include "Rodin/Variational/QuadratureRule.h"

#include "Mesh.h"
#include "SimplexTransformation.h"

#include "Simplex.h"

namespace Rodin::Geometry
{
  bool operator<(const Polytope& lhs, const Polytope& rhs)
  {
    return lhs.getIndex() < rhs.getIndex();
  }

  // ---- Simplex -----------------------------------------------------------
  Polytope::Polytope(size_t dimension, Index index, const MeshBase& mesh)
    : m_dimension(dimension), m_index(index), m_mesh(mesh)
  {}

  Attribute Polytope::getAttribute() const
  {
    return getMesh().getAttribute(getDimension(), getIndex());
  }

  Polytope::Geometry Polytope::getGeometry() const
  {
    return getMesh().getGeometry(getDimension(), getIndex());
  }

  VertexIterator Polytope::getVertex() const
  {
    const auto& vertices = getVertices();
    return VertexIterator(
        getMesh(), IteratorIndexGenerator(vertices.begin(), vertices.end()));
  }

  const Array<Index>& Polytope::getVertices() const
  {
    return m_mesh.get().getConnectivity().getPolytope(getDimension(), getIndex());
  }

  PolytopeIterator Polytope::getAdjacent() const
  {
    assert(false);
    return PolytopeIterator(0, getMesh(), EmptyIndexGenerator());
  }

  PolytopeIterator Polytope::getIncident() const
  {
    assert(false);
    return PolytopeIterator(0, getMesh(), EmptyIndexGenerator());
  }

  const PolytopeTransformation& Polytope::getTransformation() const
  {
    return m_mesh.get().getPolytopeTransformation(m_dimension, m_index);
  }

  Scalar Polytope::getVolume() const
  {
    assert(false);
    return 0;
    // const Variational::QuadratureRule& qr =
    //   Variational::QuadratureRule::get(*this, trans.OrderJ());
    // Scalar volume = 0.0;
    // for (size_t i = 0; i < qr.size(); i++)
    //   volume += qr.getWeight(i) * trans.Weight();
    // return volume;
  }

  // ---- Element -----------------------------------------------------------
  Element::Element(Index index, const MeshBase& mesh)
    : Polytope(mesh.getDimension(), index, mesh)
  {}

  // ---- Face --------------------------------------------------------------
  Face::Face(Index index, const MeshBase& mesh)
    : Polytope(mesh.getDimension() - 1, index, mesh)
  {}

  bool Face::isBoundary() const
  {
    return getMesh().isBoundary(getIndex());
  }

  bool Face::isInterface() const
  {
    return getMesh().isInterface(getIndex());
  }

  // ---- Vertex -------------------------------------------------------------
  Vertex::Vertex(Index index, const MeshBase& mesh)
    : Polytope(0, index, mesh)
  {}

  Eigen::Map<const Math::Vector> Vertex::getCoordinates() const
  {
    return getMesh().getVertexCoordinates(getIndex());
  }

  // ---- Point --------------------------------------------------------------
  Point::Point(const Polytope& simplex, const PolytopeTransformation& trans, const Math::Vector& rc)
    : m_polytope(simplex), m_trans(trans), m_rc(rc)
  {}

  Point::Point(const Polytope& simplex, const PolytopeTransformation& trans, const Math::Vector& rc, const Math::Vector& pc)
    : m_polytope(simplex), m_trans(trans), m_rc(rc), m_pc(pc)
  {}

  const Math::Vector& Point::getCoordinates(Coordinates coords) const
  {
    switch (coords)
    {
      case Coordinates::Physical:
      {
        if (!m_pc.has_value())
        {
          assert(m_rc.has_value());
          m_pc.emplace(m_trans.get().transform(m_rc.value()));
        }
        assert(m_pc.has_value());
        return m_pc.value();
      }
      case Coordinates::Reference:
      {
        if (!m_rc.has_value())
        {
          assert(m_pc.has_value());
          m_rc.emplace(m_trans.get().inverse(m_pc.value()));
        }
        assert(m_rc.has_value());
        return m_rc.value();
      }
    }

    return m_pc.value(); // Some compilers complain, so return any value
  }

  const Math::Matrix& Point::getJacobian() const
  {
    if (!m_jacobian.has_value())
    {
      assert(m_rc.has_value());
      m_jacobian.emplace(m_trans.get().jacobian(m_rc.value()));
    }
    assert(m_jacobian.has_value());
    return m_jacobian.value();
  }

  const Math::Matrix& Point::getJacobianInverse() const
  {
    if (!m_jacobianInverse.has_value())
    {
      assert(m_rc.has_value());
      const size_t rdim = Polytope::getGeometryDimension(m_polytope.get().getGeometry());
      const size_t sdim = m_polytope.get().getMesh().getSpaceDimension();
      assert(rdim <= sdim);
      if (rdim == sdim)
      {
        switch (rdim)
        {
          case 1:
          {
            Math::Matrix inv(1, 1);
            inv.coeffRef(0, 0) = 1 / getJacobian().coeff(0, 0);
            m_jacobianDeterminant.emplace(getJacobian().coeff(0, 0));
            m_jacobianInverse.emplace(std::move(inv));
            break;
          }
          case 2:
          {
            const auto& jac = getJacobian();
            const Scalar a = jac.coeff(0, 0);
            const Scalar b = jac.coeff(0, 1);
            const Scalar c = jac.coeff(1, 0);
            const Scalar d = jac.coeff(1, 1);
            const Scalar det = a * d - b * c;
            m_jacobianDeterminant.emplace(det);
            assert(det > 0);
            Math::Matrix inv(2, 2);
            inv.coeffRef(0, 0) = d / det;
            inv.coeffRef(0, 1) = -b / det;
            inv.coeffRef(1, 0) = -c / det;
            inv.coeffRef(1, 1) = a / det;
            m_jacobianInverse.emplace(std::move(inv));
            break;
          }
          case 3:
          {
            const auto& jac = getJacobian();
            const Scalar a = jac.coeff(0, 0);
            const Scalar b = jac.coeff(0, 1);
            const Scalar c = jac.coeff(0, 2);
            const Scalar d = jac.coeff(1, 0);
            const Scalar e = jac.coeff(1, 1);
            const Scalar f = jac.coeff(1, 2);
            const Scalar g = jac.coeff(2, 0);
            const Scalar h = jac.coeff(2, 1);
            const Scalar i = jac.coeff(2, 2);

            const Scalar A = e * i - f * h;
            const Scalar B = -(d * i - f * g);
            const Scalar C = d * h - e * g;
            const Scalar D = -(b * i - c * h);
            const Scalar E = a * i - c * g;
            const Scalar F = -(a * h - b * g);
            const Scalar G = b * f - c * e;
            const Scalar H = - (a * f  - c * d);
            const Scalar I = a * e - b * d;

            const Scalar det = a * A + b * B + c * C;
            m_jacobianDeterminant.emplace(det);

            assert(det > 0);
            Math::Matrix inv(3, 3);
            inv.coeffRef(0, 0) = A / det;
            inv.coeffRef(0, 1) = D / det;
            inv.coeffRef(0, 2) = G / det;
            inv.coeffRef(1, 0) = B / det;
            inv.coeffRef(1, 1) = E / det;
            inv.coeffRef(1, 2) = H / det;
            inv.coeffRef(2, 0) = C / det;
            inv.coeffRef(2, 1) = F / det;
            inv.coeffRef(2, 2) = I / det;
            m_jacobianInverse.emplace(std::move(inv));
            break;
          }
          default:
          {
            m_jacobianInverse.emplace(getJacobian().inverse());
            break;
          }
        }
      }
      else
      {
        assert(false); // Not handled yet
      }
    }
    assert(m_jacobianInverse.has_value());
    return m_jacobianInverse.value();
  }

  Scalar Point::getJacobianDeterminant() const
  {
    if (!m_jacobianDeterminant.has_value())
    {
      assert(m_rc.has_value());
      const auto& jac = getJacobian();
      const auto rows = jac.rows();
      const auto cols = jac.cols();
      if (rows == cols)
      {
        switch (rows)
        {
          case 1:
          {
            m_jacobianDeterminant.emplace(jac.coeff(0, 0));
            break;
          }
          case 2:
          {
            const Scalar a = jac.coeff(0, 0);
            const Scalar b = jac.coeff(0, 1);
            const Scalar c = jac.coeff(1, 0);
            const Scalar d = jac.coeff(1, 1);
            m_jacobianDeterminant.emplace(a * d - b * c);
            break;
          }
          case 3:
          {
            const Scalar a = jac.coeff(0, 0);
            const Scalar b = jac.coeff(0, 1);
            const Scalar c = jac.coeff(0, 2);
            const Scalar d = jac.coeff(1, 0);
            const Scalar e = jac.coeff(1, 1);
            const Scalar f = jac.coeff(1, 2);
            const Scalar g = jac.coeff(2, 0);
            const Scalar h = jac.coeff(2, 1);
            const Scalar i = jac.coeff(2, 2);
            const Scalar A = e * i - f * h;
            const Scalar B = -(d * i - f * g);
            const Scalar C = d * h - e * g;
            m_jacobianDeterminant.emplace(a * A + b * B + c * C);
            break;
          }
          default:
          {
            m_jacobianDeterminant.emplace(jac.determinant());
            break;
          }
        }
      }
      else
      {
        assert(false); // Not handled yet
      }
    }
    assert(m_jacobianDeterminant.has_value());
    return m_jacobianDeterminant.value();
  }

  Scalar Point::getDistortion() const
  {
    if (!m_distortion.has_value())
    {
      assert(m_rc.has_value());
      const auto& jac = getJacobian();
      const auto rows = jac.rows();
      const auto cols = jac.cols();
      if (rows == cols)
      {
        switch (rows)
        {
          case 1:
          {
            m_distortion.emplace(Math::abs(jac.coeff(0, 0)));
            break;
          }
          case 2:
          {
            m_distortion.emplace(Math::abs(getJacobianDeterminant()));
            break;
          }
          case 3:
          {
            m_distortion.emplace(Math::abs(getJacobianDeterminant()));
            break;
          }
          default:
          {
            m_distortion.emplace(Math::sqrt(Math::abs((jac.transpose() * jac).determinant())));
            break;
          }
        }
      }
      else
      {
        assert(false); // Not handled yet
      }
    }
    assert(m_distortion.has_value());
    return m_distortion.value();
  }

  size_t Point::getDimension(Coordinates coords) const
  {
    switch (coords)
    {
      case Coordinates::Physical:
        return m_polytope.get().getMesh().getSpaceDimension();
      case Coordinates::Reference:
        return m_polytope.get().getMesh().getDimension();
      default:
      {
        assert(false);
        return 0;
      }
    }
  }
}
