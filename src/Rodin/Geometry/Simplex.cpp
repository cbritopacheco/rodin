#include <Eigen/Cholesky>

#include "Rodin/Variational/QuadratureRule.h"

#include "Mesh.h"
#include "SimplexTransformation.h"

#include "Simplex.h"

namespace Rodin::Geometry
{

  const GeometryIndexed<Math::Matrix> Polytope::s_vertices =
  {
    { Polytope::Geometry::Point,
      Math::Matrix{{0}} },
    { Polytope::Geometry::Segment,
      Math::Matrix{{0, 1}} },
    { Polytope::Geometry::Triangle,
      Math::Matrix{{0, 1, 0},
                   {0, 0, 1}} },
    { Polytope::Geometry::Quadrilateral,
      Math::Matrix{{0, 1, 0, 1},
                   {0, 0, 1, 1}} },
    { Polytope::Geometry::Tetrahedron,
      Math::Matrix{{0, 1, 0, 0},
                   {0, 0, 1, 0},
                   {0, 0, 0, 1}} },
  };

  bool operator<(const Polytope& lhs, const Polytope& rhs)
  {
    return lhs.getIndex() < rhs.getIndex();
  }

  // ---- Polytope -----------------------------------------------------------
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

  const Math::Matrix& Polytope::getVertices(Polytope::Geometry g)
  {
    return s_vertices[g];
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

  Scalar Polytope::getMeasure() const
  {
    const auto& mesh = getMesh();
    switch (getGeometry())
    {
      case Geometry::Point:
      {
        return 0;
      }
      case Geometry::Segment:
      {
        const auto& vertices = getVertices();
        const auto& a = mesh.getVertexCoordinates(vertices.coeff(0));
        const auto& b = mesh.getVertexCoordinates(vertices.coeff(1));
        const Scalar x0 = a.x();
        const Scalar y0 = a.y();
        const Scalar x1 = b.x();
        const Scalar y1 = b.y();
        return Math::sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
      }
      case Geometry::Triangle:
      {
        const auto& vertices = getVertices();
        const auto& a = mesh.getVertexCoordinates(vertices.coeff(0));
        const auto& b = mesh.getVertexCoordinates(vertices.coeff(1));
        const auto& c = mesh.getVertexCoordinates(vertices.coeff(2));
        const Scalar x0 = a.x();
        const Scalar y0 = a.y();
        const Scalar x1 = b.x();
        const Scalar y1 = b.y();
        const Scalar x2 = c.x();
        const Scalar y2 = c.y();
        return (1.0 / 2.0) * Math::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
      }
      case Geometry::Quadrilateral:
      {
        const auto& vertices = getVertices();
        const auto& a = mesh.getVertexCoordinates(vertices.coeff(0));
        const auto& b = mesh.getVertexCoordinates(vertices.coeff(1));
        const auto& c = mesh.getVertexCoordinates(vertices.coeff(3));
        const auto& d = mesh.getVertexCoordinates(vertices.coeff(2));
        const Scalar x0 = a.x();
        const Scalar y0 = a.y();
        const Scalar x1 = b.x();
        const Scalar y1 = b.y();
        const Scalar x2 = c.x();
        const Scalar y2 = c.y();
        const Scalar x3 = d.x();
        const Scalar y3 = d.y();
        return (1.0 / 2.0) * Math::abs((x0 * y1 - x1 * y0) + (x1 * y2 - x2 * y1) + (x2 * y3 - x3 * y2) + (x3 * y0 - x0 * y3));
      }
      case Geometry::Tetrahedron:
      {
        Eigen::Matrix<Scalar, 4, 4> pm;
        const auto& vertices = getVertices();
        const auto& a = mesh.getVertexCoordinates(vertices.coeff(0));
        const auto& b = mesh.getVertexCoordinates(vertices.coeff(1));
        const auto& c = mesh.getVertexCoordinates(vertices.coeff(3));
        const auto& d = mesh.getVertexCoordinates(vertices.coeff(2));
        pm << a.x(), a.y(), a.z(), 1,
              b.x(), b.y(), b.z(), 1,
              c.x(), c.y(), c.z(), 1,
              d.x(), d.y(), d.z(), 1;
        return (1.0 / 6.0) * pm.determinant();
      }
    }

    assert(false);
    return NAN;
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

  Eigen::Map<const Math::SpatialVector> Vertex::getCoordinates() const
  {
    return getMesh().getVertexCoordinates(getIndex());
  }

  // ---- Point --------------------------------------------------------------
  PointBase::PointBase(const Polytope& simplex, const PolytopeTransformation& trans)
    : m_polytope(simplex), m_trans(trans)
  {}

  PointBase::PointBase(const Polytope& simplex, const PolytopeTransformation& trans, const Math::SpatialVector& pc)
    : m_polytope(simplex), m_trans(trans), m_pc(pc)
  {}

  const Math::SpatialVector& PointBase::getCoordinates(Coordinates coords) const
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

  const Math::SpatialVector& PointBase::getPhysicalCoordinates() const
  {
    if (!m_pc.has_value())
      m_pc.emplace(m_trans.get().transform(getReferenceCoordinates()));
    assert(m_pc.has_value());
    return m_pc.value();
  }

  const Math::SpatialMatrix& PointBase::getJacobian() const
  {
    if (!m_jacobian.has_value())
      m_jacobian.emplace(m_trans.get().jacobian(getReferenceCoordinates()));
    assert(m_jacobian.has_value());
    return m_jacobian.value();
  }

  const Math::SpatialMatrix& PointBase::getJacobianInverse() const
  {
    if (!m_jacobianInverse.has_value())
    {
      const size_t rdim = Polytope::getGeometryDimension(m_polytope.get().getGeometry());
      const size_t sdim = m_polytope.get().getMesh().getSpaceDimension();
      assert(rdim <= sdim);
      if (rdim == sdim)
      {
        return m_jacobianInverse.emplace(getJacobian().inverse());
        switch (rdim)
        {
          case 1:
          {
            Math::Matrix inv(1, 1);
            inv.coeffRef(0, 0) = 1 / getJacobian().coeff(0, 0);
            m_jacobianDeterminant.emplace(getJacobian().coeff(0, 0));
            return m_jacobianInverse.emplace(std::move(inv));
          }
          case 2:
          {
            const auto& jac = getJacobian();
            const Scalar a = jac.coeff(0, 0);
            const Scalar b = jac.coeff(0, 1);
            const Scalar c = jac.coeff(1, 0);
            const Scalar d = jac.coeff(1, 1);
            const Scalar det = a * d - b * c;
            assert(det != 0);
            m_jacobianDeterminant.emplace(det);
            Math::Matrix inv(2, 2);
            inv.coeffRef(0, 0) = d / det;
            inv.coeffRef(0, 1) = -b / det;
            inv.coeffRef(1, 0) = -c / det;
            inv.coeffRef(1, 1) = a / det;
            return m_jacobianInverse.emplace(std::move(inv));
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

            assert(det != 0);
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
            return m_jacobianInverse.emplace(std::move(inv));
          }
          default:
          {
            return m_jacobianInverse.emplace(getJacobian().inverse());
          }
        }
      }
      else
      {
        const auto& jac = getJacobian();
        return m_jacobianInverse.emplace(jac.completeOrthogonalDecomposition().pseudoInverse());
      }
    }
    assert(m_jacobianInverse.has_value());
    return m_jacobianInverse.value();
  }

  Scalar PointBase::getJacobianDeterminant() const
  {
    if (!m_jacobianDeterminant.has_value())
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
            return m_jacobianDeterminant.emplace(jac.coeff(0, 0));
          }
          case 2:
          {
            const Scalar a = jac.coeff(0, 0);
            const Scalar b = jac.coeff(0, 1);
            const Scalar c = jac.coeff(1, 0);
            const Scalar d = jac.coeff(1, 1);
            return m_jacobianDeterminant.emplace(a * d - b * c);
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
        assert(false); // Not handled yet
      }
    }
    assert(m_jacobianDeterminant.has_value());
    return m_jacobianDeterminant.value();
  }

  Scalar PointBase::getDistortion() const
  {
    if (!m_distortion.has_value())
    {
      const auto& jac = getJacobian();
      const auto rows = jac.rows();
      const auto cols = jac.cols();
      if (rows == cols)
      {
        return m_distortion.emplace(getJacobianDeterminant());
      }
      else
      {
        if (jac.rows() == 2 && jac.cols() == 1)
        {
          return m_distortion.emplace(std::sqrt(jac.coeff(0) * jac.coeff(0) + jac.coeff(1) * jac.coeff(1)));
        }
        else
        {
          return m_distortion.emplace(Math::sqrt(Math::abs((jac.transpose() * jac).determinant())));
        }
      }
    }
    assert(m_distortion.has_value());
    return m_distortion.value();
  }

  size_t PointBase::getDimension(Coordinates coords) const
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
