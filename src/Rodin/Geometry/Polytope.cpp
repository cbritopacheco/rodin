#include <Eigen/Cholesky>

#include "Rodin/Configure.h"

#include "Rodin/Variational/QuadratureRule.h"

#include "Mesh.h"
#include "PolytopeTransformation.h"

#include "Polytope.h"

namespace Rodin::Geometry
{
  const GeometryIndexed<Math::PointMatrix> Polytope::s_vertices =
  {
    { Polytope::Type::Point,
      Math::PointMatrix{{0}} },
    { Polytope::Type::Segment,
      Math::PointMatrix{{0, 1}} },
    { Polytope::Type::Triangle,
      Math::PointMatrix{{0, 1, 0},
                        {0, 0, 1}} },
    { Polytope::Type::Quadrilateral,
      Math::PointMatrix{{0, 1, 0, 1},
                        {0, 0, 1, 1}} },
    { Polytope::Type::Tetrahedron,
      Math::PointMatrix{{0, 1, 0, 0},
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

  Polytope::Type Polytope::getGeometry() const
  {
    return getMesh().getGeometry(getDimension(), getIndex());
  }

  VertexIterator Polytope::getVertex() const
  {
    const auto& vertices = getVertices();
    return VertexIterator(
        getMesh(), IteratorIndexGenerator(vertices.begin(), vertices.end()));
  }

  const Math::PointMatrix& Polytope::getVertices(Polytope::Type g)
  {
    return s_vertices[g];
  }

  const Array<Index>& Polytope::getVertices() const
  {
    return m_mesh.get().getConnectivity().getPolytope(getDimension(), getIndex());
  }

  PolytopeIterator Polytope::getAdjacent() const
  {
    const size_t d = getDimension();
    const auto& mesh = m_mesh.get();
    const auto& conn = mesh.getConnectivity();
    const auto& inc = conn.getIncidence(d, d);
    if (inc.size() == 0)
    {
      Alert::MemberFunctionException(*this, __func__)
        << Alert::Notation::Incidence(d, d)
        << " has not been computed and is required to use this function."
        << Alert::Raise;
    }
    const auto& adj = inc.at(getIndex());
    return PolytopeIterator(
        d, getMesh(), IteratorIndexGenerator(adj.begin(), adj.end()));
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
      case Type::Point:
      {
        return 0;
      }
      case Type::Segment:
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
      case Type::Triangle:
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
      case Type::Quadrilateral:
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
      case Type::Tetrahedron:
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

  bool Polytope::isCell() const
  {
    return getDimension() == getMesh().getDimension();
  }

  bool Polytope::isFace() const
  {
    return getDimension() == getMesh().getDimension() - 1;
  }

  bool Polytope::isVertex() const
  {
    return getDimension() == 0;
  }

  // ---- Element -----------------------------------------------------------
  Cell::Cell(Index index, const MeshBase& mesh)
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
}
