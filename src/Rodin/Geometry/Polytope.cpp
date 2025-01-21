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
    { Polytope::Type::TriangularPrism,
      Math::PointMatrix{{0, 1, 0, 0, 1, 0},
                        {0, 0, 1, 0, 0, 1},
                        {0, 0, 0, 1, 1, 1}} },
  };

  bool operator==(const Polytope& lhs, const Polytope& rhs)
  {
    bool res = true;
    res = res && (&lhs.getMesh() == &rhs.getMesh());
    res = res && (lhs.getDimension() == rhs.getDimension());
    res = res && (lhs.getIndex() == rhs.getIndex());
    return res;
  }

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

  Real Polytope::getMeasure() const
  {
    Real res = 0;
    QF::GenericPolytopeQuadrature qf(getTransformation().getJacobianOrder(), getGeometry());
    for (size_t i = 0; i < qf.getSize(); i++)
    {
      const Geometry::Point p(*this, getTransformation(), std::cref(qf.getPoint(i)));
      res += qf.getWeight(i) * p.getDistortion();
    }
    return res;
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

  Eigen::Map<const Math::SpatialVector<Real>> Vertex::getCoordinates() const
  {
    return getMesh().getVertexCoordinates(getIndex());
  }
}
