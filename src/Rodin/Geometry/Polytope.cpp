#include <Eigen/Cholesky>

#include "Rodin/Configure.h"

#include "Rodin/Variational/QuadratureRule.h"

#include "Mesh.h"
#include "PolytopeTransformation.h"

#include "Polytope.h"

namespace Rodin::Geometry
{
  Polytope::Traits::Traits(Type g)
    : m_g(g)
  {}

  bool Polytope::Traits::isSimplex()
  {
    switch (m_g)
    {
      case Type::Point:
      case Type::Segment:
      case Type::Triangle:
      case Type::Tetrahedron:
        return true;
      case Type::Quadrilateral:
      case Type::Wedge:
        return false;
    }
    assert(false);
    return false;
  }

  size_t Polytope::Traits::getDimension()
  {
    switch (m_g)
    {
      case Type::Point:
        return 0;
      case Type::Segment:
        return 1;
      case Type::Triangle:
      case Type::Quadrilateral:
        return 2;
      case Type::Tetrahedron:
      case Type::Wedge:
        return 3;
    }
    assert(false);
    return 0;
  }

  size_t Polytope::Traits::getVertexCount()
  {
    switch (m_g)
    {
      case Type::Point:
        return 1;
      case Type::Segment:
        return 2;
      case Type::Triangle:
        return 3;
      case Type::Quadrilateral:
      case Type::Tetrahedron:
        return 4;
      case Type::Wedge:
        return 6;
    }
    assert(false);
    return 0;
  }

  const Math::SpatialPoint& Polytope::Traits::getVertex(size_t i) const
  {
    switch (m_g)
    {
      case Type::Point:
      {
        static thread_local const Math::SpatialPoint s_node{{ 0 }};
        return s_node;
      }
      case Type::Segment:
      {
        static thread_local const std::vector<Math::SpatialPoint> s_nodes =
        {
          Math::SpatialPoint{{ 0 }},
          Math::SpatialPoint{{ 1 }}
        };
        return s_nodes[i];
      }
      case Type::Triangle:
      {
        static thread_local const std::vector<Math::SpatialPoint> s_nodes =
        {
          Math::SpatialPoint{{ 0, 0 }},
          Math::SpatialPoint{{ 1, 0 }},
          Math::SpatialPoint{{ 0, 1 }}
        };
        return s_nodes[i];
      }
      case Type::Quadrilateral:
      {
        static thread_local const std::vector<Math::SpatialPoint> s_nodes =
        {
          Math::SpatialPoint{{ 0, 0 }},
          Math::SpatialPoint{{ 1, 0 }},
          Math::SpatialPoint{{ 0, 1 }},
          Math::SpatialPoint{{ 1, 1 }}
        };
        return s_nodes[i];
      }
      case Type::Tetrahedron:
      {
        static thread_local const std::vector<Math::SpatialPoint> s_nodes =
        {
          Math::SpatialPoint{{ 0, 0, 0 }},
          Math::SpatialPoint{{ 1, 0, 0 }},
          Math::SpatialPoint{{ 0, 1, 0 }},
          Math::SpatialPoint{{ 0, 0, 1 }}
        };
        return s_nodes[i];
      }
      case Type::Wedge:
      {
        static thread_local const std::vector<Math::SpatialPoint> s_nodes =
        {
          Math::SpatialPoint{{ 0, 0, 0 }},
          Math::SpatialPoint{{ 1, 0, 0 }},
          Math::SpatialPoint{{ 0, 1, 0 }},
          Math::SpatialPoint{{ 0, 0, 1 }},
          Math::SpatialPoint{{ 1, 0, 1 }},
          Math::SpatialPoint{{ 0, 1, 1 }}
        };
        return s_nodes[i];
      }
    }
    assert(false);
    static thread_local const Math::SpatialPoint s_null{{}};
    return s_null;
  }

  const Polytope::Traits::HalfSpace& Polytope::Traits::getHalfSpace() const
  {
    switch (m_g)
    {
      case Type::Point:
      {
        static thread_local const HalfSpace s_hs =
        {
          Math::SpatialMatrix<Real>{},
          Math::SpatialVector<Real>{}
        };
        return s_hs;
      }
      case Type::Segment:
      {
        static thread_local const HalfSpace s_hs =
        {
          Math::SpatialMatrix<Real>{
            { 1 },
            { -1 }},
          Math::SpatialVector<Real>{{ 1, 0 }}
        };
        return s_hs;
      }
      case Type::Triangle:
      {
        static thread_local const HalfSpace s_hs =
        {
          Math::SpatialMatrix<Real>{
            { -1,  0 },
            {  0, -1 },
            {  1,  1 }},
          Math::SpatialVector<Real>{{ 0, 0, 1 }}
        };
        return s_hs;
      }
      case Type::Quadrilateral:
      {
        static thread_local const HalfSpace s_hs =
        {
          Math::SpatialMatrix<Real>{
            {  0, -1 },
            {  1,  0 },
            {  0,  1 },
            { -1,  0 }},
          Math::SpatialVector<Real>{{ 0, 1, 1, 0 }}
        };
        return s_hs;
      }
      case Type::Tetrahedron:
      {
        static thread_local const HalfSpace s_hs =
        {
          Math::SpatialMatrix<Real>{
            {  0, -1,  0 },
            {  0,  0, -1 },
            { -1,  0,  0 },
            {  1,  1,  1 }},
          Math::SpatialVector<Real>{{ 0, 0, 0, 1 }}
        };
        return s_hs;
      }
      case Type::Wedge:
      {
        static thread_local const HalfSpace s_hs =
        {
          Math::SpatialMatrix<Real>{
            {  0,  0, -1 },
            {  0, -1,  0 },
            {  1,  1,  0 },
            { -1,  0,  0 },
            {  0,  0,  1 }},
          Math::SpatialVector<Real>{{ 0, 0, 1, 0, 1 }}
        };
        return s_hs;
      }
    }
    assert(false);
    static thread_local const HalfSpace s_null;
    return s_null;
  }

  std::ostream& operator<<(std::ostream& os, const Polytope::Type& p)
  {
    switch (p)
    {
      case Polytope::Type::Point:
      {
        os << "Point";
        break;
      }
      case Polytope::Type::Segment:
      {
        os << "Segment";
        break;
      }
      case Polytope::Type::Triangle:
      {
        os << "Triangle";
        break;
      }
      case Polytope::Type::Quadrilateral:
      {
        os << "Quadrilateral";
        break;
      }
      case Polytope::Type::Tetrahedron:
      {
        os << "Tetrahedron";
        break;
      }
      case Polytope::Type::Wedge:
      {
        os << "Wedge";
        break;
      }
    }
    return os;
  }

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
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, d);
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
      const Geometry::Point p(*this, qf.getPoint(i));
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

  Eigen::Map<const Math::SpatialPoint> Vertex::getCoordinates() const
  {
    return getMesh().getVertexCoordinates(getIndex());
  }
}
