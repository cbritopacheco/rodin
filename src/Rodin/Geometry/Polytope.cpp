#include <Eigen/Cholesky>

#include "Rodin/Configure.h"

#include "Rodin/QF/GenericPolytopeQuadrature.h"

#include "Mesh.h"
#include "PolytopeTransformation.h"

#include "Polytope.h"

namespace Rodin::Geometry
{
  Polytope::Traits::Traits(Type g)
    : m_g(g)
  {}

  bool Polytope::Traits::isSimplex() const
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

  size_t Polytope::Traits::getDimension() const
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

  size_t Polytope::Traits::getVertexCount() const
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
    assert(i < getVertexCount());
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
          Math::SpatialPoint{{ 1, 1 }},
          Math::SpatialPoint{{ 0, 1 }}
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
          Math::Matrix<Real>{},
          Math::Vector<Real>{}
        };
        return s_hs;
      }
      case Type::Segment:
      {
        static thread_local const HalfSpace s_hs =
        {
          // local 0: x=0  -> -x <= 0
          // local 1: x=1  ->  x <= 1
          Math::Matrix<Real>{{ -1 }, { 1 }},
          Math::Vector<Real>{{ 0, 1 }}
        };
        return s_hs;
      }
      case Type::Triangle:
      {
        static thread_local const HalfSpace s_hs =
        {
          // local 0: y=0                ->  -y <= 0
          // local 1: x+y=1              ->  (x+y)/√2 <= 1/√2
          // local 2: x=0                ->  -x <= 0
          Math::Matrix<Real>{
            {  0, -1 },
            {  1 / std::sqrt(2.0),  1 / std::sqrt(2.0) },
            { -1,  0 }
          },
          Math::Vector<Real>{{ 0, 1 / std::sqrt(2.0), 0 }}
        };
        return s_hs;
      }
      case Type::Quadrilateral:
      {
        static thread_local const HalfSpace s_hs =
        {
          // local 0: y=0   -> -y <= 0
          // local 1: x=1   ->  x <= 1
          // local 2: y=1   ->  y <= 1
          // local 3: x=0   -> -x <= 0
          Math::Matrix<Real>{
            {  0, -1 },
            {  1,  0 },
            {  0,  1 },
            { -1,  0 }
          },
          Math::Vector<Real>{{ 0, 1, 1, 0 }}
        };
        return s_hs;
      }
      case Type::Tetrahedron:
      {
        static thread_local const HalfSpace s_hs =
        {
          // local 0: y=0   -> -y <= 0
          // local 1: z=0   -> -z <= 0
          // local 2: x=0   -> -x <= 0
          // local 3: x+y+z=1 -> (x+y+z)/√3 <= 1/√3
          Math::Matrix<Real>{
            {  0, -1,  0 },
            {  0,  0, -1 },
            { -1,  0,  0 },
            {  1 / std::sqrt(3.0),  1 / std::sqrt(3.0),  1 / std::sqrt(3.0) }
          },
          Math::Vector<Real>{{ 0, 0, 0, 1 / std::sqrt(3.0) }}
        };
        return s_hs;
      }
      case Type::Wedge:
      {
        static thread_local const HalfSpace s_hs =
        {
          // local 0: z=0                ->  -z <= 0
          // local 1: y=0                ->  -y <= 0
          // local 2: x+y=1              ->  (x+y)/√2 <= 1/√2
          // local 3: x=0                ->  -x <= 0
          // local 4: z=1                ->   z <= 1
          Math::Matrix<Real>{
            {  0,  0, -1 },
            {  0, -1,  0 },
            {  1 / std::sqrt(2.0),  1 / std::sqrt(2.0), 0 },
            { -1,  0,  0 },
            {  0,  0,  1 }
          },
          Math::Vector<Real>{{ 0, 0, 1 / std::sqrt(2.0), 0, 1 }}
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
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, d);
    const auto& inc = conn.getIncidence(d, d);
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

  Polytope::Project::Project(Type g)
    : m_g(g)
  {}

  void Polytope::Project::cell(Math::SpatialPoint& out, const Math::SpatialPoint& rc) const
  {
    assert(rc.size() >= 0);
    assert(static_cast<size_t>(rc.size()) == Polytope::Traits(m_g).getDimension());

    out.resize(rc.size());

    switch (m_g)
    {
      case Type::Point:
      {
        return;
      }
      case Type::Segment:
      {
        out[0] = std::clamp(rc[0], Real(0), Real(1));
        return;
      }
      case Type::Triangle:
      {
        Real x = std::max(rc[0], Real(0));
        Real y = std::max(rc[1], Real(0));
        const Real s = x + y;
        if (s <= Real(1)) { out[0] = x; out[1] = y; return; }
        const Real t = 0.5 * (s - Real(1));
        x = std::max(x - t, Real(0));
        y = std::max(y - t, Real(0));
        out[0] = std::clamp(x, Real(0), Real(1));
        out[1] = std::clamp(y, Real(0), Real(1));
        return;
      }
      case Type::Quadrilateral:
      {
        out[0] = std::clamp(rc[0], Real(0), Real(1));
        out[1] = std::clamp(rc[1], Real(0), Real(1));
        return;
      }
      case Type::Tetrahedron:
      {
        Real x = std::max(rc[0], Real(0));
        Real y = std::max(rc[1], Real(0));
        Real z = std::max(rc[2], Real(0));
        const Real s = x + y + z;
        if (s <= Real(1)) { out[0] = x; out[1] = y; out[2] = z; return; }

        Real u0 = x, u1 = y, u2 = z;
        if (u1 > u0) std::swap(u0, u1);
        if (u2 > u1) std::swap(u1, u2);
        if (u1 > u0) std::swap(u0, u1);
        const Real s1 = u0, s2 = u0 + u1, s3 = s2 + u2;
        Real theta = s1 - Real(1);
        if (u1 - (s2 - Real(1)) / 2.0 > Real(0)) theta = (s2 - Real(1)) / 2.0;
        if (u2 - (s3 - Real(1)) / 3.0 > Real(0)) theta = (s3 - Real(1)) / 3.0;

        out[0] = std::max(x - theta, Real(0));
        out[1] = std::max(y - theta, Real(0));
        out[2] = std::max(z - theta, Real(0));
        return;
      }
      case Type::Wedge:
      {
        Real x = std::max(rc[0], Real(0));
        Real y = std::max(rc[1], Real(0));
        const Real s = x + y;
        if (s > Real(1))
        {
          const Real t = 0.5 * (s - Real(1));
          x = std::max(x - t, Real(0));
          y = std::max(y - t, Real(0));
        }
        const Real z = std::clamp(rc[2], Real(0), Real(1));
        out[0] = x; out[1] = y; out[2] = z;
        return;
      }
    }
    assert(false);
  }

  void Polytope::Project::boundary(Math::SpatialPoint& out, const Math::SpatialPoint& rc) const
  {
    assert(rc.size() >= 0);
    assert(static_cast<size_t>(rc.size()) == Polytope::Traits(m_g).getDimension());

    out.resize(rc.size());

    switch (m_g)
    {
      case Type::Point:
      {
        return;
      }
      case Type::Segment:
      {
        out[0] = (rc[0] <= 0.5) ? 0 : 1;
        return;
      }
      case Type::Triangle:
      {
        const Real e0x = std::clamp(rc[0], Real(0), Real(1));
        const Real e0y = Real(0);

        const Real dx = rc[0] - Real(1);
        const Real dy = rc[1];
        Real t = 0.5 * (-dx + dy);
        t = std::clamp(t, Real(0), Real(1));
        const Real e1x = Real(1) - t;
        const Real e1y = t;

        const Real e2x = Real(0);
        const Real e2y = std::clamp(rc[1], Real(0), Real(1));

        const Real d0 = (rc[0] - e0x) * (rc[0] - e0x) + (rc[1] - e0y) * (rc[1] - e0y);
        const Real d1 = (rc[0] - e1x) * (rc[0] - e1x) + (rc[1] - e1y) * (rc[1] - e1y);
        const Real d2 = (rc[0] - e2x) * (rc[0] - e2x) + (rc[1] - e2y) * (rc[1] - e2y);

        if (d0 <= d1 && d0 <= d2) { out[0] = e0x; out[1] = e0y; return; }
        if (d1 <= d2)             { out[0] = e1x; out[1] = e1y; return; }
        out[0] = e2x; out[1] = e2y; return;
      }
      case Type::Quadrilateral:
      {
        const Real e0x = std::clamp(rc[0], Real(0), Real(1)), e0y = Real(0);
        const Real e1x = Real(1), e1y = std::clamp(rc[1], Real(0), Real(1));
        const Real e2x = std::clamp(rc[0], Real(0), Real(1)), e2y = Real(1);
        const Real e3x = Real(0), e3y = std::clamp(rc[1], Real(0), Real(1));

        const Real d0 = (rc[0] - e0x)*(rc[0] - e0x) + (rc[1] - e0y)*(rc[1] - e0y);
        const Real d1 = (rc[0] - e1x)*(rc[0] - e1x) + (rc[1] - e1y)*(rc[1] - e1y);
        const Real d2 = (rc[0] - e2x)*(rc[0] - e2x) + (rc[1] - e2y)*(rc[1] - e2y);
        const Real d3 = (rc[0] - e3x)*(rc[0] - e3x) + (rc[1] - e3y)*(rc[1] - e3y);

        if (d0 <= d1 && d0 <= d2 && d0 <= d3) { out[0] = e0x; out[1] = e0y; return; }
        if (d1 <= d2 && d1 <= d3)             { out[0] = e1x; out[1] = e1y; return; }
        if (d2 <= d3)                         { out[0] = e2x; out[1] = e2y; return; }
        out[0] = e3x; out[1] = e3y; return;
      }
      case Type::Tetrahedron:
      {
        Real f0x = std::max(rc[0], Real(0));
        Real f0z = std::max(rc[2], Real(0));
        Real s = f0x + f0z;
        if (s > Real(1))
        {
          const Real t0 = 0.5 * (s - Real(1));
          f0x = std::max(f0x - t0, Real(0));
          f0z = std::max(f0z - t0, Real(0));
        }
        const Real f0y = Real(0);

        Real f1x = std::max(rc[0], Real(0));
        Real f1y = std::max(rc[1], Real(0));
        s = f1x + f1y;
        if (s > Real(1))
        {
          const Real t1 = 0.5 * (s - Real(1));
          f1x = std::max(f1x - t1, Real(0));
          f1y = std::max(f1y - t1, Real(0));
        }
        const Real f1z = Real(0);

        Real f2y = std::max(rc[1], Real(0));
        Real f2z = std::max(rc[2], Real(0));
        s = f2y + f2z;
        if (s > Real(1))
        {
          const Real t2 = 0.5 * (s - Real(1));
          f2y = std::max(f2y - t2, Real(0));
          f2z = std::max(f2z - t2, Real(0));
        }
        const Real f2x = Real(0);

        Real a = rc[0], b = rc[1], c = rc[2];
        {
          Real u0 = a, u1 = b, u2 = c; int i0 = 0, i1 = 1, i2 = 2;
          if (u1 > u0) { std::swap(u0, u1); std::swap(i0, i1); }
          if (u2 > u1) { std::swap(u1, u2); std::swap(i1, i2); }
          if (u1 > u0) { std::swap(u0, u1); std::swap(i0, i1); }
          const Real s1 = u0, s2 = u0 + u1, s3 = s2 + u2;
          Real theta = s1 - Real(1);
          if (u1 - (s2 - Real(1)) / 2.0 > Real(0)) theta = (s2 - Real(1)) / 2.0;
          if (u2 - (s3 - Real(1)) / 3.0 > Real(0)) theta = (s3 - Real(1)) / 3.0;
          Real w0 = std::max(u0 - theta, Real(0));
          Real w1 = std::max(u1 - theta, Real(0));
          Real w2 = std::max(u2 - theta, Real(0));
          Real r[3]; r[i0] = w0; r[i1] = w1; r[i2] = w2;
          a = r[0]; b = r[1]; c = r[2];
        }

        const Real d0 = (rc[0]-f0x)*(rc[0]-f0x) + (rc[1]-f0y)*(rc[1]-f0y) + (rc[2]-f0z)*(rc[2]-f0z);
        const Real d1 = (rc[0]-f1x)*(rc[0]-f1x) + (rc[1]-f1y)*(rc[1]-f1y) + (rc[2]-f1z)*(rc[2]-f1z);
        const Real d2 = (rc[0]-f2x)*(rc[0]-f2x) + (rc[1]-f2y)*(rc[1]-f2y) + (rc[2]-f2z)*(rc[2]-f2z);
        const Real d3 = (rc[0]-a  )*(rc[0]-a  ) + (rc[1]-b  )*(rc[1]-b  ) + (rc[2]-c  )*(rc[2]-c  );

        if (d0 <= d1 && d0 <= d2 && d0 <= d3) { out[0] = f0x; out[1] = f0y; out[2] = f0z; return; }
        if (d1 <= d2 && d1 <= d3)             { out[0] = f1x; out[1] = f1y; out[2] = f1z; return; }
        if (d2 <= d3)                         { out[0] = f2x; out[1] = f2y; out[2] = f2z; return; }
        out[0] = a; out[1] = b; out[2] = c; return;
      }
      case Type::Wedge:
      {
        Real g0x = std::max(rc[0], Real(0));
        Real g0y = std::max(rc[1], Real(0));
        Real s = g0x + g0y;
        if (s > Real(1))
        {
          const Real t0 = 0.5 * (s - Real(1));
          g0x = std::max(g0x - t0, Real(0));
          g0y = std::max(g0y - t0, Real(0));
        }
        const Real g0z = Real(0);

        const Real g1x = std::clamp(rc[0], Real(0), Real(1));
        const Real g1y = Real(0);
        const Real g1z = std::clamp(rc[2], Real(0), Real(1));

        const Real dx = rc[0] - Real(1);
        const Real dy = rc[1];
        Real t = 0.5 * (-dx + dy);
        t = std::clamp(t, Real(0), Real(1));
        const Real g2x = Real(1) - t;
        const Real g2y = t;
        const Real g2z = std::clamp(rc[2], Real(0), Real(1));

        const Real g3x = Real(0);
        const Real g3y = std::clamp(rc[1], Real(0), Real(1));
        const Real g3z = std::clamp(rc[2], Real(0), Real(1));

        Real g4x = std::max(rc[0], Real(0));
        Real g4y = std::max(rc[1], Real(0));
        s = g4x + g4y;
        if (s > Real(1))
        {
          const Real t1 = 0.5 * (s - Real(1));
          g4x = std::max(g4x - t1, Real(0));
          g4y = std::max(g4y - t1, Real(0));
        }
        const Real g4z = Real(1);

        const Real d0 = (rc[0] - g0x) * (rc[0] - g0x) + (rc[1] - g0y) * (rc[1] - g0y) + (rc[2] - g0z) * (rc[2] - g0z);
        const Real d1 = (rc[0] - g1x) * (rc[0] - g1x) + (rc[1] - g1y) * (rc[1] - g1y) + (rc[2] - g1z) * (rc[2] - g1z);
        const Real d2 = (rc[0] - g2x) * (rc[0] - g2x) + (rc[1] - g2y) * (rc[1] - g2y) + (rc[2] - g2z) * (rc[2] - g2z);
        const Real d3 = (rc[0] - g3x) * (rc[0] - g3x) + (rc[1] - g3y) * (rc[1] - g3y) + (rc[2] - g3z) * (rc[2] - g3z);
        const Real d4 = (rc[0] - g4x) * (rc[0] - g4x) + (rc[1] - g4y) * (rc[1] - g4y) + (rc[2] - g4z) * (rc[2] - g4z);

        if (d0 <= d1 && d0 <= d2 && d0 <= d3 && d0 <= d4) { out[0] = g0x; out[1] = g0y; out[2] = g0z; return; }
        if (d1 <= d2 && d1 <= d3 && d1 <= d4)             { out[0] = g1x; out[1] = g1y; out[2] = g1z; return; }
        if (d2 <= d3 && d2 <= d4)                         { out[0] = g2x; out[1] = g2y; out[2] = g2z; return; }
        if (d3 <= d4)                                     { out[0] = g3x; out[1] = g3y; out[2] = g3z; return; }
        out[0]=g4x; out[1]=g4y; out[2]=g4z; return;
      }
    }
    assert(false);
  }

  void Polytope::Project::face(size_t local, Math::SpatialPoint& out, const Math::SpatialPoint& rc) const
  {
    assert(rc.size() >= 0);
    assert(static_cast<size_t>(rc.size()) == Polytope::Traits(m_g).getDimension());

    out.resize(rc.size());

    switch (m_g)
    {
      case Type::Point:
      {
        assert(local == 0);
        return;
      }
      case Type::Segment:
      {
        assert(local < 2);
        out[0] = local ? 1 : 0;
        return;
      }
      case Type::Triangle:
      {
        assert(local < 3);
        if (local == 0)
        {
          const Real x = std::clamp(rc[0], Real(0), Real(1));
          out[0] = x;
          out[1] = Real(0);
          return;
        }

        if (local == 1)
        {
          const Real dx = rc[0] - Real(1), dy = rc[1];
          Real t = 0.5 * (-dx + dy);
          t = std::clamp(t, Real(0), Real(1));
          out[0] = Real(1) - t;
          out[1] = t;
          return;
        }

        const Real y = std::clamp(rc[1], Real(0), Real(1));
        out[0] = Real(0);
        out[1] = Real(1) - y;
        return;
      }
      case Type::Quadrilateral:
      {
        assert(local < 4);
        if (local == 0) {                    // y = 0,  (0,0) → (1,0)
          const Real x = std::clamp(rc[0], Real(0), Real(1));
          out[0] = x; out[1] = Real(0); return;
        }
        if (local == 1) {                    // x = 1,  (1,0) → (1,1)
          const Real y = std::clamp(rc[1], Real(0), Real(1));
          out[0] = Real(1); out[1] = y; return;
        }
        if (local == 2) {                    // y = 1,  (1,1) → (0,1)
          const Real x = std::clamp(rc[0], Real(0), Real(1));
          out[0] = Real(1) - x; out[1] = Real(1); return;
        }
        {                                    // x = 0,  (0,1) → (0,0)
          const Real y = std::clamp(rc[1], Real(0), Real(1));
          out[0] = Real(0); out[1] = Real(1) - y; return;
        }
      }
      case Type::Tetrahedron:
      {
        assert(local < 4);
        if (local == 0)
        {
          Real x = std::max(rc[0], Real(0));
          Real z = std::max(rc[2], Real(0));
          const Real s = x + z;
          if (s > Real(1))
          {
            const Real t = 0.5 * (s - Real(1));
            x = std::max(x - t, Real(0));
            z = std::max(z - t, Real(0));
          }
          out[0] = x; out[1] = Real(0); out[2] = z; return;
        }
        if (local == 1)
        {
          Real x = std::max(rc[0], Real(0));
          Real y = std::max(rc[1], Real(0));
          const Real s = x + y;
          if (s > Real(1))
          {
            const Real t = 0.5 * (s - Real(1));
            x = std::max(x - t, Real(0));
            y = std::max(y - t, Real(0));
          }
          out[0] = x; out[1] = y; out[2] = Real(0); return;
        }
        if (local == 2)
        {
          Real y = std::max(rc[1], Real(0));
          Real z = std::max(rc[2], Real(0));
          const Real s = y + z;
          if (s > Real(1))
          {
            const Real t = 0.5 * (s - Real(1));
            y = std::max(y - t, Real(0));
            z = std::max(z - t, Real(0));
          }
          out[0] = Real(0); out[1] = y; out[2] = z; return;
        }
        {
          Real a = rc[0], b = rc[1], c = rc[2];
          Real u0 = a, u1 = b, u2 = c; int i0 = 0, i1 = 1, i2 = 2;
          if (u1 > u0) { std::swap(u0, u1); std::swap(i0, i1); }
          if (u2 > u1) { std::swap(u1, u2); std::swap(i1, i2); }
          if (u1 > u0) { std::swap(u0, u1); std::swap(i0, i1); }
          const Real s1 = u0, s2 = u0 + u1, s3 = s2 + u2;
          Real theta = s1 - Real(1);
          if (u1 - (s2 - Real(1)) / 2.0 > Real(0)) theta = (s2 - Real(1)) / 2.0;
          if (u2 - (s3 - Real(1)) / 3.0 > Real(0)) theta = (s3 - Real(1)) / 3.0;
          Real w0 = std::max(u0 - theta, Real(0));
          Real w1 = std::max(u1 - theta, Real(0));
          Real w2 = std::max(u2 - theta, Real(0));
          Real r[3]; r[i0] = w0; r[i1] = w1; r[i2] = w2;
          out[0] = r[0]; out[1] = r[1]; out[2] = r[2]; return;
        }
      }
      case Type::Wedge:
      {
        assert(local < 5);
        if (local == 0)
        {
          Real x = std::max(rc[0], Real(0));
          Real y = std::max(rc[1], Real(0));
          const Real s = x + y;
          if (s > Real(1))
          {
            const Real t = 0.5 * (s - Real(1));
            x = std::max(x - t, Real(0));
            y = std::max(y - t, Real(0));
          }
          out[0] = x; out[1] = y; out[2] = Real(0); return;
        }
        if (local == 1)
        {
          const Real x = std::clamp(rc[0], Real(0), Real(1));
          const Real z = std::clamp(rc[2], Real(0), Real(1));
          out[0] = x; out[1] = Real(0); out[2] = z; return;
        }
        if (local == 2)
        {
          const Real dx = rc[0] - Real(1), dy = rc[1];
          Real t = 0.5 * (-dx + dy);
          t = std::clamp(t, Real(0), Real(1));
          const Real x = Real(1) - t, y = t;
          const Real z = std::clamp(rc[2], Real(0), Real(1));
          out[0] = x; out[1] = y; out[2] = z; return;
        }
        if (local == 3)
        {
          const Real y = std::clamp(rc[1], Real(0), Real(1));
          const Real z = std::clamp(rc[2], Real(0), Real(1));
          out[0] = Real(0); out[1] = y; out[2] = z; return;
        }
        {
          Real x = std::max(rc[0], Real(0));
          Real y = std::max(rc[1], Real(0));
          const Real s = x + y;
          if (s > Real(1))
          {
            const Real t = 0.5 * (s - Real(1));
            x = std::max(x - t, Real(0));
            y = std::max(y - t, Real(0));
          }
          out[0] = x; out[1] = y; out[2] = Real(1); return;
        }
      }
    }
    assert(false);
  }
}
