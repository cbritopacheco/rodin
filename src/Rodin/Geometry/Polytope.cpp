#include <Eigen/Cholesky>

#include "Rodin/Configure.h"

#include "Rodin/Math/SpatialVector.h"
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
      case Type::Hexahedron:
      case Type::Wedge:
        return false;
    }
    assert(false);
    return false;
  }

  bool Polytope::Traits::isTensorProduct() const
  {
    switch (m_g)
    {
      case Type::Point:
      case Type::Segment:
      case Type::Quadrilateral:
      case Type::Hexahedron:
      case Type::Wedge:
        return true;
      case Type::Triangle:
      case Type::Tetrahedron:
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
      case Type::Hexahedron:
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
      case Type::Hexahedron:
        return 8;
      case Type::Wedge:
        return 6;
    }

    assert(false);
    return 0;
  }

  Math::SpatialPoint Polytope::Traits::getCentroid() const
  {
    switch (m_g)
    {
      case Type::Point:
      {
        return Math::SpatialPoint{ 0 };
      }
      case Type::Segment:
      {
        return Math::SpatialPoint{ 0.5 };
      }
      case Type::Triangle:
      {
        return Math::SpatialPoint{{ 1.0 / 3.0, 1.0 / 3.0 }};
      }
      case Type::Quadrilateral:
      {
        return Math::SpatialPoint{{ 0.5, 0.5 }};
      }
      case Type::Tetrahedron:
      {
        return Math::SpatialPoint{{ 0.25, 0.25, 0.25 }};
      }
      case Type::Hexahedron:
      {
        return Math::SpatialPoint{{ 0.5, 0.5, 0.5 }};
      }
      case Type::Wedge:
      {
        return Math::SpatialPoint{{ 1.0 / 3.0, 1.0 / 3.0, 0.5 }};
      }
    }
    return Math::SpatialPoint{};
  }

  const Math::SpatialPoint& Polytope::Traits::getVertex(size_t i) const
  {
    assert(i < getVertexCount());
    switch (m_g)
    {
      case Type::Point:
      {
        static thread_local const Math::SpatialPoint s_node{ 0 };
        return s_node;
      }
      case Type::Segment:
      {
        static thread_local const std::vector<Math::SpatialPoint> s_nodes =
        {
          Math::SpatialPoint{ 0 },
          Math::SpatialPoint{ 1 }
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
      case Type::Hexahedron:
      {
        static thread_local const std::vector<Math::SpatialPoint> s_nodes =
        {
          Math::SpatialPoint{{ 0, 0, 0 }}, // 0: (x=0,y=0,z=0)
          Math::SpatialPoint{{ 1, 0, 0 }}, // 1: (1,0,0)
          Math::SpatialPoint{{ 1, 1, 0 }}, // 2: (1,1,0)
          Math::SpatialPoint{{ 0, 1, 0 }}, // 3: (0,1,0)
          Math::SpatialPoint{{ 0, 0, 1 }}, // 4: (0,0,1)
          Math::SpatialPoint{{ 1, 0, 1 }}, // 5: (1,0,1)
          Math::SpatialPoint{{ 1, 1, 1 }}, // 6: (1,1,1)
          Math::SpatialPoint{{ 0, 1, 1 }}  // 7: (0,1,1)
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
      case Type::Hexahedron:
      {
        static thread_local const HalfSpace s_hs =
        {
          // local 0: z=0   -> -z <= 0
          // local 1: y=0   -> -y <= 0
          // local 2: x=1   ->  x <= 1
          // local 3: y=1   ->  y <= 1
          // local 4: x=0   -> -x <= 0
          // local 5: z=1   ->  z <= 1
          Math::Matrix<Real>{
            {  0,  0, -1 },
            {  0, -1,  0 },
            {  1,  0,  0 },
            {  0,  1,  0 },
            { -1,  0,  0 },
            {  0,  0,  1 }
          },
          Math::Vector<Real>{{ 0, 0, 1, 1, 0, 1 }}
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
          // Must match Connectivity::getSubPolytopes tetra face order:
          //
          // face 0: (1,2,3)  -> x+y+z=1   ->  (x+y+z)/sqrt(3) <= 1/sqrt(3)
          // face 1: (0,3,2)  -> x=0       ->  -x <= 0
          // face 2: (0,1,3)  -> y=0       ->  -y <= 0
          // face 3: (0,2,1)  -> z=0       ->  -z <= 0
          Math::Matrix<Real>{
            {  1 / std::sqrt(3.0),  1 / std::sqrt(3.0),  1 / std::sqrt(3.0) }, // x+y+z <= 1
            { -1,  0,  0 },                                                    // x >= 0
            {  0, -1,  0 },                                                    // y >= 0
            {  0,  0, -1 }                                                     // z >= 0
          },
          Math::Vector<Real>{{ 1 / std::sqrt(3.0), 0, 0, 0 }}
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
      case Polytope::Type::Hexahedron:
      {
        os << "Hexahedron";
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
    QF::GenericPolytopeQuadrature qf(getTransformation().getOrder(), getGeometry());
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

  Math::SpatialPoint Vertex::getCoordinates() const
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
      case Type::Hexahedron:
      {
        // Unit cube [0,1]^3: clamp each coordinate independently,
        // consistent with Segment / Quadrilateral tensor-product logic.
        out[0] = std::clamp(rc[0], Real(0), Real(1));
        out[1] = std::clamp(rc[1], Real(0), Real(1));
        out[2] = std::clamp(rc[2], Real(0), Real(1));
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
        // Inline projection to boundary of reference tetra:
        // T = { x>=0, y>=0, z>=0, x+y+z<=1 }.
        // Boundary faces (must match your HalfSpace order):
        //   face 0: x+y+z = 1
        //   face 1: x = 0
        //   face 2: y = 0
        //   face 3: z = 0

        assert(rc.size() == 3);
        out.resize(3);

        Real bestd = std::numeric_limits<Real>::infinity();

        // Candidate best point
        Real bx = Real(0), by = Real(0), bz = Real(0);

        // ---------------------------
        // Face 0: x+y+z = 1, x,y,z >= 0
        // Project (rc0,rc1,rc2) onto simplex {sum=1, nonneg} using a fixed sorting network.
        {
          const Real a = rc[0];
          const Real b = rc[1];
          const Real c = rc[2];

          // sort u0>=u1>=u2 of (a,b,c) with indices, no std::sort
          Real u0 = a, u1 = b, u2 = c;
          int  i0 = 0, i1 = 1, i2 = 2;

          if (u1 > u0) { std::swap(u0, u1); std::swap(i0, i1); }
          if (u2 > u0) { std::swap(u0, u2); std::swap(i0, i2); }
          if (u2 > u1) { std::swap(u1, u2); std::swap(i1, i2); }

          // Compute theta (standard simplex projection, equality sum=1)
          Real theta = (u0 - Real(1)); // k=1
          if (u1 - (u0 + u1 - Real(1)) / Real(2) > Real(0))
            theta = (u0 + u1 - Real(1)) / Real(2); // k=2
          if (u2 - (u0 + u1 + u2 - Real(1)) / Real(3) > Real(0))
            theta = (u0 + u1 + u2 - Real(1)) / Real(3); // k=3

          Real x = std::max(a - theta, Real(0));
          Real y = std::max(b - theta, Real(0));
          Real z = std::max(c - theta, Real(0));

          // Distance in reference coords
          const Real dx = rc[0] - x;
          const Real dy = rc[1] - y;
          const Real dz = rc[2] - z;
          const Real d  = dx*dx + dy*dy + dz*dz;

          if (d < bestd)
          {
            bestd = d;
            bx = x; by = y; bz = z;
          }
        }

        // ---------------------------
        // Face 1: x = 0, y>=0, z>=0, y+z <= 1
        {
          Real y = std::max(rc[1], Real(0));
          Real z = std::max(rc[2], Real(0));
          const Real s = y + z;
          if (s > Real(1))
          {
            const Real t = Real(0.5) * (s - Real(1));
            y = std::max(y - t, Real(0));
            z = std::max(z - t, Real(0));
          }
          const Real x = Real(0);

          const Real dx = rc[0] - x;
          const Real dy = rc[1] - y;
          const Real dz = rc[2] - z;
          const Real d  = dx*dx + dy*dy + dz*dz;

          if (d < bestd)
          {
            bestd = d;
            bx = x; by = y; bz = z;
          }
        }

        // ---------------------------
        // Face 2: y = 0, x>=0, z>=0, x+z <= 1
        {
          Real x = std::max(rc[0], Real(0));
          Real z = std::max(rc[2], Real(0));
          const Real s = x + z;
          if (s > Real(1))
          {
            const Real t = Real(0.5) * (s - Real(1));
            x = std::max(x - t, Real(0));
            z = std::max(z - t, Real(0));
          }
          const Real y = Real(0);

          const Real dx = rc[0] - x;
          const Real dy = rc[1] - y;
          const Real dz = rc[2] - z;
          const Real d  = dx*dx + dy*dy + dz*dz;

          if (d < bestd)
          {
            bestd = d;
            bx = x; by = y; bz = z;
          }
        }

        // ---------------------------
        // Face 3: z = 0, x>=0, y>=0, x+y <= 1
        {
          Real x = std::max(rc[0], Real(0));
          Real y = std::max(rc[1], Real(0));
          const Real s = x + y;
          if (s > Real(1))
          {
            const Real t = Real(0.5) * (s - Real(1));
            x = std::max(x - t, Real(0));
            y = std::max(y - t, Real(0));
          }
          const Real z = Real(0);

          const Real dx = rc[0] - x;
          const Real dy = rc[1] - y;
          const Real dz = rc[2] - z;
          const Real d  = dx*dx + dy*dy + dz*dz;

          if (d < bestd)
          {
            bestd = d;
            bx = x; by = y; bz = z;
          }
        }

        out[0] = bx;
        out[1] = by;
        out[2] = bz;
        return;
      }
      case Type::Hexahedron:
      {
        // Project to each of the 6 faces of the unit cube [0,1]^3
        const Real cx = std::clamp(rc[0], Real(0), Real(1));
        const Real cy = std::clamp(rc[1], Real(0), Real(1));
        const Real cz = std::clamp(rc[2], Real(0), Real(1));

        // z = 0
        const Real h0x = cx, h0y = cy, h0z = Real(0);
        // z = 1
        const Real h1x = cx, h1y = cy, h1z = Real(1);
        // y = 0
        const Real h2x = cx, h2y = Real(0), h2z = cz;
        // y = 1
        const Real h3x = cx, h3y = Real(1), h3z = cz;
        // x = 0
        const Real h4x = Real(0), h4y = cy, h4z = cz;
        // x = 1
        const Real h5x = Real(1), h5y = cy, h5z = cz;

        const auto sq = [](Real v) { return v * v; };

        const Real d0 = sq(rc[0]-h0x) + sq(rc[1]-h0y) + sq(rc[2]-h0z);
        const Real d1 = sq(rc[0]-h1x) + sq(rc[1]-h1y) + sq(rc[2]-h1z);
        const Real d2 = sq(rc[0]-h2x) + sq(rc[1]-h2y) + sq(rc[2]-h2z);
        const Real d3 = sq(rc[0]-h3x) + sq(rc[1]-h3y) + sq(rc[2]-h3z);
        const Real d4 = sq(rc[0]-h4x) + sq(rc[1]-h4y) + sq(rc[2]-h4z);
        const Real d5 = sq(rc[0]-h5x) + sq(rc[1]-h5y) + sq(rc[2]-h5z);

        Real bx = h0x, by = h0y, bz = h0z;
        Real bd = d0;

        if (d1 < bd) { bd = d1; bx = h1x; by = h1y; bz = h1z; }
        if (d2 < bd) { bd = d2; bx = h2x; by = h2y; bz = h2z; }
        if (d3 < bd) { bd = d3; bx = h3x; by = h3y; bz = h3z; }
        if (d4 < bd) { bd = d4; bx = h4x; by = h4y; bz = h4z; }
        if (d5 < bd) {           bx = h5x; by = h5y; bz = h5z; }

        out[0] = bx; out[1] = by; out[2] = bz;
        return;
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

        // --- Helpers (no std::sort)
        auto proj_edge01 = [](Real u, Real v, Real& ou, Real& ov)
        {
          u = std::max(u, Real(0));
          v = std::max(v, Real(0));
          const Real s = u + v;
          if (s <= Real(1)) { ou = u; ov = v; return; }
          const Real t = Real(0.5) * (s - Real(1));
          ou = std::max(u - t, Real(0));
          ov = std::max(v - t, Real(0));
        };

        // Project (a,b,c) onto {x>=0,y>=0,z>=0,x+y+z=1} without std::sort:
        // Uses a 3-element sorting network (comparisons + swaps) and the standard simplex projection formula.
        auto proj_sum1 = [](Real a, Real b, Real c, Real& x, Real& y, Real& z)
        {
          // Keep original order in v0,v1,v2; sort copies u0>=u1>=u2 with indices.
          Real u0 = a, u1 = b, u2 = c;
          int  i0 = 0, i1 = 1, i2 = 2;

          auto swap_uv = [](Real& uA, Real& uB, int& iA, int& iB)
          {
            if (uB > uA) { std::swap(uA, uB); std::swap(iA, iB); }
          };

          // Sorting network for 3 elements (descending)
          swap_uv(u0, u1, i0, i1);
          swap_uv(u0, u2, i0, i2);
          swap_uv(u1, u2, i1, i2);

          // Find rho and theta
          Real theta = Real(0);

          // k=1 candidate
          const Real t1 = (u0 - Real(1)) / Real(1);
          if (u0 - t1 > Real(0))
          {
            theta = t1;

            // k=2 candidate
            const Real t2 = (u0 + u1 - Real(1)) / Real(2);
            if (u1 - t2 > Real(0))
            {
              theta = t2;

              // k=3 candidate
              const Real t3 = (u0 + u1 + u2 - Real(1)) / Real(3);
              if (u2 - t3 > Real(0))
                theta = t3;
            }
          }

          // Apply to original components (no reordering needed)
          x = std::max(a - theta, Real(0));
          y = std::max(b - theta, Real(0));
          z = std::max(c - theta, Real(0));

          // (Optional) tiny renormalization to fight roundoff:
          // Real s = x+y+z; if (s > 0) { const Real inv = Real(1)/s; x*=inv; y*=inv; z*=inv; }
        };

        // --- Faces must match your HalfSpace order:
        // face 0: (1,2,3) -> x+y+z = 1
        // face 1: (0,3,2) -> x = 0
        // face 2: (0,1,3) -> y = 0
        // face 3: (0,2,1) -> z = 0

        if (local == 0)
        {
          Real x, y, z;
          proj_sum1(rc[0], rc[1], rc[2], x, y, z);
          out[0] = x; out[1] = y; out[2] = z;
          return;
        }
        if (local == 1)
        {
          Real y, z;
          proj_edge01(rc[1], rc[2], y, z);          // y,z >=0, y+z<=1
          out[0] = Real(0); out[1] = y; out[2] = z; // x=0
          return;
        }
        if (local == 2)
        {
          Real x, z;
          proj_edge01(rc[0], rc[2], x, z);          // x,z >=0, x+z<=1
          out[0] = x; out[1] = Real(0); out[2] = z; // y=0
          return;
        }
        // local == 3
        {
          Real x, y;
          proj_edge01(rc[0], rc[1], x, y);          // x,y >=0, x+y<=1
          out[0] = x; out[1] = y; out[2] = Real(0); // z=0
          return;
        }
      }
      case Type::Hexahedron:
      {
        // 6 faces, ordered as in Connectivity:
        // 0: bottom z=0   (0,1,2,3)
        // 1: y=0          (0,1,5,4)
        // 2: x=1          (1,2,6,5)
        // 3: y=1          (2,3,7,6)
        // 4: x=0          (3,0,4,7)
        // 5: top z=1      (4,5,6,7)
        assert(local < 6);

        const Real cx = std::clamp(rc[0], Real(0), Real(1));
        const Real cy = std::clamp(rc[1], Real(0), Real(1));
        const Real cz = std::clamp(rc[2], Real(0), Real(1));

        switch (local)
        {
          case 0: // z = 0
            out[0] = cx; out[1] = cy; out[2] = Real(0); return;
          case 1: // y = 0
            out[0] = cx; out[1] = Real(0); out[2] = cz; return;
          case 2: // x = 1
            out[0] = Real(1); out[1] = cy; out[2] = cz; return;
          case 3: // y = 1
            out[0] = cx; out[1] = Real(1); out[2] = cz; return;
          case 4: // x = 0
            out[0] = Real(0); out[1] = cy; out[2] = cz; return;
          default: // 5: z = 1
            out[0] = cx; out[1] = cy; out[2] = Real(1); return;
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
