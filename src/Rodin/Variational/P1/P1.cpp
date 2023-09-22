/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "P1.h"

namespace Rodin::Variational
{
  const Geometry::GeometryIndexed<ScalarP1Element>
  P1<Scalar, Context::Serial, Geometry::Mesh<Context::Serial>>::s_elements =
  {
    { Geometry::Polytope::Type::Point, ScalarP1Element(Geometry::Polytope::Type::Point) },
    { Geometry::Polytope::Type::Segment, ScalarP1Element(Geometry::Polytope::Type::Segment) },
    { Geometry::Polytope::Type::Triangle, ScalarP1Element(Geometry::Polytope::Type::Triangle) },
    { Geometry::Polytope::Type::Quadrilateral, ScalarP1Element(Geometry::Polytope::Type::Quadrilateral) },
    { Geometry::Polytope::Type::Tetrahedron, ScalarP1Element(Geometry::Polytope::Type::Tetrahedron) }
  };

  P1<Scalar, Context::Serial, Geometry::Mesh<Context::Serial>>
  ::P1(const Geometry::Mesh<Context>& mesh)
    : m_mesh(mesh)
  {}

  P1<Math::Vector, Context::Serial, Geometry::Mesh<Context::Serial>>
  ::P1(const Geometry::Mesh<Context>& mesh, size_t vdim)
    : m_mesh(mesh), m_vdim(vdim)
  {
    const size_t vn = mesh.getVertexCount();
    assert(m_vdim < RODIN_P1_MAX_VECTOR_DIMENSION);
    m_dofs.resize(mesh.getDimension() + 1);
    for (size_t d = 1; d < mesh.getDimension() + 1; d++)
    {
      const size_t n = mesh.getConnectivity().getCount(d);
      m_dofs[d].reserve(n);
      for (size_t i = 0; i < n; i++)
      {
        const auto& polytope = mesh.getConnectivity().getPolytope(d, i);
        const size_t count = polytope.size();
        auto& dofs = m_dofs[d].emplace_back(count * vdim);
        for (size_t local = 0; local < count * vdim; local++)
        {
          const size_t q = local / vdim;
          const size_t r = local % vdim;
          assert(q < count);
          dofs.coeffRef(local) = polytope(q) + r * vn;
        }
      }
    }
  }

  namespace Internal
  {
    std::array<Geometry::GeometryIndexed<VectorP1Element>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Elements()
    {
      std::array<Geometry::GeometryIndexed<VectorP1Element>, RODIN_P1_MAX_VECTOR_DIMENSION> res;
      for (size_t i = 0; i < RODIN_P1_MAX_VECTOR_DIMENSION; i++)
      {
        res[i] =
        {
          { Geometry::Polytope::Type::Point, VectorP1Element(i, Geometry::Polytope::Type::Point) },
          { Geometry::Polytope::Type::Segment, VectorP1Element(i, Geometry::Polytope::Type::Segment) },
          { Geometry::Polytope::Type::Triangle, VectorP1Element(i, Geometry::Polytope::Type::Triangle) },
          { Geometry::Polytope::Type::Quadrilateral, VectorP1Element(i, Geometry::Polytope::Type::Quadrilateral) },
          { Geometry::Polytope::Type::Tetrahedron, VectorP1Element(i, Geometry::Polytope::Type::Tetrahedron) }
        };
      }
      return res;
    }
  }

  const std::array<Geometry::GeometryIndexed<VectorP1Element>, RODIN_P1_MAX_VECTOR_DIMENSION>
  P1<Math::Vector, Context::Serial, Geometry::Mesh<Context::Serial>>::s_elements = Internal::initVectorP1Elements();
}
