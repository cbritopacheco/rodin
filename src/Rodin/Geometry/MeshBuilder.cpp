/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "Mesh.h"

namespace Rodin::Geometry
{
  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::initialize(size_t sdim)
  {
    m_sdim = sdim;

    // Set indexes
    m_attrs.initialize(m_sdim);
    m_transformations.initialize(m_sdim);

    // Emplace empty connectivity objects
    m_connectivity.initialize(m_sdim);

    return *this;
  }

  Mesh<Context::Serial>::Builder& Mesh<Context::Serial>::Builder::nodes(size_t n)
  {
    m_nodes = 0;
    m_vertices.resize(m_sdim, n);
    m_connectivity.nodes(n);
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::vertex(Math::Vector&& x)
  {
    assert(x.size() >= 0);
    m_vertices.col(m_nodes++) = std::move(x);
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::vertex(const Math::Vector& x)
  {
    assert(x.size() >= 0);
    m_vertices.col(m_nodes++) = x;
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::attribute(size_t d, Index idx, Attribute attr)
  {
    m_attrs.track(d, idx, attr);
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::polytope(Polytope::Geometry t, const Array<Index>& vs)
  {
    m_connectivity.polytope(t, vs);
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::polytope(Polytope::Geometry t, Array<Index>&& vs)
  {
    m_connectivity.polytope(t, std::move(vs));
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::reserve(size_t d, size_t count)
  {
    m_attrs.reserve(d, count);
    m_connectivity.reserve(d, count);
    m_transformations.reserve(d, count);
    return *this;
  }

  Mesh<Context::Serial> Mesh<Context::Serial>::Builder::finalize()
  {
    Mesh res;
    res.m_sdim = m_sdim;
    res.m_connectivity = std::move(m_connectivity);
    res.m_vertices = std::move(m_vertices);
    res.m_attrs = std::move(m_attrs);
    res.m_transformations = std::move(m_transformations);

    return res;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::setConnectivity(MeshConnectivity&& connectivity)
  {
    m_connectivity = std::move(connectivity);
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::setVertices(Math::Matrix&& vertices)
  {
    m_vertices = std::move(vertices);
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::setAttributeIndex(AttributeIndex&& attrs)
  {
    m_attrs = std::move(attrs);
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::setTransformationIndex(TransformationIndex&& transformations)
  {
    m_transformations = std::move(transformations);
    return *this;
  }
}
