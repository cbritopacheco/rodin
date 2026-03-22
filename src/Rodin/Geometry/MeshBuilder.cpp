/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert/MemberFunctionException.h"

#include "Mesh.h"
#include "Rodin/Geometry/PointCloud.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Geometry
{
  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::initialize(size_t sdim)
  {
    m_sdim = sdim;
    m_connectivity.initialize(m_sdim);
    m_attributes.initialize(m_sdim);
    m_transformations.initialize(m_sdim);
    m_initialized = true;
    return *this;
  }

  Mesh<Context::Local>::Builder& Mesh<Context::Local>::Builder::nodes(size_t n)
  {
    assert(m_initialized);
    if (!m_initialized)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "initialize(size_t) must be called before setting the number of nodes.";
    }
    m_nodes = 0;
    m_connectivity.nodes(n);
    m_vertices.resize(m_sdim, n);
    m_attributes.resize(0, n);
    m_transformations.resize(0, n);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::vertex(std::initializer_list<Real> l)
  {
    assert(m_vertices.cols() > 0);
    assert(m_nodes < static_cast<size_t>(m_vertices.cols()));
    assert(l.size() == m_sdim);
    for (size_t i = 0; i < m_sdim; ++i)
      m_vertices(i, m_nodes) = *(l.begin() + i);
    m_nodes++;
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::vertex(const Real* data)
  {
    assert(m_vertices.cols() > 0);
    assert(m_nodes < static_cast<size_t>(m_vertices.cols()));
    return vertex(Eigen::Map<const Math::Vector<Real>>(data, m_sdim));
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::vertex(const Math::SpatialPoint& x)
  {
    assert(m_vertices.cols() > 0);
    assert(m_nodes < static_cast<size_t>(m_vertices.cols()));
    assert(x.size() >= 0);
    assert(static_cast<size_t>(x.size()) == m_sdim);
    for (size_t i = 0; i < m_sdim; ++i)
      m_vertices(i, m_nodes) = x(i);
    m_nodes++;
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::attribute(const std::pair<size_t, Index>& p, const Optional<Attribute>& attr)
  {
    m_attributes.set(p, attr);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::polytope(Polytope::Type t, const Array<Index>& vs)
  {
    m_connectivity.polytope(t, vs);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::polytope(Polytope::Type t, Array<Index>&& vs)
  {
    m_connectivity.polytope(t, std::move(vs));
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::reserve(size_t d, size_t count)
  {
    if (count > 0)
    {
      m_connectivity.reserve(d, count);
      m_attributes.resize(d, count);
      m_transformations.resize(d, count);
    }
    return *this;
  }

  Mesh<Context::Local> Mesh<Context::Local>::Builder::finalize()
  {
    Mesh res;
    res.m_sdim = m_sdim;
    res.m_vertices = std::move(m_vertices);
    res.m_connectivity = std::move(m_connectivity);
    res.m_attributes = std::move(m_attributes);
    res.m_transformations = std::move(m_transformations);

    return res;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::setConnectivity(Connectivity<Context>&& connectivity)
  {
    m_connectivity = std::move(connectivity);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::setVertices(const PointCloud& vertices)
  {
    m_nodes = vertices.getCount();
    m_vertices = vertices;
    m_sdim = vertices.getDimension();
    m_connectivity.nodes(m_nodes);
    m_attributes.resize(0, m_nodes);
    m_transformations.resize(0, m_nodes);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::setVertices(PointCloud&& vertices)
  {
    m_nodes = vertices.getCount();
    m_sdim = vertices.getDimension();
    m_vertices = std::move(vertices);
    m_connectivity.nodes(m_nodes);
    m_attributes.resize(0, m_nodes);
    m_transformations.resize(0, m_nodes);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::setAttributeIndex(AttributeIndex&& attrs)
  {
    m_attributes = std::move(attrs);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::setTransformationIndex(
      PolytopeTransformationIndex&& transformations)
  {
    m_transformations = std::move(transformations);
    return *this;
  }
}
