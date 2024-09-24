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
  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::initialize(size_t sdim)
  {
    m_sdim = sdim;

    m_attributes.resize(m_sdim + 1);

    // Emplace empty connectivity objects
    m_connectivity.initialize(m_sdim + 1);

    // Set indexes
    m_attributeIndex.initialize(m_sdim + 1);

    m_transformationIndex.resize(sdim + 1);
    return *this;
  }

  Mesh<Context::Local>::Builder& Mesh<Context::Local>::Builder::nodes(size_t n)
  {
    m_nodes = 0;
    m_vertices.resize(m_sdim, n);
    m_connectivity.nodes(n);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::vertex(std::initializer_list<Real> l)
  {
    assert(m_vertices.cols() > 0);
    assert(m_nodes < static_cast<size_t>(m_vertices.cols()));
    assert(l.size() == m_sdim);
    std::copy(l.begin(), l.end(), m_vertices.col(m_nodes++).begin());
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
  Mesh<Context::Local>::Builder::vertex(const Eigen::Map<const Math::Vector<Real>>& x)
  {
    assert(m_vertices.cols() > 0);
    assert(m_nodes < static_cast<size_t>(m_vertices.cols()));
    assert(x.size() >= 0);
    assert(static_cast<size_t>(x.size()) == m_sdim);
    m_vertices.col(m_nodes++) = x;
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::vertex(Math::Vector<Real>&& x)
  {
    assert(m_vertices.cols() > 0);
    assert(m_nodes < static_cast<size_t>(m_vertices.cols()));
    assert(x.size() >= 0);
    assert(static_cast<size_t>(x.size()) == m_sdim);
    m_vertices.col(m_nodes++) = std::move(x);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::vertex(const Math::Vector<Real>& x)
  {
    assert(m_vertices.cols() > 0);
    assert(m_nodes < static_cast<size_t>(m_vertices.cols()));
    assert(x.size() >= 0);
    assert(static_cast<size_t>(x.size()) == m_sdim);
    m_vertices.col(m_nodes++) = x;
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::attribute(const std::pair<size_t, Index>& p, Attribute attr)
  {
    m_attributeIndex.track(p, attr);
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
      m_attributeIndex.reserve(d, count);
      m_transformationIndex[d].write([&](auto& obj){ obj.reserve(count); });
    }
    return *this;
  }

  Mesh<Context::Local> Mesh<Context::Local>::Builder::finalize()
  {
    Mesh res;
    res.m_sdim = m_sdim;
    res.m_vertices = std::move(m_vertices);
    res.m_attributes = std::move(m_attributes);
    res.m_connectivity = std::move(m_connectivity);
    res.m_attributeIndex = std::move(m_attributeIndex);
    res.m_transformationIndex = std::move(m_transformationIndex);
    return res;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::setConnectivity(Connectivity<Context>&& connectivity)
  {
    m_connectivity = std::move(connectivity);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::setVertices(Math::Matrix<Real>&& vertices)
  {
    m_vertices = std::move(vertices);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::setAttributeIndex(AttributeIndex&& attrs)
  {
    m_attributeIndex = std::move(attrs);
    return *this;
  }

  Mesh<Context::Local>::Builder&
  Mesh<Context::Local>::Builder::setTransformationIndex(TransformationIndex&& transformations)
  {
    m_transformationIndex = std::move(transformations);
    return *this;
  }
}
