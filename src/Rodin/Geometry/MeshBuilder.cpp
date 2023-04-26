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
  Mesh<Context::Serial>::Builder::Builder()
  {}

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::setReference(Mesh<Context::Serial>& mesh)
  {
    m_sdim = mesh.getSpaceDimension();

    // Track the object
    m_ref.emplace(std::ref(mesh));

    // Set indexes
    m_attrs.initialize(m_sdim);
    m_transformations.initialize(m_sdim);

    // Emplace empty connectivity objects
    m_connectivity.initialize(m_sdim);

    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::vertex(Math::Vector&& x)
  {
    assert(x.size() >= 0);
    m_vertices.push_back(std::move(x));
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::vertex(const Math::Vector& x)
  {
    assert(x.size() >= 0);
    m_vertices.push_back(x);
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
    if (d == 0) m_vertices.reserve(count);
    return *this;
  }

  void Mesh<Context::Serial>::Builder::finalize()
  {
    m_connectivity.nodes(m_vertices.size());

    assert(m_ref.has_value());
    auto& ref = m_ref->get();

    ref.m_connectivity = std::move(m_connectivity);
    ref.m_vertices = std::move(m_vertices);
    ref.m_attrs = std::move(m_attrs);
    ref.m_transformations = std::move(m_transformations);
  }
}
