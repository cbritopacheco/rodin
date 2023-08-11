/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Mesh.h"
#include "SimplexIterator.h"

namespace Rodin::Geometry
{
  // ---- SimplexIterator ---------------------------------------------------
  PolytopeIterator::PolytopeIterator(
      size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen)
    : m_dimension(dimension), m_mesh(mesh), m_gen(std::move(gen).move()), m_dirty(false)
  {}

  bool PolytopeIterator::end() const
  {
    return getIndexGenerator().end();
  }

  PolytopeIterator& PolytopeIterator::operator++()
  {
    ++getIndexGenerator();
    m_dirty = true;
    return *this;
  }

  Polytope& PolytopeIterator::operator*() const noexcept
  {
    if (!m_polytope || m_dirty)
      m_polytope.reset(generate());
    m_dirty = false;
    return *m_polytope;
  }

  Polytope* PolytopeIterator::operator->() const noexcept
  {
    if (!m_polytope || m_dirty)
      m_polytope.reset(generate());
    m_dirty = false;
    return m_polytope.get();
  }

  Polytope* PolytopeIterator::generate() const
  {
    if (end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    const auto& dimension = m_dimension;
    const auto& mesh = getMesh();
    return new Polytope(dimension, index, mesh);
  }

  // ---- ElementIterator ---------------------------------------------------
  ElementIterator::ElementIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
    : m_mesh(mesh), m_gen(std::move(gen).move()), m_dirty(false)
  {}

  bool ElementIterator::end() const
  {
    return getIndexGenerator().end();
  }

  ElementIterator& ElementIterator::operator++()
  {
    ++getIndexGenerator();
    m_dirty = true;
    return *this;
  }

  Element& ElementIterator::operator*() const noexcept
  {
    if (!m_polytope || m_dirty)
      m_polytope.reset(generate());
    m_dirty = false;
    return *m_polytope;
  }

  Element* ElementIterator::operator->() const noexcept
  {
    if (!m_polytope || m_dirty)
      m_polytope.reset(generate());
    m_dirty = false;
    return m_polytope.get();
  }

  size_t ElementIterator::getDimension() const
  {
    return getMesh().getDimension();
  }

  Element* ElementIterator::generate() const
  {
    if (end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    const auto& mesh = getMesh();
    return new Element(index, mesh);
  }

  // ---- FaceIterator ------------------------------------------------------
  FaceIterator::FaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
    : m_mesh(mesh), m_gen(std::move(gen).move()), m_dirty(false)
  {}

  bool FaceIterator::end() const
  {
    return getIndexGenerator().end();
  }

  FaceIterator& FaceIterator::operator++()
  {
    ++getIndexGenerator();
    m_dirty = true;
    return *this;
  }

  Face& FaceIterator::operator*() const noexcept
  {
    if (!m_polytope || m_dirty)
      m_polytope.reset(generate());
    m_dirty = false;
    return *m_polytope;
  }

  Face* FaceIterator::operator->() const noexcept
  {
    if (!m_polytope || m_dirty)
      m_polytope.reset(generate());
    m_dirty = false;
    return m_polytope.get();
  }

  size_t FaceIterator::getDimension() const
  {
    return getMesh().getDimension() - 1;
  }

  Face* FaceIterator::generate() const
  {
    if (end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    const auto& mesh = getMesh();
    return new Face(index, mesh);
  }

  // ---- VertexIterator ------------------------------------------------------
  VertexIterator::VertexIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
    : m_mesh(mesh), m_gen(std::move(gen).move()), m_dirty(false)
  {}

  bool VertexIterator::end() const
  {
    return getIndexGenerator().end();
  }

  VertexIterator& VertexIterator::operator++()
  {
    ++getIndexGenerator();
    m_dirty = true;
    return *this;
  }

  Vertex& VertexIterator::operator*() const noexcept
  {
    if (!m_vertex || m_dirty)
      m_vertex.reset(generate());
    m_dirty = false;
    return *m_vertex;
  }

  Vertex* VertexIterator::operator->() const noexcept
  {
    if (!m_vertex || m_dirty)
      m_vertex.reset(generate());
    m_dirty = false;
    return m_vertex.get();
  }

  constexpr size_t VertexIterator::getDimension() const
  {
    return 0;
  }

  Vertex* VertexIterator::generate() const
  {
    if (end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    const auto& mesh = getMesh();
    return new Vertex(index, mesh);
  }
}
