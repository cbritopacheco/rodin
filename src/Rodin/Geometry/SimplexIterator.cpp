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
    :  m_dimension(dimension), m_mesh(mesh), m_gen(std::move(gen).move()),
      m_dirty(false)
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
    if (!m_simplex || m_dirty)
      m_simplex.reset(generate());
    m_dirty = false;
    return *m_simplex;
  }

  Polytope* PolytopeIterator::operator->() const noexcept
  {
    if (!m_simplex || m_dirty)
      m_simplex.reset(generate());
    m_dirty = false;
    return m_simplex.get();
  }

  Polytope* PolytopeIterator::generate() const
  {
    if (end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    const auto& dimension = m_dimension;
    const auto& mesh = m_mesh.get();
    if (dimension == mesh.getDimension())
    {
      return new Element(index, mesh);
    }
    else if (dimension == mesh.getDimension() - 1)
    {
      return new Face(index, mesh);
    }
    else if (dimension == 0)
    {
      assert(false);
      return nullptr;
    }
    else
    {
      assert(false);
      return nullptr;
    }
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
    if (!m_simplex || m_dirty)
      m_simplex.reset(generate());
    m_dirty = false;
    return *m_simplex;
  }

  Element* ElementIterator::operator->() const noexcept
  {
    if (!m_simplex || m_dirty)
      m_simplex.reset(generate());
    m_dirty = false;
    return m_simplex.get();
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
    if (!m_simplex || m_dirty)
      m_simplex.reset(generate());
    m_dirty = false;
    return *m_simplex;
  }

  Face* FaceIterator::operator->() const noexcept
  {
    if (!m_simplex || m_dirty)
      m_simplex.reset(generate());
    m_dirty = false;
    return m_simplex.get();
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
    if (!m_simplex || m_dirty)
      m_simplex.reset(generate());
    m_dirty = false;
    return *m_simplex;
  }

  Vertex* VertexIterator::operator->() const noexcept
  {
    if (!m_simplex || m_dirty)
      m_simplex.reset(generate());
    m_dirty = false;
    return m_simplex.get();
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
    return new Vertex(index, mesh, mesh.getVertexCoordinates(index));
  }
}
