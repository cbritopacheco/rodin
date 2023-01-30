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
  SimplexIterator::SimplexIterator(
      size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen)
    :  m_dimension(dimension), m_mesh(mesh), m_gen(std::move(gen).move()),
      m_dirty(false)
  {}

  bool SimplexIterator::end() const
  {
    return getIndexGenerator().end();
  }

  SimplexIterator& SimplexIterator::operator++()
  {
    ++getIndexGenerator();
    m_dirty = true;
    return *this;
  }

  Simplex& SimplexIterator::operator*() const noexcept
  {
    if (!m_simplex || m_dirty)
      m_simplex.reset(generate());
    m_dirty = false;
    return *m_simplex;
  }

  Simplex* SimplexIterator::operator->() const noexcept
  {
    if (!m_simplex || m_dirty)
      m_simplex.reset(generate());
    m_dirty = false;
    return m_simplex.get();
  }

  Simplex* SimplexIterator::generate() const
  {
    if (end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    const auto& dimension = m_dimension;
    const auto& mesh = m_mesh.get();
    if (dimension == mesh.getDimension())
    {
      const auto& attribute = mesh.getAttribute(dimension, index);
      mfem::Array<int> vs;
      mesh.getHandle().GetElementVertices(index, vs);
      assert(vs.Size() > 0);
      return new Element(index, mesh, std::vector<Index>(vs.begin(), vs.end()), attribute);
    }
    else if (dimension == mesh.getDimension() - 1)
    {
      const auto& attribute = mesh.getAttribute(dimension, index);
      mfem::Array<int> vs;
      mesh.getHandle().GetFaceVertices(index, vs);
      assert(vs.Size() > 0);
      return new Face(index, mesh, std::vector<Index>(vs.begin(), vs.end()), attribute);
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
    const auto& dimension = getDimension();
    const auto& attribute = mesh.getAttribute(dimension, index);
    mfem::Array<int> vs;
    mesh.getHandle().GetElementVertices(index, vs);
    assert(vs.Size() > 0);
    return new Element(index, mesh, std::vector<Index>(vs.begin(), vs.end()), attribute);
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
    const auto& dimension = getDimension();
    const auto& attribute = mesh.getAttribute(dimension, index);
    mfem::Array<int> vs;
    mesh.getHandle().GetFaceVertices(index, vs);
    assert(vs.Size() > 0);
    return new Face(index, mesh, std::vector<Index>(vs.begin(), vs.end()), attribute);
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
    const auto& dimension = getDimension();
    const auto& attribute = mesh.getAttribute(dimension, index);
    Math::Vector coordinates(mesh.getSpaceDimension());
    assert(coordinates.size() > 0);
    std::copy_n(mesh.getHandle().GetVertex(index), coordinates.size(), coordinates.begin());
    return new Vertex(index, mesh, coordinates, attribute);
  }
}
