/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Mesh.h"
#include "PolytopeIterator.h"

namespace Rodin::Geometry
{
  // ---- PolytopeIterator ---------------------------------------------------
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

  const Polytope& PolytopeIterator::operator*() const noexcept
  {
    if (!m_polytope.read() || m_dirty.read())
    {
      Polytope* polytope = generate();
      m_polytope.write(
          [&](auto& obj)
          {
            obj.reset(polytope);
          });
    }
    m_dirty.write(
        [](auto& obj)
        {
          obj = false;
        });
    return *(m_polytope.read());
  }

  const Polytope* PolytopeIterator::operator->() const noexcept
  {
    if (!m_polytope.read() || m_dirty.read())
    {
      Polytope* polytope = generate();
      m_polytope.write(
          [&](auto& obj)
          {
            obj.reset(polytope);
          });
    }
    m_dirty.write(
        [](auto& obj)
        {
          obj = false;
        });
    return m_polytope.read().get();
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

  Polytope* PolytopeIterator::release()
  {
    if (!m_polytope.read() || m_dirty.read())
    {
      return generate();
    }
    else
    {
      Polytope* polytope = nullptr;
      m_polytope.write(
          [&](auto& obj)
          {
            polytope = obj.release();
          });
      return polytope;
    }
  }

  // ---- CellIterator -------------------------------------------------------
  CellIterator::CellIterator(const MeshBase& mesh, IndexGeneratorBase&& gen)
    : m_mesh(mesh), m_gen(std::move(gen).move()), m_dirty(false)
  {}

  bool CellIterator::end() const
  {
    return getIndexGenerator().end();
  }

  CellIterator& CellIterator::operator++()
  {
    ++getIndexGenerator();
    m_dirty = true;
    return *this;
  }

  const Cell& CellIterator::operator*() const noexcept
  {
    if (!m_polytope.read() || m_dirty.read())
    {
      Cell* polytope = generate();
      m_polytope.write(
          [&](auto& obj)
          {
            obj.reset(polytope);
          });
    }
    m_dirty.write(
        [](auto& obj)
        {
          obj = false;
        });
    return *(m_polytope.read());
  }

  const Cell* CellIterator::operator->() const noexcept
  {
    if (!m_polytope.read() || m_dirty.read())
    {
      Cell* polytope = generate();
      m_polytope.write(
          [&](auto& obj)
          {
            obj.reset(polytope);
          });
    }
    m_dirty.write(
        [](auto& obj)
        {
          obj = false;
        });
    return m_polytope.read().get();
  }

  size_t CellIterator::getDimension() const
  {
    return getMesh().getDimension();
  }

  Cell* CellIterator::generate() const
  {
    if (end()) return nullptr;
    const auto& gen = getIndexGenerator();
    const auto& index = *gen;
    const auto& mesh = getMesh();
    return new Cell(index, mesh);
  }

  Cell* CellIterator::release()
  {
    if (!m_polytope.read() || m_dirty.read())
    {
      return generate();
    }
    else
    {
      Cell* polytope = nullptr;
      m_polytope.write(
          [&](auto& obj)
          {
            polytope = obj.release();
          });
      return polytope;
    }
  }

  // ---- FaceIterator -------------------------------------------------------
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

  const Face& FaceIterator::operator*() const noexcept
  {
    if (!m_polytope.read() || m_dirty.read())
    {
      Face* polytope = generate();
      m_polytope.write(
          [&](auto& obj)
          {
            obj.reset(polytope);
          });
    }
    m_dirty.write(
        [](auto& obj)
        {
          obj = false;
        });
    return *(m_polytope.read());
  }

  const Face* FaceIterator::operator->() const noexcept
  {
    if (!m_polytope.read() || m_dirty.read())
    {
      Face* polytope = generate();
      m_polytope.write(
          [&](auto& obj)
          {
            obj.reset(polytope);
          });
    }
    m_dirty.write(
        [](auto& obj)
        {
          obj = false;
        });
    return m_polytope.read().get();
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

  Face* FaceIterator::release()
  {
    if (!m_polytope.read() || m_dirty.read())
    {
      return generate();
    }
    else
    {
      Face* polytope = nullptr;
      m_polytope.write(
          [&](auto& obj)
          {
            polytope = obj.release();
          });
      return polytope;
    }
  }

  // ---- VertexIterator -----------------------------------------------------
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

  const Vertex& VertexIterator::operator*() const noexcept
  {
    if (!m_vertex.read() || m_dirty.read())
    {
      Vertex* polytope = generate();
      m_vertex.write(
          [&](auto& obj)
          {
            obj.reset(polytope);
          });
    }
    m_dirty.write(
        [](auto& obj)
        {
          obj = false;
        });
    return *(m_vertex.read());
  }

  const Vertex* VertexIterator::operator->() const noexcept
  {
    if (!m_vertex.read() || m_dirty.read())
    {
      Vertex* polytope = generate();
      m_vertex.write(
          [&](auto& obj)
          {
            obj.reset(polytope);
          });
    }
    m_dirty.write(
        [](auto& obj)
        {
          obj = false;
        });
    return m_vertex.read().get();
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

  Vertex* VertexIterator::release()
  {
    if (!m_vertex.read() || m_dirty.read())
    {
      return generate();
    }
    else
    {
      Vertex* polytope = nullptr;
      m_vertex.write(
          [&](auto& obj)
          {
            polytope = obj.release();
          });
      return polytope;
    }
  }
}
