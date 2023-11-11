/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POLYTOPEITERATOR_H
#define RODIN_GEOMETRY_POLYTOPEITERATOR_H

#include <memory>
#include <utility>

#include "Rodin/Threads/Unsafe.h"

#include "ForwardDecls.h"

#include "Polytope.h"
#include "IndexGenerator.h"

namespace Rodin::Geometry
{
  /**
   * @brief Represents an iterator over a set of polytopes of a mesh.
   *
   * @warning The PolytopeIterator class is not thread safe, i.e. only one
   * thread should have access to one particular instance of the iterator at a
   * time.
   */
  class PolytopeIterator
  {
    public:
      PolytopeIterator() = default;

      PolytopeIterator(size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen);

      PolytopeIterator(const PolytopeIterator&) = delete;

      PolytopeIterator(PolytopeIterator&&) = default;

      PolytopeIterator& operator=(PolytopeIterator&&) = default;

      inline
      operator bool() const
      {
        return !end();
      }

      bool end() const;

      PolytopeIterator& operator++();

      const Polytope& operator*() const noexcept;

      const Polytope* operator->() const noexcept;

      Polytope* release();

      inline
      constexpr
      size_t getDimension() const
      {
        return m_dimension;
      }

      inline
      const MeshBase& getMesh() const
      {
        assert(m_mesh);
        return m_mesh->get();
      }

      inline
      const IndexGeneratorBase& getIndexGenerator() const
      {
        assert(m_gen);
        return *m_gen;
      }

    private:
      Polytope* generate() const;

      IndexGeneratorBase& getIndexGenerator()
      {
        assert(m_gen);
        return *m_gen;
      }

      size_t m_dimension;
      std::optional<std::reference_wrapper<const MeshBase>> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable Threads::Unsafe<bool> m_dirty;
      mutable Threads::Unsafe<std::unique_ptr<Polytope>> m_polytope;
  };

  /**
   * @brief Represents an iterator over a set of cells of a mesh.
   *
   * @warning The CellIterator class is not thread safe, i.e. only one thread
   * should have access to one particular instance of the iterator at a time.
   */
  class CellIterator
  {
    public:
      CellIterator() = default;

      CellIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      CellIterator(const CellIterator&) = delete;

      CellIterator(CellIterator&& other) = default;

      CellIterator& operator=(CellIterator&&) = default;

      inline
      operator bool() const
      {
        return !end();
      }

      bool end() const;

      Cell* release();

      CellIterator& operator++();

      const Cell& operator*() const noexcept;

      const Cell* operator->() const noexcept;

      size_t getDimension() const;

      inline
      const MeshBase& getMesh() const
      {
        assert(m_mesh);
        return m_mesh->get();
      }

      inline
      const IndexGeneratorBase& getIndexGenerator() const
      {
        assert(m_gen);
        return *m_gen;
      }

    private:
      Cell* generate() const;

      IndexGeneratorBase& getIndexGenerator()
      {
        assert(m_gen);
        return *m_gen;
      }

      std::optional<std::reference_wrapper<const MeshBase>> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable Threads::Unsafe<bool> m_dirty;
      mutable Threads::Unsafe<std::unique_ptr<Cell>> m_polytope;
  };

  /**
   * @brief Represents an iterator over a set of faces of a mesh.
   *
   * @warning The FaceIterator class is not thread safe, i.e. only one thread
   * should have access to one particular instance of the iterator at a time.
   */
  class FaceIterator
  {
    public:
      FaceIterator() = default;

      FaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      FaceIterator(const FaceIterator&) = delete;

      FaceIterator(FaceIterator&&) = default;

      FaceIterator& operator=(FaceIterator&&) = default;

      inline
      operator bool() const
      {
        return !end();
      }

      bool end() const;

      Face* release();

      FaceIterator& operator++();

      const Face& operator*() const noexcept;

      const Face* operator->() const noexcept;

      size_t getDimension() const;

      inline
      const MeshBase& getMesh() const
      {
        assert(m_mesh);
        return m_mesh->get();
      }

      inline
      const IndexGeneratorBase& getIndexGenerator() const
      {
        assert(m_gen);
        return *m_gen;
      }

    private:
      Face* generate() const;

      IndexGeneratorBase& getIndexGenerator()
      {
        assert(m_gen);
        return *m_gen;
      }

      std::optional<std::reference_wrapper<const MeshBase>> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable Threads::Unsafe<bool> m_dirty;
      mutable Threads::Unsafe<std::unique_ptr<Face>> m_polytope;
  };

  /**
   * @brief Represents an iterator over a set of vertices of a mesh.
   *
   * @warning The VertexIterator class is not thread safe, i.e. only one thread
   * should have access to one particular instance of the iterator at a time.
   */
  class VertexIterator
  {
    public:
      VertexIterator() = default;

      VertexIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      VertexIterator(const VertexIterator&) = delete;

      VertexIterator(VertexIterator&&) = default;

      VertexIterator& operator=(VertexIterator&&) = default;

      inline
      operator bool() const
      {
        return !end();
      }

      bool end() const;

      Vertex* release();

      VertexIterator& operator++();

      const Vertex& operator*() const noexcept;

      const Vertex* operator->() const noexcept;

      constexpr size_t getDimension() const;

      inline
      const MeshBase& getMesh() const
      {
        assert(m_mesh);
        return m_mesh->get();
      }

      inline
      const IndexGeneratorBase& getIndexGenerator() const
      {
        assert(m_gen);
        return *m_gen;
      }

    private:
      Vertex* generate() const;

      IndexGeneratorBase& getIndexGenerator()
      {
        assert(m_gen);
        return *m_gen;
      }

      std::optional<std::reference_wrapper<const MeshBase>> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable Threads::Unsafe<bool> m_dirty;
      mutable Threads::Unsafe<std::unique_ptr<Vertex>> m_vertex;
  };
}

#endif
