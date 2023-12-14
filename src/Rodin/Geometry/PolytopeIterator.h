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
  template <class T, class Derived>
  class PolytopeIteratorBase
  {
    public:
      PolytopeIteratorBase() = default;

      PolytopeIteratorBase(size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen)
        : m_dimension(dimension), m_mesh(mesh), m_gen(std::move(gen).move()), m_dirty(false)
      {}

      PolytopeIteratorBase(const PolytopeIterator&) = delete;

      PolytopeIteratorBase(PolytopeIteratorBase&&) = default;

      PolytopeIteratorBase& operator=(PolytopeIteratorBase&&) = default;

      inline
      operator bool() const
      {
        return !end();
      }

      inline
      bool end() const
      {
        return getIndexGenerator().end();
      }

      Derived& operator++()
      {
        ++getIndexGenerator();
        m_dirty = true;
        return static_cast<Derived&>(*this);
      }

      const T& operator*() const noexcept
      {
        if (!m_polytope.read() || m_dirty.read())
        {
          T* polytope = generate();
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

      const T* operator->() const noexcept
      {
        if (!m_polytope.read() || m_dirty.read())
        {
          T* polytope = generate();
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

      T* release()
      {
        if (!m_polytope.read() || m_dirty.read())
        {
          return generate();
        }
        else
        {
          T* polytope = nullptr;
          m_polytope.write(
              [&](auto& obj)
              {
                polytope = obj.release();
              });
          return polytope;
        }
      }

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

      virtual T* generate() const = 0;

      IndexGeneratorBase& getIndexGenerator()
      {
        assert(m_gen);
        return *m_gen;
      }

    private:
      size_t m_dimension;
      std::optional<std::reference_wrapper<const MeshBase>> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable Threads::Unsafe<bool> m_dirty;
      mutable Threads::Unsafe<std::unique_ptr<T>> m_polytope;
  };

  /**
   * @brief Represents an iterator over a set of polytopes of a mesh.
   *
   * @warning The PolytopeIterator class is not thread safe, i.e. only one
   * thread should have access to one particular instance of the iterator at a
   * time.
   */
  class PolytopeIterator : public PolytopeIteratorBase<Polytope, PolytopeIterator>
  {
    public:
      using Parent = PolytopeIteratorBase<Polytope, PolytopeIterator>;

      PolytopeIterator() = default;

      PolytopeIterator(size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen);

      PolytopeIterator(const PolytopeIterator&) = delete;

      PolytopeIterator(PolytopeIterator&& other) = default;

      PolytopeIterator& operator=(PolytopeIterator&&) = default;

      Polytope* generate() const override;
  };

  /**
   * @brief Represents an iterator over a set of cells of a mesh.
   *
   * @warning The CellIterator class is not thread safe, i.e. only one thread
   * should have access to one particular instance of the iterator at a time.
   */
  class CellIterator : public PolytopeIteratorBase<Cell, CellIterator>
  {
    public:
      using Parent = PolytopeIteratorBase<Cell, CellIterator>;

      CellIterator() = default;

      CellIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      CellIterator(const CellIterator&) = delete;

      CellIterator(CellIterator&& other) = default;

      CellIterator& operator=(CellIterator&&) = default;

      Cell* generate() const override;
  };

  /**
   * @brief Represents an iterator over a set of faces of a mesh.
   *
   * @warning The FaceIterator class is not thread safe, i.e. only one thread
   * should have access to one particular instance of the iterator at a time.
   */
  class FaceIterator : public PolytopeIteratorBase<Face, FaceIterator>
  {
    public:
      using Parent = PolytopeIteratorBase<Face, FaceIterator>;

      FaceIterator() = default;

      FaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      FaceIterator(const FaceIterator&) = delete;

      FaceIterator(FaceIterator&&) = default;

      FaceIterator& operator=(FaceIterator&&) = default;

      Face* generate() const override;
  };

  /**
   * @brief Represents an iterator over a set of vertices of a mesh.
   *
   * @warning The VertexIterator class is not thread safe, i.e. only one thread
   * should have access to one particular instance of the iterator at a time.
   */
  class VertexIterator : public PolytopeIteratorBase<Vertex, VertexIterator>
  {
    public:
      using Parent = PolytopeIteratorBase<Vertex, VertexIterator>;

      VertexIterator() = default;

      VertexIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      VertexIterator(const VertexIterator&) = delete;

      VertexIterator(VertexIterator&&) = default;

      VertexIterator& operator=(VertexIterator&&) = default;

      Vertex* generate() const override;
  };
}

#endif
