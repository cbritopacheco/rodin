/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SIMPLEXITERATOR_H
#define RODIN_GEOMETRY_SIMPLEXITERATOR_H

#include <memory>
#include <utility>

#include "ForwardDecls.h"

#include "Simplex.h"
#include "IndexGenerator.h"

namespace Rodin::Geometry
{
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

      Polytope& operator*() const noexcept;

      Polytope* operator->() const noexcept;

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
      mutable bool m_dirty;
      mutable std::unique_ptr<Polytope> m_simplex;
  };

  class ElementIterator
  {
    public:
      ElementIterator() = default;

      ElementIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      ElementIterator(const ElementIterator&) = delete;

      ElementIterator(ElementIterator&& other) = default;

      ElementIterator& operator=(ElementIterator&&) = default;

      inline
      operator bool() const
      {
        return !end();
      }

      bool end() const;

      ElementIterator& operator++();

      Element& operator*() const noexcept;

      Element* operator->() const noexcept;

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
      Element* generate() const;

      IndexGeneratorBase& getIndexGenerator()
      {
        assert(m_gen);
        return *m_gen;
      }

      std::optional<std::reference_wrapper<const MeshBase>> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable bool m_dirty;
      mutable std::unique_ptr<Element> m_simplex;
  };

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

      FaceIterator& operator++();

      Face& operator*() const noexcept;

      Face* operator->() const noexcept;

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
      mutable bool m_dirty;
      mutable std::unique_ptr<Face> m_simplex;
  };

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

      VertexIterator& operator++();

      Vertex& operator*() const noexcept;

      Vertex* operator->() const noexcept;

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
      mutable bool m_dirty;
      mutable std::unique_ptr<Vertex> m_simplex;
  };
}

#endif
