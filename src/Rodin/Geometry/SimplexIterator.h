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

#include "IndexGenerator.h"
#include "Element.h"
#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  class SimplexIterator
  {
    public:
      SimplexIterator(size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen);

      SimplexIterator(const SimplexIterator&) = delete;

      SimplexIterator(SimplexIterator&&) = default;

      bool end() const;

      SimplexIterator& operator++();

      Simplex& operator*() const noexcept;

      Simplex* operator->() const noexcept;

      size_t getDimension() const
      {
        return m_dimension;
      }

      const MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

      const IndexGeneratorBase& getIndexGenerator() const
      {
        assert(m_gen);
        return *m_gen;
      }

    private:
      Simplex* generate() const;

      IndexGeneratorBase& getIndexGenerator()
      {
        assert(m_gen);
        return *m_gen;
      }

      const size_t m_dimension;
      std::reference_wrapper<const MeshBase> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable bool m_dirty;
      mutable std::unique_ptr<Simplex> m_simplex;
  };

  class ElementIterator
  {
    public:
      ElementIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      ElementIterator(const ElementIterator&) = delete;

      ElementIterator(ElementIterator&& other) = default;

      bool end() const;

      ElementIterator& operator++();

      Element& operator*() const noexcept;

      Element* operator->() const noexcept;

      size_t getDimension() const;

      const MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

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

      std::reference_wrapper<const MeshBase> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable bool m_dirty;
      mutable std::unique_ptr<Element> m_simplex;
  };

  class FaceIterator
  {
    public:
      FaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      FaceIterator(const FaceIterator&) = delete;

      FaceIterator(FaceIterator&&) = default;

      bool end() const;

      FaceIterator& operator++();

      Face& operator*() const noexcept;

      Face* operator->() const noexcept;

      size_t getDimension() const;

      const MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

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

      std::reference_wrapper<const MeshBase> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable bool m_dirty;
      mutable std::unique_ptr<Face> m_simplex;
  };

  class VertexIterator
  {
    public:
      VertexIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      VertexIterator(const VertexIterator&) = delete;

      VertexIterator(VertexIterator&&) = default;

      bool end() const;

      VertexIterator& operator++();

      Vertex& operator*() const noexcept;

      Vertex* operator->() const noexcept;

      constexpr size_t getDimension() const;

      const MeshBase& getMesh() const
      {
        return m_mesh.get();
      }

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

      std::reference_wrapper<const MeshBase> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable bool m_dirty;
      mutable std::unique_ptr<Vertex> m_simplex;
  };
}

#endif
