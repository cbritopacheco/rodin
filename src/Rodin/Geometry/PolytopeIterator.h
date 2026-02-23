/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POLYTOPEITERATOR_H
#define RODIN_GEOMETRY_POLYTOPEITERATOR_H

/**
 * @file
 * @brief Iterator classes for traversing mesh polytopes.
 */

#include <memory>
#include <utility>

#include "ForwardDecls.h"

#include "Polytope.h"
#include "IndexGenerator.h"

namespace Rodin::Geometry
{
  /**
   * @brief Base class template for polytope iterators.
   *
   * This CRTP (Curiously Recurring Template Pattern) base class provides
   * common iterator functionality for traversing polytopes in a mesh.
   *
   * @tparam T Type of polytope (Polytope, Cell, Face, or Vertex)
   * @tparam Derived Derived iterator class
   *
   * # Usage Example
   * @code{.cpp}
   * // Iterate over all cells in a mesh
   * for (auto it = mesh.getCell(); it; ++it)
   * {
   *   const Cell& cell = *it;
   *   // Process cell...
   * }
   * @endcode
   *
   * # Thread Safety
   * Polytope iterators are **not** thread-safe. Each thread should use its
   * own iterator instance.
   *
   * @see PolytopeIterator, CellIterator, FaceIterator, VertexIterator
   */
  template <class T, class Derived>
  class PolytopeIteratorBase
  {
    public:
      /**
       * @brief Default constructor.
       */
      PolytopeIteratorBase() = default;

      /**
       * @brief Constructs an iterator for polytopes of given dimension.
       * @param[in] dimension Dimension of polytopes to iterate
       * @param[in] mesh Mesh containing the polytopes
       * @param[in] gen Index generator for iteration order
       */
      PolytopeIteratorBase(size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen)
        : m_dimension(dimension), m_mesh(mesh), m_gen(std::move(gen).move()), m_dirty(false)
      {}

      /**
       * @brief Copy constructor (deleted).
       */
      PolytopeIteratorBase(const PolytopeIterator&) = delete;

      /**
       * @brief Move constructor.
       */
      PolytopeIteratorBase(PolytopeIteratorBase&&) = default;

      /**
       * @brief Move assignment operator.
       */
      PolytopeIteratorBase& operator=(PolytopeIteratorBase&&) = default;

      /**
       * @brief Conversion to bool (checks if iterator is valid).
       * @returns True if iterator is not at end, false otherwise
       */
      operator bool() const
      {
        return !end();
      }

      /**
       * @brief Checks if iterator has reached the end.
       * @returns True if at end, false otherwise
       */
      bool end() const
      {
        return getIndexGenerator().end();
      }

      /**
       * @brief Pre-increment operator (advances to next polytope).
       * @returns Reference to this iterator
       */
      Derived& operator++()
      {
        ++getIndexGenerator();
        m_dirty = true;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Dereferences iterator to get current polytope.
       * @returns Reference to current polytope
       */
      const T& operator*() const noexcept
      {
        if (!m_polytope || m_dirty)
          m_polytope.reset(this->generate());
        m_dirty = false;
        return *(m_polytope);
      }

      /**
       * @brief Member access operator for current polytope.
       * @returns Pointer to current polytope
       */
      const T* operator->() const noexcept
      {
        if (!m_polytope || m_dirty)
          m_polytope.reset(this->generate());
        m_dirty = false;
        return m_polytope.get();
      }

      /**
       * @brief Releases ownership of current polytope.
       * @returns Pointer to current polytope (caller takes ownership)
       *
       * @warning The returned pointer must be deleted by the caller.
       */
      T* release()
      {
        if (!m_polytope || m_dirty)
          return generate();
        else
          return m_polytope.release();
      }

      /**
       * @brief Gets the dimension of polytopes being iterated.
       * @returns Polytope dimension
       */
      constexpr
      size_t getDimension() const
      {
        return m_dimension;
      }

      /**
       * @brief Gets the mesh being iterated.
       * @returns Reference to the mesh
       */
      const MeshBase& getMesh() const
      {
        assert(m_mesh);
        return m_mesh->get();
      }

      /**
       * @brief Gets the index generator (const).
       * @returns Reference to the index generator
       */
      const IndexGeneratorBase& getIndexGenerator() const
      {
        assert(m_gen);
        return *m_gen;
      }

      /**
       * @brief Gets the index generator (mutable).
       * @returns Reference to the index generator
       */
      IndexGeneratorBase& getIndexGenerator()
      {
        assert(m_gen);
        return *m_gen;
      }

      std::pair<size_t, Index> key() const
      {
        return { this->getDimension(), *this->getIndexGenerator() };
      }

      /**
       * @brief Generates a polytope at the current iterator position.
       * @returns Pointer to newly created polytope
       *
       * Pure virtual method to be implemented by derived classes.
       */
      virtual T* generate() const = 0;

    private:
      size_t m_dimension;
      Optional<std::reference_wrapper<const MeshBase>> m_mesh;
      std::unique_ptr<IndexGeneratorBase> m_gen;
      mutable bool m_dirty;
      mutable std::unique_ptr<T> m_polytope;
  };

  /**
   * @brief Iterator over general polytopes of specified dimension.
   *
   * Provides iteration over polytopes of any dimension (0 for vertices,
   * 1 for edges, etc.) in a mesh.
   *
   * # Thread Safety
   * Not thread-safe. Each thread should use its own iterator instance.
   *
   * @see CellIterator, FaceIterator, VertexIterator
   */
  class PolytopeIterator : public PolytopeIteratorBase<Polytope, PolytopeIterator>
  {
    public:
      /**
       * @brief Parent class type.
       */
      using Parent = PolytopeIteratorBase<Polytope, PolytopeIterator>;

      /**
       * @brief Default constructor.
       */
      PolytopeIterator() = default;

      /**
       * @brief Constructs an iterator for polytopes of given dimension.
       * @param[in] dimension Dimension of polytopes to iterate
       * @param[in] mesh Mesh containing the polytopes
       * @param[in] gen Index generator for iteration order
       */
      PolytopeIterator(size_t dimension, const MeshBase& mesh, IndexGeneratorBase&& gen);

      /**
       * @brief Copy constructor (deleted).
       */
      PolytopeIterator(const PolytopeIterator&) = delete;

      /**
       * @brief Move constructor.
       */
      PolytopeIterator(PolytopeIterator&& other) = default;

      /**
       * @brief Move assignment operator.
       */
      PolytopeIterator& operator=(PolytopeIterator&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      /**
       * @brief Assigns from a CellIterator.
       * @param[in] it Cell iterator to convert from
       */
      PolytopeIterator& operator=(CellIterator it);

      /**
       * @brief Assigns from a FaceIterator.
       * @param[in] it Face iterator to convert from
       */
      PolytopeIterator& operator=(FaceIterator it);

      /**
       * @brief Assigns from a VertexIterator.
       * @param[in] it Vertex iterator to convert from
       */
      PolytopeIterator& operator=(VertexIterator it);

      /**
       * @brief Generates a polytope at current position.
       * @returns Pointer to new Polytope object
       */
      Polytope* generate() const override;
  };

  /**
   * @brief Iterator over mesh cells (highest-dimensional polytopes).
   *
   * Provides iteration over cells, which are the highest-dimensional
   * polytopes in the mesh (e.g., triangles in 2D, tetrahedra in 3D).
   *
   * # Usage Example
   * @code{.cpp}
   * for (auto it = mesh.getCell(); it; ++it)
   * {
   *   const Cell& cell = *it;
   *   std::cout << "Cell index: " << cell.getIndex() << std::endl;
   * }
   * @endcode
   *
   * # Thread Safety
   * Not thread-safe. Each thread should use its own iterator instance.
   *
   * @see PolytopeIterator, FaceIterator, VertexIterator
   */
  class CellIterator : public PolytopeIteratorBase<Cell, CellIterator>
  {
    public:
      /**
       * @brief Parent class type.
       */
      using Parent = PolytopeIteratorBase<Cell, CellIterator>;

      /**
       * @brief Default constructor.
       */
      CellIterator() = default;

      /**
       * @brief Constructs a cell iterator.
       * @param[in] mesh Mesh containing the cells
       * @param[in] gen Index generator for iteration order
       */
      CellIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      /**
       * @brief Copy constructor (deleted).
       */
      CellIterator(const CellIterator&) = delete;

      /**
       * @brief Move constructor.
       */
      CellIterator(CellIterator&& other) = default;

      /**
       * @brief Move assignment operator.
       */
      CellIterator& operator=(CellIterator&&) = default;

      /**
       * @brief Generates a cell at current position.
       * @returns Pointer to new Cell object
       */
      Cell* generate() const override;
  };

  /**
   * @brief Iterator over mesh faces (codimension-1 polytopes).
   *
   * Provides iteration over faces, which are polytopes of dimension d-1
   * where d is the mesh dimension (e.g., edges in 2D, triangles in 3D).
   *
   * # Usage Example
   * @code{.cpp}
   * for (auto it = mesh.getFace(); it; ++it)
   * {
   *   const Face& face = *it;
   *   if (face.isOnBoundary())
   *     std::cout << "Boundary face " << face.getIndex() << std::endl;
   * }
   * @endcode
   *
   * # Thread Safety
   * Not thread-safe. Each thread should use its own iterator instance.
   *
   * @see PolytopeIterator, CellIterator, VertexIterator
   */
  class FaceIterator : public PolytopeIteratorBase<Face, FaceIterator>
  {
    public:
      /**
       * @brief Parent class type.
       */
      using Parent = PolytopeIteratorBase<Face, FaceIterator>;

      /**
       * @brief Default constructor.
       */
      FaceIterator() = default;

      /**
       * @brief Constructs a face iterator.
       * @param[in] mesh Mesh containing the faces
       * @param[in] gen Index generator for iteration order
       */
      FaceIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      /**
       * @brief Copy constructor (deleted).
       */
      FaceIterator(const FaceIterator&) = delete;

      /**
       * @brief Move constructor.
       */
      FaceIterator(FaceIterator&&) = default;

      /**
       * @brief Move assignment operator.
       */
      FaceIterator& operator=(FaceIterator&&) = default;

      /**
       * @brief Generates a face at current position.
       * @returns Pointer to new Face object
       */
      Face* generate() const override;
  };

  /**
   * @brief Iterator over mesh vertices (0-dimensional polytopes).
   *
   * Provides iteration over vertices, which are the 0-dimensional
   * polytopes (nodes) in the mesh.
   *
   * # Usage Example
   * @code{.cpp}
   * for (auto it = mesh.getVertex(); it; ++it)
   * {
   *   const Vertex& vertex = *it;
   *   const auto& coords = vertex.getCoordinates();
   *   std::cout << "Vertex at (" << coords(0) << ", " << coords(1) << ")" << std::endl;
   * }
   * @endcode
   *
   * # Thread Safety
   * Not thread-safe. Each thread should use its own iterator instance.
   *
   * @see PolytopeIterator, CellIterator, FaceIterator
   */
  class VertexIterator : public PolytopeIteratorBase<Vertex, VertexIterator>
  {
    public:
      /**
       * @brief Parent class type.
       */
      using Parent = PolytopeIteratorBase<Vertex, VertexIterator>;

      /**
       * @brief Default constructor.
       */
      VertexIterator() = default;

      /**
       * @brief Constructs a vertex iterator.
       * @param[in] mesh Mesh containing the vertices
       * @param[in] gen Index generator for iteration order
       */
      VertexIterator(const MeshBase& mesh, IndexGeneratorBase&& gen);

      /**
       * @brief Copy constructor (deleted).
       */
      VertexIterator(const VertexIterator&) = delete;

      /**
       * @brief Move constructor.
       */
      VertexIterator(VertexIterator&&) = default;

      /**
       * @brief Move assignment operator.
       */
      VertexIterator& operator=(VertexIterator&&) = default;

      /**
       * @brief Generates a vertex at current position.
       * @returns Pointer to new Vertex object
       */
      Vertex* generate() const override;
  };
}

#endif
