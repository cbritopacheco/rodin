/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_GEOMETRY_MESHVIEW_H
#define RODIN_KOKKOS_GEOMETRY_MESHVIEW_H

/**
 * @file
 * @brief Borrowed, device-callable mesh topology and geometry.
 */

#include <Kokkos_Core.hpp>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/Types.h"
#include "Rodin/Types.h"

namespace Rodin::Kokkos::Geometry
{
  /**
   * @brief Device-local vector with at most three active coordinates.
   *
   * Storage is independent of host Eigen evaluation and never allocates.
   */
  class Vector
  {
    public:
      /// @brief Constructs the empty coordinate tuple of a reference point.
      KOKKOS_INLINE_FUNCTION
      Vector()
        : m_size(0),
          m_data{}
      {}

      /**
       * @brief Constructs a zero vector.
       * @param size Active dimension, between zero and three.
       */
      KOKKOS_INLINE_FUNCTION
      explicit Vector(Index size)
        : m_size(size),
          m_data{}
      {
        KOKKOS_ASSERT(size <= 3);
      }

      /// @brief Returns the active dimension. @return Coordinate count.
      KOKKOS_INLINE_FUNCTION
      Index size() const
      {
        return m_size;
      }

      /// @brief Returns the active dimension. @return Coordinate count.
      KOKKOS_INLINE_FUNCTION
      Index getDimension() const
      {
        return m_size;
      }

      /**
       * @brief Accesses a coordinate.
       * @param i Active coordinate index.
       * @return Mutable coordinate.
       */
      KOKKOS_INLINE_FUNCTION
      Real& operator()(Index i)
      {
        KOKKOS_ASSERT(i < m_size);
        return m_data[i];
      }

      /**
       * @brief Accesses a coordinate.
       * @param i Active coordinate index.
       * @return Coordinate value.
       */
      KOKKOS_INLINE_FUNCTION
      Real operator()(Index i) const
      {
        KOKKOS_ASSERT(i < m_size);
        return m_data[i];
      }

      /**
       * @brief Accesses a coordinate.
       * @param i Active coordinate index.
       * @return Mutable coordinate.
       */
      KOKKOS_INLINE_FUNCTION
      Real& operator[](Index i)
      {
        return operator()(i);
      }

      /**
       * @brief Accesses a coordinate.
       * @param i Active coordinate index.
       * @return Coordinate value.
       */
      KOKKOS_INLINE_FUNCTION
      Real operator[](Index i) const
      {
        return operator()(i);
      }

    private:
      Index m_size;
      ::Kokkos::Array<Real, 3> m_data;
  };

  /**
   * @brief Device-local matrix with at most three rows and columns.
   *
   * Active dimensions include zero for point geometry. Storage never allocates.
   */
  class Matrix
  {
    public:
      /**
       * @brief Constructs a zero matrix.
       * @param rows Active row count.
       * @param cols Active column count.
       */
      KOKKOS_INLINE_FUNCTION
      Matrix(Index rows, Index cols)
        : m_rows(rows),
          m_cols(cols),
          m_data{}
      {
        KOKKOS_ASSERT(rows <= 3 && cols <= 3);
      }

      /// @brief Returns the row count. @return Active row count.
      KOKKOS_INLINE_FUNCTION
      Index rows() const
      {
        return m_rows;
      }

      /// @brief Returns the column count. @return Active column count.
      KOKKOS_INLINE_FUNCTION
      Index cols() const
      {
        return m_cols;
      }

      /**
       * @brief Accesses a matrix entry.
       * @param row Active row index.
       * @param col Active column index.
       * @return Mutable matrix entry.
       */
      KOKKOS_INLINE_FUNCTION
      Real& operator()(Index row, Index col)
      {
        KOKKOS_ASSERT(row < m_rows && col < m_cols);
        return m_data[row * 3 + col];
      }

      /**
       * @brief Accesses a matrix entry.
       * @param row Active row index.
       * @param col Active column index.
       * @return Matrix entry value.
       */
      KOKKOS_INLINE_FUNCTION
      Real operator()(Index row, Index col) const
      {
        KOKKOS_ASSERT(row < m_rows && col < m_cols);
        return m_data[row * 3 + col];
      }

      /**
       * @brief Computes the determinant of a square matrix.
       * @return Signed determinant, with empty determinant equal to one.
       */
      KOKKOS_INLINE_FUNCTION
      Real determinant() const
      {
        KOKKOS_ASSERT(m_rows == m_cols);
        if (m_rows == 0)
          return 1;
        if (m_rows == 1)
          return operator()(0, 0);
        if (m_rows == 2)
        {
          return operator()(0, 0) * operator()(1, 1) -
                 operator()(0, 1) * operator()(1, 0);
        }
        return operator()(0, 0) *
          (operator()(1, 1) * operator()(2, 2) - operator()(1, 2) * operator()(2, 1)) -
               operator()(0, 1) *
          (operator()(1, 0) * operator()(2, 2) - operator()(1, 2) * operator()(2, 0)) +
               operator()(0, 2) *
          (operator()(1, 0) * operator()(2, 1) - operator()(1, 1) * operator()(2, 0));
      }

      /**
       * @brief Computes the inverse of a nonsingular square matrix.
       * @return Inverse matrix, including the empty inverse in dimension zero.
       */
      KOKKOS_INLINE_FUNCTION
      Matrix inverse() const
      {
        KOKKOS_ASSERT(m_rows == m_cols);
        Matrix result(m_rows, m_cols);
        const Real det = determinant();
        KOKKOS_ASSERT(det != 0);
        if (m_rows == 1)
          result(0, 0) = 1 / det;
        else if (m_rows == 2)
        {
          result(0, 0) = operator()(1, 1) / det;
          result(0, 1) = -operator()(0, 1) / det;
          result(1, 0) = -operator()(1, 0) / det;
          result(1, 1) = operator()(0, 0) / det;
        }
        else if (m_rows == 3)
        {
          for (Index i = 0; i < 3; ++i)
          {
            for (Index j = 0; j < 3; ++j)
            {
              const Index r0 = (j + 1) % 3;
              const Index r1 = (j + 2) % 3;
              const Index c0 = (i + 1) % 3;
              const Index c1 = (i + 2) % 3;
              result(i, j) = (operator()(r0, c0) * operator()(r1, c1) -
                               operator()(r0, c1) * operator()(r1, c0)) /
                det;
            }
          }
        }
        return result;
      }

    private:
      Index m_rows;
      Index m_cols;
      ::Kokkos::Array<Real, 9> m_data;
  };

  /// @brief Read-only borrowed index range in resident mesh storage.
  struct IndexView
  {
    /// @brief First entry in the range.
    const Index* data;

    /// @brief Number of entries.
    Index count;

    /// @brief Returns the range length. @return Number of entries.
    KOKKOS_INLINE_FUNCTION
    Index size() const
    {
      return count;
    }

    /**
     * @brief Reads an index.
     * @param i Position within the range.
     * @return Mesh index stored at this position.
     */
    KOKKOS_INLINE_FUNCTION
    Index operator[](Index i) const
    {
      KOKKOS_ASSERT(i < count);
      return data[i];
    }

    /**
     * @brief Reads an index.
     * @param i Position within the range.
     * @return Mesh index stored at this position.
     */
    KOKKOS_INLINE_FUNCTION
    Index operator()(Index i) const
    {
      return operator[](i);
    }
  };

  /**
   * @brief Borrowed resident arrays for one dimension-pair relation.
   * Independent CSR allocations let compute() replace a single relation.
   */
  struct IncidenceDescriptor
  {
    std::uint8_t present; ///< Whether the relation is materialized.
    const Index* offsets; ///< Row offsets, including the terminal offset.
    const Index* indices; ///< Dimension-local target indices.
    Index count; ///< Number of source rows.
  };

  /// @brief Borrowed compressed-row incidence relation.
  struct IncidenceView
  {
    /// @brief Relative row offsets, including the terminal offset.
    const Index* offsets;

    /// @brief Relation-local incident entity indices.
    const Index* indices;

    /// @brief Number of source entities.
    Index count;

    /// @brief Returns the source count. @return Number of incidence rows.
    KOKKOS_INLINE_FUNCTION
    Index size() const
    {
      return count;
    }

    /// @brief Returns the source count. @return Number of incidence rows.
    KOKKOS_INLINE_FUNCTION
    Index getCount() const
    {
      return count;
    }

    /**
     * @brief Returns incident entities for one source.
     * @param i Source entity index.
     * @return Borrowed incident entity indices.
     */
    KOKKOS_INLINE_FUNCTION
    IndexView getIncidence(Index i) const
    {
      KOKKOS_ASSERT(i < count);
      const Index begin = offsets[i];
      return {indices ? indices + begin : nullptr, offsets[i + 1] - begin};
    }

    /**
     * @brief Returns incident entities for one source.
     * @param i Source entity index.
     * @return Borrowed incident entity indices.
     */
    KOKKOS_INLINE_FUNCTION
    IndexView operator[](Index i) const
    {
      return getIncidence(i);
    }
  };

  /**
   * @brief Read-only device representation of the graded incidence complex.
   *
   * ## Architecture
   *
   * Entity arrays are flattened by dimension. Ordered vertex lists are stored
   * separately from arbitrary incidence relations, preserving reference vertex
   * orientation. Relations are independently owned CSR arrays. This view
   * neither owns storage nor computes topology. Connectivity mutation or owner
   * destruction invalidates every borrowed handle.
   */
  struct ConnectivityView
  {
    /// @brief Maximum topological dimension.
    Index dimension;

    /// @brief Dimension offsets, of length dimension plus two.
    const Index* entityOffsets;

    /// @brief Reference geometry of every flattened entity.
    const ::Rodin::Geometry::Polytope::Type* geometries;

    /// @brief Ordered vertex-list offsets per flattened entity.
    const Index* vertexOffsets;

    /// @brief Ordered mesh vertex indices.
    const Index* vertexIndices;

    /// @brief Relation descriptors, indexed by source and target dimension.
    const IncidenceDescriptor* incidenceDescriptors;

    /// @brief Returns the maximal dimension. @return Topological dimension.
    KOKKOS_INLINE_FUNCTION
    Index getDimension() const
    {
      return dimension;
    }

    /**
     * @brief Returns the entity count in a dimension.
     * @param d Entity dimension.
     * @return Number of entities.
     */
    KOKKOS_INLINE_FUNCTION
    Index getCount(Index d) const
    {
      KOKKOS_ASSERT(d <= dimension);
      return entityOffsets[d + 1] - entityOffsets[d];
    }

    /**
     * @brief Returns the reference geometry of an entity.
     * @param d Entity dimension.
     * @param i Dimension-local entity index.
     * @return Reference geometry identifier.
     */
    KOKKOS_INLINE_FUNCTION
    ::Rodin::Geometry::Polytope::Type getGeometry(Index d, Index i) const
    {
      KOKKOS_ASSERT(i < getCount(d));
      return geometries[entityOffsets[d] + i];
    }

    /**
     * @brief Returns the ordered vertices of an entity.
     * @param d Entity dimension.
     * @param i Dimension-local entity index.
     * @return Borrowed ordered vertex indices.
     */
    KOKKOS_INLINE_FUNCTION
    IndexView getPolytope(Index d, Index i) const
    {
      KOKKOS_ASSERT(i < getCount(d));
      const Index entity = entityOffsets[d] + i;
      const Index begin = vertexOffsets[entity];
      return {vertexIndices + begin, vertexOffsets[entity + 1] - begin};
    }

    /**
     * @brief Tests whether a relation was explicitly computed.
     * @param d Source dimension.
     * @param dp Target dimension.
     * @return Whether the relation is present.
     */
    KOKKOS_INLINE_FUNCTION
    bool hasIncidence(Index d, Index dp) const
    {
      KOKKOS_ASSERT(d <= dimension && dp <= dimension);
      return incidenceDescriptors[d * (dimension + 1) + dp].present != 0;
    }

    /**
     * @brief Returns an explicitly computed incidence relation.
     * @param d Source dimension.
     * @param dp Target dimension.
     * @return Borrowed compressed-row relation.
     */
    KOKKOS_INLINE_FUNCTION
    IncidenceView getIncidence(Index d, Index dp) const
    {
      KOKKOS_ASSERT(hasIncidence(d, dp));
      const auto descriptor = incidenceDescriptors[d * (dimension + 1) + dp];
      return {descriptor.offsets, descriptor.indices, descriptor.count};
    }

    /**
     * @brief Returns incident indices for a source entity.
     * @param d Source dimension.
     * @param dp Target dimension.
     * @param i Source entity index.
     * @return Borrowed incident entity indices.
     */
    KOKKOS_INLINE_FUNCTION
    IndexView getIncidence(Index d, Index dp, Index i) const
    {
      return getIncidence(d, dp).getIncidence(i);
    }
  };

  /// @brief Device-supported geometric transformation representation.
  enum class GeometryKind : std::uint8_t
  {
    ParametricP1, ///< Exact Rodin P1 geometry with independent control points.
    VertexP1, ///< Default P1 map referencing resident mesh vertices.
    Identity ///< Identity map on equal reference and physical dimensions.
  };

  struct PolytopeView;
  class PointView;

  /**
   * @brief Read-only borrowed device mesh with attached geometric maps.
   *
   * ## Architecture
   *
   * Topology, labels, coordinates and geometric control points refer to
   * allocations owned by the resident mesh. Polytopes and points are small
   * value handles into those allocations. Point evaluation uses thread-local
   * fixed-capacity numerical objects; no host objects or virtual functions are
   * accessed. Destruction or replacement of the owner invalidates this
   * view and all derived handles.
   */
  struct MeshView
  {
    /// @brief Maximum topological dimension.
    Index dimension;

    /// @brief Physical coordinate dimension.
    Index spatialDimension;

    /// @brief Dimension offsets, of length dimension plus two.
    const Index* entityOffsets;

    /// @brief Reference geometry per flattened entity.
    const ::Rodin::Geometry::Polytope::Type* geometries;

    /// @brief Ordered vertex-list offsets per flattened entity.
    const Index* vertexOffsets;

    /// @brief Ordered mesh vertex indices.
    const Index* vertexIndices;

    /// @brief Vertex-major physical coordinates.
    const Real* coordinates;

    /// @brief Attribute values per flattened entity.
    const ::Rodin::Geometry::Attribute* attributes;

    /// @brief Explicit attribute presence per flattened entity.
    const std::uint8_t* attributePresent;

    /// @brief Supported transformation kind per flattened entity.
    const GeometryKind* geometryKinds;

    /// @brief Control-point offsets per flattened entity, in point units.
    const Index* controlOffsets;

    /// @brief Point-major geometric control coordinates.
    const Real* controls;

    /// @brief Explicit relation descriptors indexed by dimension pairs.
    const IncidenceDescriptor* incidenceDescriptors;

    /// @brief Returns the mesh dimension. @return Topological dimension.
    KOKKOS_INLINE_FUNCTION
    Index getDimension() const
    {
      return dimension;
    }

    /// @brief Returns the embedding dimension. @return Coordinate dimension.
    KOKKOS_INLINE_FUNCTION
    Index getSpaceDimension() const
    {
      return spatialDimension;
    }

    /// @brief Borrows mesh topology. @return Device connectivity view.
    KOKKOS_INLINE_FUNCTION
    ConnectivityView getConnectivity() const
    {
      return {dimension, entityOffsets, geometries, vertexOffsets, vertexIndices,
        incidenceDescriptors};
    }

    /**
     * @brief Returns the number of entities in a dimension.
     * @param d Entity dimension.
     * @return Number of entities.
     */
    KOKKOS_INLINE_FUNCTION
    Index getCount(Index d) const
    {
      return getConnectivity().getCount(d);
    }

    /// @brief Returns the vertex count. @return Number of vertices.
    KOKKOS_INLINE_FUNCTION
    Index getVertexCount() const
    {
      return getCount(0);
    }

    /// @brief Returns the cell count. @return Number of maximal entities.
    KOKKOS_INLINE_FUNCTION
    Index getCellCount() const
    {
      return getCount(dimension);
    }

    /**
     * @brief Reads the physical coordinates of a mesh vertex.
     * @param i Vertex index.
     * @return Device-local coordinate vector.
     */
    KOKKOS_INLINE_FUNCTION
    Vector getVertexCoordinates(Index i) const
    {
      KOKKOS_ASSERT(i < getVertexCount());
      Vector result(spatialDimension);
      for (Index c = 0; c < spatialDimension; ++c)
        result(c) = coordinates[i * spatialDimension + c];
      return result;
    }

    /**
     * @brief Borrows a dimension-local mesh entity.
     * @param d Entity dimension.
     * @param i Entity index within its dimension.
     * @return Device polytope handle.
     */
    KOKKOS_INLINE_FUNCTION
    PolytopeView getPolytope(Index d, Index i) const;

    /**
     * @brief Borrows a maximal-dimensional entity.
     * @param i Cell index.
     * @return Device cell handle.
     */
    KOKKOS_INLINE_FUNCTION
    PolytopeView getCell(Index i) const;
  };

  /// @brief Read-only device handle to one mesh entity.
  struct PolytopeView
  {
    /// @brief Borrowed mesh storage and metadata.
    MeshView mesh;

    /// @brief Entity dimension.
    Index dimension;

    /// @brief Dimension-local entity index.
    Index index;

    /// @brief Returns the entity dimension. @return Topological dimension.
    KOKKOS_INLINE_FUNCTION
    Index getDimension() const
    {
      return dimension;
    }

    /// @brief Returns the mesh index. @return Dimension-local index.
    KOKKOS_INLINE_FUNCTION
    Index getIndex() const
    {
      return index;
    }

    /// @brief Returns the reference geometry. @return Geometry identifier.
    KOKKOS_INLINE_FUNCTION
    ::Rodin::Geometry::Polytope::Type getGeometry() const
    {
      return mesh.geometries[mesh.entityOffsets[dimension] + index];
    }

    /// @brief Returns ordered vertices. @return Borrowed vertex index range.
    KOKKOS_INLINE_FUNCTION
    IndexView getVertices() const
    {
      return mesh.getConnectivity().getPolytope(dimension, index);
    }

    /// @brief Tests attribute presence. @return Whether a label is assigned.
    KOKKOS_INLINE_FUNCTION
    bool hasAttribute() const
    {
      return mesh.attributePresent[mesh.entityOffsets[dimension] + index] != 0;
    }

    /**
     * @brief Reads an explicitly assigned attribute.
     * @return Attribute value; the entity must have an attribute.
     */
    KOKKOS_INLINE_FUNCTION
    ::Rodin::Geometry::Attribute getAttribute() const
    {
      KOKKOS_ASSERT(hasAttribute());
      return mesh.attributes[mesh.entityOffsets[dimension] + index];
    }

    /**
     * @brief Constructs a thread-local evaluation point.
     * @param reference Reference coordinates with this entity's dimension.
     * @return Device point borrowing this entity's geometry.
     */
    KOKKOS_INLINE_FUNCTION
    PointView getPoint(const Vector& reference) const;
  };

  /**
   * @brief Thread-local reference point on a device mesh entity.
   *
   * ## Architecture
   *
   * The point owns its reference coordinates and borrows immutable entity data.
   * Evaluation reproduces the source transformation, including nonaffine P1
   * tensor-product and rational pyramid geometry. Numerical results are values
   * in fixed-capacity storage, requiring neither allocation nor shared mutable
   * caches. At the pyramid apex the symmetric centreline derivative convention
   * matches Rodin's host P1 basis.
   */
  class PointView
  {
    public:
      /**
       * @brief Binds a point to an entity.
       * @param polytope Borrowed device entity.
       * @param reference Owned reference coordinate tuple.
       */
      KOKKOS_INLINE_FUNCTION
      PointView(PolytopeView polytope, const Vector& reference)
        : m_polytope(polytope),
          m_reference(reference)
      {
        KOKKOS_ASSERT(reference.size() == polytope.getDimension());
      }

      /// @brief Returns the bound entity. @return Borrowed polytope handle.
      KOKKOS_INLINE_FUNCTION
      const PolytopeView& getPolytope() const
      {
        return m_polytope;
      }

      /// @brief Returns reference coordinates. @return Owned reference tuple.
      KOKKOS_INLINE_FUNCTION
      const Vector& getReferenceCoordinates() const
      {
        return m_reference;
      }

      /**
       * @brief Evaluates the geometric map.
       * @return Physical coordinates @f$ F_K(r) @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Vector getPhysicalCoordinates() const
      {
        const auto& mesh = m_polytope.mesh;
        const Index entity = mesh.entityOffsets[m_polytope.dimension] + m_polytope.index;
        if (mesh.geometryKinds[entity] == GeometryKind::Identity)
          return m_reference;
        Vector result(mesh.spatialDimension);
        const bool vertexMap = mesh.geometryKinds[entity] == GeometryKind::VertexP1;
        const Index begin =
          vertexMap ? mesh.vertexOffsets[entity] : mesh.controlOffsets[entity];
        const Index count = vertexMap ? mesh.vertexOffsets[entity + 1] - begin
                                      : mesh.controlOffsets[entity + 1] - begin;
        for (Index a = 0; a < count; ++a)
        {
          Real value;
          Vector gradient(m_reference.size());
          evaluateBasis(a, value, gradient);
          for (Index c = 0; c < mesh.spatialDimension; ++c)
          {
            result(c) +=
              (vertexMap
                  ? mesh
                      .coordinates[mesh.vertexIndices[begin + a] * mesh.spatialDimension +
                        c]
                  : mesh.controls[(begin + a) * mesh.spatialDimension + c]) *
              value;
          }
        }
        return result;
      }

      /**
       * @brief Evaluates the geometric Jacobian.
       * @return Matrix @f$ J = \partial F_K / \partial r @f$, physical by reference.
       */
      KOKKOS_INLINE_FUNCTION
      Matrix getJacobian() const
      {
        const auto& mesh = m_polytope.mesh;
        const Index entity = mesh.entityOffsets[m_polytope.dimension] + m_polytope.index;
        Matrix result(mesh.spatialDimension, m_reference.size());
        if (mesh.geometryKinds[entity] == GeometryKind::Identity)
        {
          for (Index d = 0; d < m_reference.size(); ++d)
            result(d, d) = 1;
          return result;
        }
        const bool vertexMap = mesh.geometryKinds[entity] == GeometryKind::VertexP1;
        const Index begin =
          vertexMap ? mesh.vertexOffsets[entity] : mesh.controlOffsets[entity];
        const Index count = vertexMap ? mesh.vertexOffsets[entity + 1] - begin
                                      : mesh.controlOffsets[entity + 1] - begin;
        for (Index a = 0; a < count; ++a)
        {
          Real value;
          Vector gradient(m_reference.size());
          evaluateBasis(a, value, gradient);
          for (Index c = 0; c < mesh.spatialDimension; ++c)
          {
            const Real control = vertexMap
              ? mesh
                  .coordinates[mesh.vertexIndices[begin + a] * mesh.spatialDimension + c]
              : mesh.controls[(begin + a) * mesh.spatialDimension + c];
            for (Index d = 0; d < m_reference.size(); ++d)
              result(c, d) += control * gradient(d);
          }
        }
        return result;
      }

      /**
       * @brief Evaluates the induced reference metric.
       * @return Gram matrix @f$ J^T J @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Matrix getMetric() const
      {
        const Matrix jacobian = getJacobian();
        Matrix result(jacobian.cols(), jacobian.cols());
        for (Index i = 0; i < result.rows(); ++i)
        {
          for (Index j = 0; j < result.cols(); ++j)
          {
            for (Index c = 0; c < jacobian.rows(); ++c)
              result(i, j) += jacobian(c, i) * jacobian(c, j);
          }
        }
        return result;
      }

      /**
       * @brief Evaluates the inverse or full-column-rank pseudoinverse.
       * @return @f$ J^{-1} @f$ for square maps, otherwise @f$ (J^T J)^{-1}J^T @f$.
       *
       * Embedded entities must have full reference rank. Point geometry returns
       * the empty matrix with zero rows and physical-dimension columns.
       */
      KOKKOS_INLINE_FUNCTION
      Matrix getJacobianInverse() const
      {
        const Matrix jacobian = getJacobian();
        if (jacobian.rows() == jacobian.cols())
          return jacobian.inverse();
        const Matrix inverseMetric = getMetric().inverse();
        Matrix result(jacobian.cols(), jacobian.rows());
        for (Index i = 0; i < result.rows(); ++i)
        {
          for (Index j = 0; j < result.cols(); ++j)
          {
            for (Index k = 0; k < jacobian.cols(); ++k)
              result(i, j) += inverseMetric(i, k) * jacobian(j, k);
          }
        }
        return result;
      }

      /**
       * @brief Evaluates the signed square determinant or embedded measure.
       * @return @f$ \det J @f$ for square maps, otherwise @f$ \sqrt{\det(J^TJ)} @f$.
       */
      KOKKOS_INLINE_FUNCTION
      Real getJacobianDeterminant() const
      {
        const Matrix jacobian = getJacobian();
        if (jacobian.rows() == jacobian.cols())
          return jacobian.determinant();
        const Real det = getMetric().determinant();
        return ::Kokkos::sqrt(det >= 0 ? det : Real(0));
      }

      /**
       * @brief Evaluates the integration measure.
       * @return @f$ \sqrt{\det(J^T J)} @f$, with point measure equal to one.
       */
      KOKKOS_INLINE_FUNCTION
      Real getDistortion() const
      {
        const Matrix jacobian = getJacobian();
        const Index d = jacobian.cols();
        if (d == 0)
          return 1;
        if (jacobian.rows() == d)
          return ::Kokkos::abs(jacobian.determinant());
        if (d == 1)
        {
          Real norm = 0;
          for (Index i = 0; i < jacobian.rows(); ++i)
            norm += jacobian(i, 0) * jacobian(i, 0);
          return ::Kokkos::sqrt(norm);
        }
        if (d == 2 && jacobian.rows() == 3)
        {
          Real norm = 0;
          for (Index i = 0; i < 3; ++i)
          {
            const Index j = (i + 1) % 3;
            const Index k = (i + 2) % 3;
            const Real component =
              jacobian(j, 0) * jacobian(k, 1) - jacobian(k, 0) * jacobian(j, 1);
            norm += component * component;
          }
          return ::Kokkos::sqrt(norm);
        }
        const Real det = getMetric().determinant();
        return ::Kokkos::sqrt(det >= 0 ? det : Real(0));
      }

    private:
      /**
       * @brief Evaluates one geometry basis function and its first derivatives.
       * @param a Ordered P1 control-point index.
       * @param value Basis value at the owned reference tuple.
       * @param gradient Reference gradient at the owned reference tuple.
       */
      KOKKOS_INLINE_FUNCTION
      void evaluateBasis(Index a, Real& value, Vector& gradient) const
      {
        using Type = ::Rodin::Geometry::Polytope::Type;
        const Type type = m_polytope.getGeometry();
        const Index dimension = m_reference.size();
        if (type == Type::Point)
        {
          KOKKOS_ASSERT(a == 0);
          value = 1;
          return;
        }
        if (type == Type::Segment || type == Type::Triangle || type == Type::Tetrahedron)
        {
          KOKKOS_ASSERT(a <= dimension);
          value = a == 0 ? Real(1) : m_reference(a - 1);
          for (Index d = 0; d < dimension; ++d)
          {
            if (a == 0)
            {
              value -= m_reference(d);
              gradient(d) = -1;
            }
            else
              gradient(d) = a == d + 1 ? Real(1) : Real(0);
          }
          return;
        }
        if (type == Type::Quadrilateral || type == Type::Hexahedron)
        {
          KOKKOS_ASSERT(a < (type == Type::Quadrilateral ? 4 : 8));
          const bool upperX = a % 4 == 1 || a % 4 == 2;
          const bool upperY = a % 4 >= 2;
          const bool upperZ = a >= 4;
          const ::Kokkos::Array<bool, 3> upper{upperX, upperY, upperZ};
          value = 1;
          for (Index d = 0; d < dimension; ++d)
            value *= upper[d] ? m_reference(d) : 1 - m_reference(d);
          for (Index d = 0; d < dimension; ++d)
          {
            gradient(d) = upper[d] ? Real(1) : Real(-1);
            for (Index k = 0; k < dimension; ++k)
            {
              if (k != d)
                gradient(d) *= upper[k] ? m_reference(k) : 1 - m_reference(k);
            }
          }
          return;
        }
        if (type == Type::Wedge)
        {
          KOKKOS_ASSERT(a < 6);
          const Index triangleIndex = a % 3;
          const bool upper = a >= 3;
          const Real triangle = triangleIndex == 0 ? 1 - m_reference(0) - m_reference(1)
                                                   : m_reference(triangleIndex - 1);
          const Real segment = upper ? m_reference(2) : 1 - m_reference(2);
          value = triangle * segment;
          gradient(0) = (triangleIndex == 0      ? Real(-1)
                            : triangleIndex == 1 ? Real(1)
                                                 : Real(0)) *
            segment;
          gradient(1) = (triangleIndex == 0      ? Real(-1)
                            : triangleIndex == 2 ? Real(1)
                                                 : Real(0)) *
            segment;
          gradient(2) = triangle * (upper ? Real(1) : Real(-1));
          return;
        }
        KOKKOS_ASSERT(type == Type::Pyramid && a < 5);
        const Real x = m_reference(0);
        const Real y = m_reference(1);
        const Real z = m_reference(2);
        const Real q = 1 - z;
        constexpr Real centreline = Real(0.5);
        const Real u = q == 0 ? centreline : x / q;
        const Real v = q == 0 ? centreline : y / q;
        if (a == 4)
        {
          value = z;
          gradient(2) = 1;
          return;
        }
        const bool upperX = a == 1 || a == 2;
        const bool upperY = a >= 2;
        const Real sx = upperX ? x : q - x;
        const Real sy = upperY ? y : q - y;
        value = q == 0 ? Real(0) : sx * sy / q;
        gradient(0) = (upperX ? Real(1) : Real(-1)) * (upperY ? v : 1 - v);
        gradient(1) = (upperY ? Real(1) : Real(-1)) * (upperX ? u : 1 - u);
        gradient(2) = a == 0 ? -1 + u * v : a == 2 ? u * v : -u * v;
      }

      PolytopeView m_polytope;
      Vector m_reference;
  };

  KOKKOS_INLINE_FUNCTION
  PolytopeView MeshView::getPolytope(Index d, Index i) const
  {
    KOKKOS_ASSERT(i < getCount(d));
    return {*this, d, i};
  }

  KOKKOS_INLINE_FUNCTION
  PolytopeView MeshView::getCell(Index i) const
  {
    return getPolytope(dimension, i);
  }

  KOKKOS_INLINE_FUNCTION
  PointView PolytopeView::getPoint(const Vector& reference) const
  {
    return {*this, reference};
  }
}

#endif
