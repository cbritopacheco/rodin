/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MPI_MESH_H
#define RODIN_GEOMETRY_MPI_MESH_H

/**
 * @file
 * @brief Distributed mesh specialization for MPI contexts.
 *
 * This file defines @ref Rodin::Geometry::Mesh<Rodin::Context::MPI>, which
 * stores a rank-local shard together with distributed index mappings and
 * communication-aware geometric queries.
 */

#include <mpi.h>
#include <type_traits>

#include "Rodin/Configure.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/Shard.h"

#include "Rodin/MPI/Context/MPI.h"

#include "Rodin/Math/SpatialVector.h"

#include "Rodin/Serialization/Optional.h"

namespace Rodin::Geometry
{
  /**
   * @brief Convenience alias for the distributed mesh type.
   */
  using MPIMesh = Mesh<Context::MPI>;

  template <>
  class Mesh<Context::MPI> : public MeshBase
  {
    public:
      /**
       * @brief Builder used to construct a distributed mesh from a shard.
       */
      class Builder
      {
        public:
          /**
           * @brief Constructs a mesh builder for the given MPI context.
           *
           * The builder is used to assemble a distributed mesh from a local shard and
           * finalize it into an `MPIMesh` object.
           *
           * @param[in] context MPI context associated with the distributed mesh.
           */
          Builder(const Context::MPI& context);

          /**
           * @brief Initializes the builder with a local mesh shard.
           *
           * @param[in] shard Local shard to store in the builder.
           * @return Reference to the builder.
           */
          Builder& initialize(Shard&& shard);

          /**
           * @brief Finalizes the builder and returns the distributed mesh.
           *
           * The resulting mesh stores the shard previously provided to the builder and
           * is associated with the MPI context used to construct the builder.
           *
           * @return Distributed mesh containing the stored shard.
           */
          Mesh finalize();

        private:
          /// Rank-local shard that will be moved into the finalized mesh.
          Shard m_shard;
          /// MPI context bound to the mesh being constructed.
          Context::MPI m_context;
      };

      /**
       * @brief Constructs an empty distributed mesh associated with the given MPI context.
       *
       * @param[in] context MPI context associated with the mesh.
       */
      Mesh(const Context::MPI& context)
        : m_context(context)
      {}

      /**
       * @brief Returns the local mesh shard stored on the current MPI rank.
       *
       * The shard contains the local portion of the distributed mesh, including
       * owned and ghost polytopes, local geometry, attributes, and connectivity.
       *
       * @return Reference to the local mesh shard.
       */
      Shard& getShard();

      /**
       * @brief Returns the local mesh shard stored on the current MPI rank.
       *
       * The shard represents the local view of the distributed mesh and contains
       * all entities required for local computation, including owned and ghost
       * polytopes.
       *
       * @return Const reference to the local mesh shard.
       */
      const Shard& getShard() const;

      /**
       * @brief Returns the local shard index corresponding to a global polytope index.
       *
       * This method converts a global polytope index of the distributed mesh to the
       * corresponding local index in the shard of the current MPI rank.
       *
       * If the polytope is not present in the local shard (neither Owned nor Ghost),
       * the function returns `std::nullopt`.
       *
       * @param[in] dimension Topological dimension of the polytope.
       * @param[in] globalIdx Global index of the polytope in the distributed mesh.
       * @return Local shard index if the polytope is present, otherwise `std::nullopt`.
       */
      Optional<Index> getLocalIndex(size_t dimension, Index globalIdx) const;

      /**
       * @brief Returns the global index corresponding to a local shard polytope.
       *
       * This method converts a local polytope index of the shard into the global
       * index of the corresponding polytope in the distributed mesh.
       *
       * @param[in] dimension Topological dimension of the polytope.
       * @param[in] localIdx Local index of the polytope in the shard.
       * @return Global index of the polytope in the distributed mesh.
       */
      Index getGlobalIndex(size_t dimension, Index localIdx) const;

      /**
       * @brief Scales the coordinates of all vertices in the distributed mesh.
       *
       * This operation scales the geometry of the mesh by a factor @f$ c @f$.
       * The scaling is applied to the local shard stored on the current MPI rank.
       *
       * Each vertex coordinate @f$ x @f$ is transformed according to
       *
       * @f[
       *   x \leftarrow c\,x .
       * @f]
       *
       * Since the mesh is distributed, each rank scales the coordinates of the
       * vertices present in its local shard. Ghost vertices are updated consistently
       * through the shard representation.
       *
       * @param[in] c Scaling factor applied to the vertex coordinates.
       * @return Reference to the mesh.
       */
      Mesh& scale(Real c) override;

      /**
       * @brief Flushes pending updates on the local shard representation.
       *
       * This forwards to the underlying local shard flush logic.
       */
      void flush() override;

      /**
       * @brief Indicates whether this mesh is a submesh.
       *
       * The distributed mesh represents a full partitioned mesh and therefore
       * is not considered a submesh.
       *
       * @return Always returns `false`.
       */
      bool isSubMesh() const override;

      /**
       * @brief Returns the topological dimension of the distributed mesh.
       *
       * The dimension corresponds to the dimension of the underlying mesh
       * shard, which is identical across all MPI ranks.
       *
       * @return Topological dimension of the mesh.
       */
      size_t getDimension() const override;

      /**
       * @brief Returns the embedding space dimension of the distributed mesh.
       *
       * This is the dimension of the ambient space in which the mesh vertices
       * are embedded. The value is obtained from the local shard and is
       * identical across all MPI ranks.
       *
       * @return Embedding space dimension of the mesh.
       */
      size_t getSpaceDimension() const override;

      /**
       * @brief Returns the MPI execution context associated with this mesh.
       *
       * The context provides access to the MPI communicator and environment
       * used by the distributed mesh.
       *
       * @return Reference to the MPI context.
       */
      const Context::MPI& getContext() const override;

      /**
       * @brief Returns the global number of polytopes of a given dimension.
       *
       * This method returns the total number of polytopes of topological dimension
       * @p d in the distributed mesh.
       *
       * Each MPI rank counts the number of polytopes of dimension @p d that it
       * owns in its local shard. Ownership guarantees that every polytope of the
       * distributed mesh is counted exactly once across all ranks.
       *
       * The global count is then obtained through an MPI reduction:
       *
       * @f[
       *   N_d = \sum_{p=0}^{P-1} N_d^{(p)},
       * @f]
       *
       * where
       * - @f$ N_d^{(p)} @f$ is the number of owned polytopes of dimension @f$ d @f$
       *   on rank @f$ p @f$, and
       * - @f$ P @f$ is the number of MPI ranks.
       *
       * @param[in] d Topological dimension of the polytopes.
       * @return Global number of polytopes of dimension @p d.
       */
      size_t getPolytopeCount(size_t d) const override;

      /**
       * @brief Returns the global number of polytopes of a given geometry.
       *
       * This method returns the total number of polytopes whose geometry matches
       * @p g in the distributed mesh.
       *
       * Each MPI rank counts the number of locally owned polytopes whose geometry
       * is @p g. The global count is then obtained through an MPI reduction over
       * all ranks.
       *
       * Ownership ensures that every polytope of the distributed mesh contributes
       * exactly once to the global count.
       *
       * @param[in] g Geometry type of the polytopes.
       * @return Global number of polytopes with geometry @p g.
       */
      size_t getPolytopeCount(Polytope::Type g) const override;

      /**
       * @brief Returns an iterator over local polytopes of the given dimension.
       *
       * The iterator traverses polytopes of dimension @p dimension in the local mesh
       * shard, starting at the local index @p localIdx.
       *
       * @param[in] dimension Topological dimension of the polytopes.
       * @param[in] localIdx Local starting index in the shard.
       * @return Iterator over local shard polytopes of dimension @p dimension.
       */
      PolytopeIterator getPolytope(size_t dimension, Index localIdx) const override;

      /**
       * @brief Returns an iterator over local cells of the mesh shard.
       *
       * The iterator starts at the local cell index @p localIdx.
       *
       * @param[in] localIdx Local starting cell index in the shard.
       * @return Iterator over local shard cells.
       */
      CellIterator getCell(Index localIdx) const override;

      /**
       * @brief Returns an iterator over local faces of the mesh shard.
       *
       * The iterator starts at the local face index @p localIdx.
       *
       * @param[in] localIdx Local starting face index in the shard.
       * @return Iterator over local shard faces.
       */
      FaceIterator getFace(Index localIdx) const override;

      /**
       * @brief Returns an iterator over local vertices of the mesh shard.
       *
       * The iterator starts at the local vertex index @p localIdx.
       *
       * @param[in] localIdx Local starting vertex index in the shard.
       * @return Iterator over local shard vertices.
       */
      VertexIterator getVertex(Index localIdx) const override;

      /**
       * @brief Returns an iterator over local polytopes of the given dimension.
       *
       * The iterator traverses polytopes of dimension @p dimension in the local mesh
       * shard, starting at local index @f$ 0 @f$.
       *
       * @param[in] dimension Topological dimension of the polytopes.
       * @return Iterator over local shard polytopes of dimension @p dimension.
       */
      PolytopeIterator getPolytope(size_t dimension) const override;

      /**
       * @brief Returns an iterator over local cells of the mesh shard.
       *
       * The iterator starts at local cell index @f$ 0 @f$.
       *
       * @return Iterator over local shard cells.
       */
      CellIterator getCell() const override;

      /**
       * @brief Returns an iterator over local faces of the mesh shard.
       *
       * The iterator starts at local face index @f$ 0 @f$.
       *
       * @return Iterator over local shard faces.
       */
      FaceIterator getFace() const override;

      /**
       * @brief Returns an iterator over local vertices of the mesh shard.
       *
       * The iterator starts at local vertex index @f$ 0 @f$.
       *
       * @return Iterator over local shard vertices.
       */
      VertexIterator getVertex() const override;

      /**
       * @brief Returns an iterator over the boundary faces of the distributed mesh.
       *
       * Faces are indexed using the local indexing of the mesh shard.
       *
       * A face @f$ f @f$ is considered a boundary face of the distributed mesh if
       *
       * @f[
       *   f \text{ is owned by the current rank and } f \text{ is marked as a boundary face in the shard topology.}
       * @f]
       *
       * @note Since shard boundary classification is computed on the full
       * shard topology, including ghost entities, this corresponds to the true
       * boundary of the distributed mesh rather than the boundary of the cut
       * shard alone. Restricting to owned faces ensures that each global
       * boundary face of the distributed mesh appears exactly once across all
       * MPI ranks.
       *
       * @return Iterator over locally owned boundary faces.
       */
      FaceIterator getBoundary() const override;

      /**
       * @brief Returns an iterator over the interface faces of the distributed mesh.
       *
       * Faces are indexed using the local indexing of the mesh shard.
       *
       * A face @f$ f @f$ is considered an interface face if
       *
       * @f[
       *   f \text{ is owned by the current rank and } f \text{ is marked as an interface face in the shard.}
       * @f]
       *
       * Interface faces correspond to faces shared between mesh partitions. The
       * ownership restriction ensures that each global interface face is represented
       * exactly once across all MPI ranks.
       *
       * @return Iterator over locally owned interface faces.
       */
      FaceIterator getInterface() const override;

      /**
       * @brief Tests whether a local face lies on a partition interface.
       *
       * The index @p faceIdx refers to the local index of the face in the mesh shard.
       *
       * A face @f$ f @f$ is considered an interface face of the distributed mesh if
       *
       * @f[
       *   f \text{ is owned by the current rank and } f \text{ is marked as an interface face in the shard.}
       * @f]
       *
       * @param faceIdx Local index of the face in the mesh shard.
       * @return `true` if the face is an owned interface face, `false` otherwise.
       */
      bool isInterface(Index faceIdx) const override;

      /**
       * @brief Tests whether a local face lies on the boundary of the distributed mesh.
       *
       * The index @p faceIdx refers to the local index of the face in the mesh shard.
       *
       * A face @f$ f @f$ is considered a boundary face of the distributed mesh if
       *
       * @f[
       *   f \text{ is owned by the current rank and } f \text{ is marked as a boundary face in the shard topology.}
       * @f]
       *
       * @note Since shard boundary classification is computed on the full
       * shard topology, including ghost entities, this corresponds to the true
       * boundary of the distributed mesh.
       *
       * @param faceIdx Local index of the face in the mesh shard.
       * @return `true` if the face is an owned boundary face, `false` otherwise.
       */
      bool isBoundary(Index faceIdx) const override;

      /**
       * @brief Returns the total volume of the distributed mesh.
       *
       * Each MPI rank accumulates the volume of its locally owned cells, and the
       * global volume is obtained by summing these contributions across all ranks.
       *
       * @return Total volume of the distributed mesh.
       */
      Real getVolume() const override;

      /**
       * @brief Returns the total volume of the distributed mesh restricted to a given attribute.
       *
       * Each MPI rank accumulates the volume of its locally owned cells whose
       * attribute is equal to @p attr, and the global volume is obtained by an MPI
       * reduction.
       *
       * @param[in] attr Attribute used to filter the cells.
       * @return Total volume of the matching cells.
       */
      Real getVolume(Attribute attr) const override;

      /**
       * @brief Returns the total volume of the distributed mesh restricted to a set of attributes.
       *
       * Each MPI rank accumulates the volume of its locally owned cells whose
       * attribute belongs to @p attrs, and the global volume is obtained by an MPI
       * reduction.
       *
       * @param[in] attrs Set of attributes used to filter the cells.
       * @return Total volume of the matching cells.
       */
      Real getVolume(const FlatSet<Attribute>& attrs) const override;

      /**
       * @brief Returns the total perimeter of the distributed mesh.
       *
       * Each MPI rank accumulates the measure of its locally owned boundary faces,
       * and the global perimeter is obtained by summing these contributions across
       * all ranks.
       *
       * @return Total perimeter of the distributed mesh.
       */
      Real getPerimeter() const override;

      /**
       * @brief Returns the total perimeter of the distributed mesh restricted to a given attribute.
       *
       * Each MPI rank accumulates the measure of its locally owned boundary faces
       * whose attribute is equal to @p attr, and the global perimeter is obtained by
       * an MPI reduction.
       *
       * @param[in] attr Attribute used to filter the boundary faces.
       * @return Total perimeter of the matching boundary faces.
       */
      Real getPerimeter(Attribute attr) const override;

      /**
       * @brief Returns the total perimeter of the distributed mesh restricted to a set of attributes.
       *
       * Each MPI rank accumulates the measure of its locally owned boundary faces
       * whose attribute belongs to @p attrs, and the global perimeter is obtained by
       * an MPI reduction.
       *
       * @param[in] attrs Set of attributes used to filter the boundary faces.
       * @return Total perimeter of the matching boundary faces.
       */
      Real getPerimeter(const FlatSet<Attribute>& attrs) const override;

      /**
       * @brief Returns the total area of the distributed mesh.
       *
       * Each MPI rank accumulates the area of its locally owned two-dimensional
       * polytopes, and the global area is obtained by summing these contributions
       * across all ranks.
       *
       * @return Total area of the distributed mesh.
       */
      Real getArea() const override;

      /**
       * @brief Returns the total area of the distributed mesh restricted to a given attribute.
       *
       * Each MPI rank accumulates the area of its locally owned two-dimensional
       * polytopes whose attribute is equal to @p attr, and the global area is
       * obtained by an MPI reduction.
       *
       * @param[in] attr Attribute used to filter the polytopes.
       * @return Total area of the matching polytopes.
       */
      Real getArea(Attribute attr) const override;

      /**
       * @brief Returns the total area of the distributed mesh restricted to a set of attributes.
       *
       * Each MPI rank accumulates the area of its locally owned two-dimensional
       * polytopes whose attribute belongs to @p attrs, and the global area is
       * obtained by an MPI reduction.
       *
       * @param[in] attrs Set of attributes used to filter the polytopes.
       * @return Total area of the matching polytopes.
       */
      Real getArea(const FlatSet<Attribute>& attrs) const override;

      /**
       * @brief Returns the total measure of the distributed mesh in dimension @p d.
       *
       * Each MPI rank accumulates the measure of its locally owned polytopes of
       * dimension @p d, and the global measure is obtained by summing these
       * contributions across all ranks.
       *
       * @param[in] d Topological dimension of the polytopes.
       * @return Total measure of the polytopes of dimension @p d.
       */
      Real getMeasure(size_t d) const override;

      /**
       * @brief Returns the total measure of the distributed mesh in dimension @p d
       * restricted to a given attribute.
       *
       * Each MPI rank accumulates the measure of its locally owned polytopes of
       * dimension @p d whose attribute is equal to @p attr, and the global measure
       * is obtained by an MPI reduction.
       *
       * @param[in] d Topological dimension of the polytopes.
       * @param[in] attr Attribute used to filter the polytopes.
       * @return Total measure of the matching polytopes.
       */
      Real getMeasure(size_t d, Attribute attr) const override;

      /**
       * @brief Returns the total measure of the distributed mesh in dimension @p d
       * restricted to a set of attributes.
       *
       * Each MPI rank accumulates the measure of its locally owned polytopes of
       * dimension @p d whose attribute belongs to @p attrs, and the global measure
       * is obtained by an MPI reduction.
       *
       * @param[in] d Topological dimension of the polytopes.
       * @param[in] attrs Set of attributes used to filter the polytopes.
       * @return Total measure of the matching polytopes.
       */
      Real getMeasure(size_t d, const FlatSet<Attribute>& attrs) const override;

      /**
       * @brief Returns the geometry type of a local polytope.
       *
       * @param[in] dimension Topological dimension of the polytope.
       * @param[in] localIdx Local index of the polytope in the shard.
       * @return Geometry type of the polytope.
       */
      Polytope::Type getGeometry(size_t dimension, Index localIdx) const override;

      /**
       * @brief Returns the attribute of a local polytope.
       *
       * @param[in] dimension Topological dimension of the polytope.
       * @param[in] localIdx Local index of the polytope in the shard.
       * @return Attribute of the polytope, or `std::nullopt` if no attribute is assigned.
       */
      Optional<Attribute> getAttribute(size_t dimension, Index localIdx) const override;

      /**
       * @brief Sets the attribute of a local polytope.
       *
       * The key @f$ ( d, i ) @f$ is interpreted using the local indexing of the mesh
       * shard.
       *
       * @param[in] p Pair @f$ ( d, i ) @f$ identifying the local polytope.
       * @param[in] attr Attribute to assign to the polytope.
       * @return Reference to the mesh.
       */
      Mesh& setAttribute(const std::pair<size_t, Index>& p, const Optional<Attribute>& attr) override;

      /**
       * @brief Sets one coordinate of a local vertex.
       *
       * @param[in] localIdx Local index of the vertex in the shard.
       * @param[in] s Value assigned to the selected coordinate.
       * @param[in] i Coordinate index in the embedding space.
       * @return Reference to the mesh.
       */
      Mesh& setVertexCoordinates(Index localIdx, Real s, size_t i) override;

      /**
       * @brief Sets the coordinates of a local vertex.
       *
       * @param[in] localIdx Local index of the vertex in the shard.
       * @param[in] coords New coordinates of the vertex.
       * @return Reference to the mesh.
       */
      Mesh& setVertexCoordinates(Index localIdx, const Math::SpatialPoint& coords) override;

      /**
       * @brief Returns the transformation associated with a local polytope.
       *
       * @param[in] dimension Topological dimension of the polytope.
       * @param[in] localIdx Local index of the polytope in the shard.
       * @return Transformation associated with the polytope.
       */
      virtual const PolytopeTransformation& getPolytopeTransformation(size_t dimension, Index localIdx) const override;

      /**
       * @brief Sets the transformation associated with a local polytope.
       *
       * The key @f$ ( d, i ) @f$ is interpreted using the local indexing of the mesh
       * shard.
       *
       * @param[in] p Pair @f$ ( d, i ) @f$ identifying the local polytope.
       * @param[in] trans Transformation assigned to the polytope.
       * @return Reference to the mesh.
       */
      virtual Mesh& setPolytopeTransformation(const std::pair<size_t, Index> p, PolytopeTransformation* trans) override;

      /**
       * @brief Loads the local shard of the distributed mesh from a file.
       *
       * Each MPI rank loads its own shard from the given file.
       *
       * @param[in] filename Path to the input file.
       * @param[in] fmt File format used to read the mesh.
       * @return Reference to the mesh.
       */
      MPIMesh& load(
        const boost::filesystem::path& filename, IO::FileFormat fmt) override;

      template <class Filename, typename = std::enable_if_t<std::is_invocable_v<Filename, int>>>
      /**
       * @brief Rank-dependent overload for loading distributed shard files.
       *
       * The callable is invoked with the current rank and must return the
       * path from which that rank should load its local shard.
       *
       * @tparam Filename Callable type invocable as `filename(int rank)`.
       * @param[in] filename Rank-dependent path generator.
       * @param[in] fmt File format used to read the shard file.
       * @return Reference to the mesh.
       */
      MPIMesh& load(const Filename& filename, IO::FileFormat fmt)
      {
        const auto& comm = m_context.getCommunicator();
        const boost::filesystem::path& path(filename(comm.rank()));
        return this->load(path, fmt);
      }

      /**
       * @brief Saves the local shard of the distributed mesh to a file.
       *
       * Each MPI rank writes its own shard to @p filename using @p fmt.
       *
       * @param[in] filename Path to the output file.
       * @param[in] fmt File format used to write the mesh.
       */
      void save(const boost::filesystem::path& filename, IO::FileFormat fmt) const override;

      template <class Filename, typename = std::enable_if_t<std::is_invocable_v<Filename, int>>>
      /**
       * @brief Rank-dependent overload for saving distributed shard files.
       *
       * The callable is invoked with the current rank and must return the
       * path to which that rank should save its local shard.
       *
       * @tparam Filename Callable type invocable as `filename(int rank)`.
       * @param[in] filename Rank-dependent path generator.
       * @param[in] fmt File format used to write the shard file.
       */
      void save(const Filename& filename, IO::FileFormat fmt)
      {
        const auto& comm = m_context.getCommunicator();
        const boost::filesystem::path& path(filename(comm.rank()));
        this->save(path, fmt);
      }

      /**
       * @brief Returns the coordinates of a local vertex.
       *
       * @param[in] localIdx Local index of the vertex in the shard.
       * @return Coordinates of the vertex in the embedding space.
       */
      Math::SpatialPoint getVertexCoordinates(Index localIdx) const override;

      /**
       * @brief Submesh downcast is unsupported for distributed meshes.
       * @throws std::runtime_error always.
       */
      SubMeshBase& asSubMesh() override
      {
        throw std::runtime_error("asSubMesh() not implemented");
      }

      /**
       * @brief Const submesh downcast is unsupported for distributed meshes.
       * @throws std::runtime_error always.
       */
      const SubMeshBase& asSubMesh() const override
      {
        throw std::runtime_error("asSubMesh() not implemented");
      }

      /**
       * @brief Returns mutable connectivity data of the local shard.
       */
      Connectivity<Context::Local>& getConnectivity() override;

      /**
       * @brief Returns const connectivity data of the local shard.
       */
      const Connectivity<Context::Local>& getConnectivity() const override;

      /**
       * @brief Reconciles entities of dimension @p d across MPI ranks.
       *
       * After a local connectivity computation, entities of dimension @p d
       * discovered independently on neighboring shards may not yet have
       * consistent distributed indices or ownership. This method resolves that
       * inconsistency.
       *
       * The reconciliation step:
       * - identifies entities shared across ranks,
       * - assigns a unique distributed index,
       * - determines a unique owner rank,
       * - updates shard metadata (ownership flags, owner and halo maps, and
       *   shard-to-distributed index mapping).
       *
       * This method is typically called after local connectivity discovery:
       *
       * @code{.cpp}
       * mesh.getConnectivity().compute(d, dp);
       * mesh.reconcile(d);
       * @endcode
       *
       * Shared entities are identified through their vertex sets, therefore
       * vertices must already have consistent distributed indices.
       *
       * Communication is restricted to neighboring ranks that may share
       * entities of dimension @p d. Interior entities remain purely local.
       *
       * Only shard metadata is modified; the local mesh topology is unchanged.
       *
       * @param[in] d Topological dimension of the entities to reconcile.
       * @return Reference to the mesh.
       */
      Mesh& reconcile(size_t d);

    private:
      /// MPI execution context associated with this distributed mesh.
      Context::MPI m_context;
      /// Rank-local shard containing geometry, topology, and ownership metadata.
      Shard m_shard;
  };
}

#endif
