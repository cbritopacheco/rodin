/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_MPI_GEOMETRY_MESH_H
#define RODIN_KOKKOS_MPI_GEOMETRY_MESH_H

/**
 * @file
 * @brief Resident Kokkos storage of a distributed mesh.
 */

#include <cstdint>
#include <functional>
#include <vector>

#include "Rodin/MPI/Geometry/Mesh.h"

#include "Rodin/Kokkos/Geometry/Mesh.h"

namespace Rodin::Geometry
{
  /**
   * @brief Shard of a distributed mesh resident in a Kokkos memory space.
   *
   * ## Architecture
   *
   * The mesh holds the shard of this process as a local resident mesh, the
   * ownership of each of its entities, and a reference to the distributed
   * host mesh it was built from. Every query is answered in the local
   * numbering of the shard, which the kernels of a local resident mesh
   * expect. The boundary, the interface and the measures are those of the
   * distributed mesh: they are read from the host mesh, which owns the
   * entities across processes.
   *
   * The device view marks the entities this process owns. Integrators skip
   * the others, so that each polytope is integrated once over the
   * processes. The host mesh must outlive this mesh, and its shard must not
   * change after this mesh is built.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class Mesh<Rodin::Context::Kokkos<Rodin::Context::MPI, Backend>> final : public MeshBase
  {
    public:
      /// @brief Parent class.
      using Parent = MeshBase;

      /// @brief Distributed execution context.
      using ContextType = Rodin::Context::Kokkos<Rodin::Context::MPI, Backend>;

      /// @brief Distributed execution context.
      using Context = ContextType;

      /// @brief Distributed host mesh.
      using HostType = Mesh<Rodin::Context::MPI>;

      /// @brief Resident mesh of the shard.
      using LocalType = Mesh<Rodin::Context::Kokkos<Rodin::Context::Local, Backend>>;

      /// @brief Execution space of the kernels.
      using ExecutionSpace = typename ContextType::ExecutionSpaceType;

      /// @brief Memory space of the resident data.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// @brief Resident connectivity of the shard.
      using ConnectivityType = typename LocalType::ConnectivityType;

      /// @brief Borrowed device view.
      using View = Rodin::Kokkos::Geometry::MeshView;

      /**
       * @param context Execution association; Kokkos must already be initialized.
       * @param host Distributed mesh, which must outlive this mesh.
       */
      Mesh(ContextType context, const HostType& host)
        : m_context(std::move(context)),
          m_host(host),
          m_shard(typename LocalType::ContextType(
                    Rodin::Context::Local(), m_context.getExecutionSpace()),
            host.getShard())
      {
        const auto& shard = host.getShard();
        const size_t dimension = shard.getDimension();
        std::vector<std::uint8_t> owned;
        for (size_t d = 0; d <= dimension; ++d)
        {
          for (Index i = 0; i < shard.getPolytopeCount(d); ++i)
            owned.push_back(shard.isOwned(d, i) ? 1 : 0);
        }
        m_owned = ::Kokkos::View<std::uint8_t*, MemorySpace>(
          ::Kokkos::view_alloc(m_context.getExecutionSpace(),
            ::Kokkos::WithoutInitializing, std::string("Rodin mesh ownership")),
          owned.size());
        ::Kokkos::deep_copy(m_context.getExecutionSpace(), m_owned,
          ::Kokkos::View<const std::uint8_t*, ::Kokkos::HostSpace,
            ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(owned.data(), owned.size()));
        m_context.getExecutionSpace().fence("Rodin mesh ownership");
      }

      /**
       * @brief Copies the mesh, sharing the host mesh and the resident ownership flags.
       * @param other Mesh to copy.
       */
      Mesh(const Mesh& other) = default;

      /**
       * @brief Moves the mesh.
       * @param other Mesh to move.
       */
      Mesh(Mesh&& other) = default;

      /// @returns Distributed execution context.
      const ContextType& getContext() const override
      {
        return m_context;
      }

      /// @returns Distributed host mesh.
      const HostType& getHost() const
      {
        return m_host.get();
      }

      /// @returns Host shard of this process.
      const Shard& getShard() const
      {
        return getHost().getShard();
      }

      /// @returns Resident mesh of the shard.
      const LocalType& getLocal() const
      {
        return m_shard;
      }

      /**
       * @brief Gets the resident copies of arrays read by kernels on this
       * mesh, such as the coefficients of grid functions.
       * @returns Transfer cache of the shard.
       */
      Rodin::Kokkos::Internal::TransferCache<MemorySpace>& getTransferCache() const
      {
        return m_shard.getTransferCache();
      }

      /**
       * @returns Borrowed resident data with the ownership of the entities;
       * kernels must finish before owner replacement.
       */
      View getView() const
      {
        View res = m_shard.getView();
        res.owned = m_owned.data();
        return res;
      }

      /**
       * @brief Gets the revision of the entities and of their attributes.
       *
       * The entities are those of the resident shard, whose revision this
       * forwards; the ownership is fixed when this mesh is built.
       *
       * @returns Revision of the resident shard.
       */
      std::uint64_t getRevision() const
      {
        return m_shard.getRevision();
      }

      /// @param d Entity dimension. @param i Local entity index. @returns Whether this process owns the entity.
      bool isOwned(size_t d, Index i) const
      {
        return getShard().isOwned(d, i);
      }

      /// @brief Gets the topological dimension. @returns Dimension of the cells of the shard.
      size_t getDimension() const override
      {
        return m_shard.getDimension();
      }

      /// @brief Gets the dimension of the ambient space. @returns Number of coordinates of each vertex.
      size_t getSpaceDimension() const override
      {
        return m_shard.getSpaceDimension();
      }

      /// @brief Gets the connectivity of the shard. @returns Resident connectivity, in the local numbering.
      ConnectivityType& getConnectivity() override
      {
        return m_shard.getConnectivity();
      }

      /// @brief Gets the connectivity of the shard. @returns Read-only resident connectivity, in the local numbering.
      const ConnectivityType& getConnectivity() const override
      {
        return m_shard.getConnectivity();
      }

      /// @param d Entity dimension. @returns Local entity count, ghosts included.
      size_t getPolytopeCount(size_t d) const override
      {
        return m_shard.getPolytopeCount(d);
      }

      /// @param g Reference family. @returns Local entity count, ghosts included.
      size_t getPolytopeCount(Polytope::Type g) const override
      {
        return m_shard.getPolytopeCount(g);
      }

      /// @param d Entity dimension. @param i Local entity index. @returns Reference family of the entity.
      Polytope::Type getGeometry(size_t d, Index i) const override
      {
        return m_shard.getGeometry(d, i);
      }

      /// @brief Whether the mesh is a submesh. @returns False.
      bool isSubMesh() const override
      {
        return false;
      }

      /// @brief Raises an exception, the mesh not being a submesh. @returns Never returns.
      SubMeshBase& asSubMesh() override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "This Kokkos mesh is not a submesh." << Alert::Raise;
        std::abort();
      }

      /// @brief Raises an exception, the mesh not being a submesh. @returns Never returns.
      const SubMeshBase& asSubMesh() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "This Kokkos mesh is not a submesh." << Alert::Raise;
        std::abort();
      }

      /**
       * @brief Raises an exception, the mesh being immutable.
       * @param[in] c Requested scaling factor, ignored.
       * @returns Never returns.
       */
      Mesh& scale([[maybe_unused]] Real c) override
      {
        immutable(__func__);
        return *this;
      }

      /// @brief Flushes the cached data of the resident shard.
      void flush() override
      {
        m_shard.flush();
      }

      /// @param vertex Local vertex index. @returns Coordinates of the vertex.
      ::Rodin::Math::SpatialPoint getVertexCoordinates(Index vertex) const override
      {
        return m_shard.getVertexCoordinates(vertex);
      }

      /**
       * @brief Raises an exception, the mesh being immutable.
       * @param[in] idx Local vertex index, ignored.
       * @param[in] s Requested coordinate, ignored.
       * @param[in] i Coordinate index, ignored.
       * @returns Never returns.
       */
      Mesh& setVertexCoordinates([[maybe_unused]] Index idx, [[maybe_unused]] Real s,
        [[maybe_unused]] size_t i) override
      {
        immutable(__func__);
        return *this;
      }

      /**
       * @brief Raises an exception, the mesh being immutable.
       * @param[in] idx Local vertex index, ignored.
       * @param[in] coords Requested coordinates, ignored.
       * @returns Never returns.
       */
      Mesh& setVertexCoordinates([[maybe_unused]] Index idx,
        [[maybe_unused]] const ::Rodin::Math::SpatialPoint& coords) override
      {
        immutable(__func__);
        return *this;
      }

      /// @param d Entity dimension. @param i Local entity index. @returns Attribute of the entity, if any.
      Optional<Attribute> getAttribute(size_t d, Index i) const override
      {
        return m_shard.getAttribute(d, i);
      }

      /**
       * @brief Raises an exception, the mesh being immutable.
       * @param[in] p Pair of (dimension, local index) of the polytope, ignored.
       * @param[in] attr Requested attribute, ignored.
       * @returns Never returns.
       */
      Mesh& setAttribute([[maybe_unused]] const std::pair<size_t, Index>& p,
        [[maybe_unused]] const Optional<Attribute>& attr) override
      {
        immutable(__func__);
        return *this;
      }

      /**
       * @brief Raises an exception, the mesh being immutable.
       * @param[in] p Pair of (dimension, local index) of the polytope, ignored.
       * @param[in] trans Requested transformation, ignored and not taken over.
       * @returns Never returns.
       */
      Mesh& setPolytopeTransformation([[maybe_unused]] const std::pair<size_t, Index> p,
        [[maybe_unused]] PolytopeTransformation* trans) override
      {
        immutable(__func__);
        return *this;
      }

      /// @param d Entity dimension. @param i Local entity index. @returns Geometric transformation of the entity.
      const PolytopeTransformation& getPolytopeTransformation(
        size_t d, Index i) const override
      {
        return m_shard.getPolytopeTransformation(d, i);
      }

      /// @param d Entity dimension. @param i Local entity index. @param formula Reference quadrature formula. @returns Quadrature of the formula on the entity.
      const PolytopeQuadrature& getQuadrature(
        size_t d, Index i, const QF::QuadratureFormula& formula) const override
      {
        return m_shard.getQuadrature(d, i, formula);
      }

      /**
       * @brief Raises an exception, the mesh being immutable.
       * @param[in] filename Path to the mesh file, ignored.
       * @param[in] fmt File format, ignored.
       * @returns Never returns.
       */
      Mesh& load([[maybe_unused]] const boost::filesystem::path& filename,
        [[maybe_unused]] IO::FileFormat fmt) override
      {
        immutable(__func__);
        return *this;
      }

      /// @param path Destination. @param format File format. Saves the distributed host mesh.
      void save(const boost::filesystem::path& path, IO::FileFormat format) const override
      {
        getHost().save(path, format);
      }

      /// @returns Owned boundary faces, in the local numbering.
      FaceIterator getBoundary() const override
      {
        return getHost().getBoundary();
      }

      /// @returns Owned interface faces, in the local numbering.
      FaceIterator getInterface() const override
      {
        return getHost().getInterface();
      }

      /// @returns Iterator over the local cells, ghosts included.
      CellIterator getCell() const override
      {
        return m_shard.getCell();
      }

      /// @returns Iterator over the local faces, ghosts included.
      FaceIterator getFace() const override
      {
        return m_shard.getFace();
      }

      /// @returns Iterator over the local vertices, ghosts included.
      VertexIterator getVertex() const override
      {
        return m_shard.getVertex();
      }

      /// @param dimension Entity dimension. @returns Iterator over the local entities of that dimension, ghosts included.
      PolytopeIterator getPolytope(size_t dimension) const override
      {
        return m_shard.getPolytope(dimension);
      }

      /// @param idx Local cell index. @returns Iterator starting at the cell.
      CellIterator getCell(Index idx) const override
      {
        return m_shard.getCell(idx);
      }

      /// @param idx Local face index. @returns Iterator starting at the face.
      FaceIterator getFace(Index idx) const override
      {
        return m_shard.getFace(idx);
      }

      /// @param idx Local vertex index. @returns Iterator starting at the vertex.
      VertexIterator getVertex(Index idx) const override
      {
        return m_shard.getVertex(idx);
      }

      /// @param dimension Entity dimension. @param idx Local entity index. @returns Iterator starting at the entity.
      PolytopeIterator getPolytope(size_t dimension, Index idx) const override
      {
        return m_shard.getPolytope(dimension, idx);
      }

      /// @param face Local face index. @returns Whether the face lies between two cells of the distributed mesh.
      bool isInterface(Index face) const override
      {
        return getHost().isInterface(face);
      }

      /// @param face Local face index. @returns Whether the face lies on the boundary of the distributed mesh.
      bool isBoundary(Index face) const override
      {
        return getHost().isBoundary(face);
      }

      /// @brief Gets the volume of the distributed mesh. @returns Measure of the three-dimensional entities over all processes.
      Real getVolume() const override
      {
        return getHost().getVolume();
      }

      /// @param attribute Attribute. @returns Measure of the three-dimensional entities with that attribute over all processes.
      Real getVolume(Attribute attribute) const override
      {
        return getHost().getVolume(attribute);
      }

      /// @param attributes Attributes. @returns Measure of the three-dimensional entities with one of the attributes over all processes.
      Real getVolume(const FlatSet<Attribute>& attributes) const override
      {
        return getHost().getVolume(attributes);
      }

      /// @brief Gets the area of the distributed mesh. @returns Measure of the two-dimensional entities over all processes.
      Real getArea() const override
      {
        return getHost().getArea();
      }

      /// @param attribute Attribute. @returns Measure of the two-dimensional entities with that attribute over all processes.
      Real getArea(Attribute attribute) const override
      {
        return getHost().getArea(attribute);
      }

      /// @param attributes Attributes. @returns Measure of the two-dimensional entities with one of the attributes over all processes.
      Real getArea(const FlatSet<Attribute>& attributes) const override
      {
        return getHost().getArea(attributes);
      }

      /// @param d Entity dimension. @returns Measure of the entities of that dimension over all processes.
      Real getMeasure(size_t d) const override
      {
        return getHost().getMeasure(d);
      }

      /// @param d Entity dimension. @param attribute Attribute. @returns Measure of the entities of that dimension with the attribute over all processes.
      Real getMeasure(size_t d, Attribute attribute) const override
      {
        return getHost().getMeasure(d, attribute);
      }

      /// @param d Entity dimension. @param attributes Attributes. @returns Measure of the entities of that dimension with one of the attributes over all processes.
      Real getMeasure(size_t d, const FlatSet<Attribute>& attributes) const override
      {
        return getHost().getMeasure(d, attributes);
      }

      /// @brief Gets the perimeter of the distributed mesh. @returns Measure of the boundary over all processes.
      Real getPerimeter() const override
      {
        return getHost().getPerimeter();
      }

      /// @param attribute Attribute. @returns Measure of the boundary faces with that attribute over all processes.
      Real getPerimeter(Attribute attribute) const override
      {
        return getHost().getPerimeter(attribute);
      }

      /// @param attributes Attributes. @returns Measure of the boundary faces with one of the attributes over all processes.
      Real getPerimeter(const FlatSet<Attribute>& attributes) const override
      {
        return getHost().getPerimeter(attributes);
      }

    private:
      /**
       * @brief Raises the exception reporting that the mesh is immutable.
       * @param[in] function Name of the rejected member function.
       */
      [[noreturn]] void immutable(const char* function) const
      {
        Alert::MemberFunctionException(*this, function)
          << "A distributed Kokkos mesh mirrors its host mesh and cannot be modified."
          << Alert::Raise;
        std::abort();
      }

      ContextType m_context;
      std::reference_wrapper<const HostType> m_host;
      LocalType m_shard;
      ::Kokkos::View<std::uint8_t*, MemorySpace> m_owned;
  };

  /**
   * @brief Deduces the distributed resident mesh of a host mesh.
   * @param[in] context Kokkos MPI execution context.
   * @param[in] host Distributed host mesh.
   */
  template <class Backend>
  Mesh(Rodin::Context::Kokkos<Rodin::Context::MPI, Backend> context,
    const Mesh<Rodin::Context::MPI>& host)
    -> Mesh<Rodin::Context::Kokkos<Rodin::Context::MPI, Backend>>;
}

#endif
