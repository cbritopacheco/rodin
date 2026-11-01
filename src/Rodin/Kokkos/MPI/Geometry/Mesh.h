/*
 *          Copyright Carlos BRITO PACHECO 2026.
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
      using DeviceView = Rodin::Kokkos::Geometry::MeshView;

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

      Mesh(const Mesh& other) = default;

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
       * @returns Borrowed resident data with the ownership of the entities;
       * kernels must finish before owner replacement.
       */
      DeviceView getDeviceView() const
      {
        DeviceView res = m_shard.getDeviceView();
        res.owned = m_owned.data();
        return res;
      }

      /// @param d Entity dimension. @param i Local entity index. @returns Whether this process owns the entity.
      bool isOwned(size_t d, Index i) const
      {
        return getShard().isOwned(d, i);
      }

      size_t getDimension() const override
      {
        return m_shard.getDimension();
      }

      size_t getSpaceDimension() const override
      {
        return m_shard.getSpaceDimension();
      }

      ConnectivityType& getConnectivity() override
      {
        return m_shard.getConnectivity();
      }

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

      Polytope::Type getGeometry(size_t d, Index i) const override
      {
        return m_shard.getGeometry(d, i);
      }

      bool isSubMesh() const override
      {
        return false;
      }

      SubMeshBase& asSubMesh() override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "This Kokkos mesh is not a submesh." << Alert::Raise;
        std::abort();
      }

      const SubMeshBase& asSubMesh() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "This Kokkos mesh is not a submesh." << Alert::Raise;
        std::abort();
      }

      Mesh& scale(Real) override
      {
        immutable(__func__);
        return *this;
      }

      void flush() override
      {
        m_shard.flush();
      }

      Math::SpatialPoint getVertexCoordinates(Index vertex) const override
      {
        return m_shard.getVertexCoordinates(vertex);
      }

      Mesh& setVertexCoordinates(Index, Real, size_t) override
      {
        immutable(__func__);
        return *this;
      }

      Mesh& setVertexCoordinates(Index, const Math::SpatialPoint&) override
      {
        immutable(__func__);
        return *this;
      }

      Optional<Attribute> getAttribute(size_t d, Index i) const override
      {
        return m_shard.getAttribute(d, i);
      }

      Mesh& setAttribute(
        const std::pair<size_t, Index>&, const Optional<Attribute>&) override
      {
        immutable(__func__);
        return *this;
      }

      Mesh& setPolytopeTransformation(
        const std::pair<size_t, Index>, PolytopeTransformation*) override
      {
        immutable(__func__);
        return *this;
      }

      const PolytopeTransformation& getPolytopeTransformation(
        size_t d, Index i) const override
      {
        return m_shard.getPolytopeTransformation(d, i);
      }

      const PolytopeQuadrature& getQuadrature(
        size_t d, Index i, const QF::QuadratureFormula& formula) const override
      {
        return m_shard.getQuadrature(d, i, formula);
      }

      Mesh& load(const boost::filesystem::path&, IO::FileFormat) override
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

      CellIterator getCell() const override
      {
        return m_shard.getCell();
      }

      FaceIterator getFace() const override
      {
        return m_shard.getFace();
      }

      VertexIterator getVertex() const override
      {
        return m_shard.getVertex();
      }

      PolytopeIterator getPolytope(size_t dimension) const override
      {
        return m_shard.getPolytope(dimension);
      }

      CellIterator getCell(Index idx) const override
      {
        return m_shard.getCell(idx);
      }

      FaceIterator getFace(Index idx) const override
      {
        return m_shard.getFace(idx);
      }

      VertexIterator getVertex(Index idx) const override
      {
        return m_shard.getVertex(idx);
      }

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

      Real getVolume() const override
      {
        return getHost().getVolume();
      }

      Real getVolume(Attribute attribute) const override
      {
        return getHost().getVolume(attribute);
      }

      Real getVolume(const FlatSet<Attribute>& attributes) const override
      {
        return getHost().getVolume(attributes);
      }

      Real getArea() const override
      {
        return getHost().getArea();
      }

      Real getArea(Attribute attribute) const override
      {
        return getHost().getArea(attribute);
      }

      Real getArea(const FlatSet<Attribute>& attributes) const override
      {
        return getHost().getArea(attributes);
      }

      Real getMeasure(size_t d) const override
      {
        return getHost().getMeasure(d);
      }

      Real getMeasure(size_t d, Attribute attribute) const override
      {
        return getHost().getMeasure(d, attribute);
      }

      Real getMeasure(size_t d, const FlatSet<Attribute>& attributes) const override
      {
        return getHost().getMeasure(d, attributes);
      }

      Real getPerimeter() const override
      {
        return getHost().getPerimeter();
      }

      Real getPerimeter(Attribute attribute) const override
      {
        return getHost().getPerimeter(attribute);
      }

      Real getPerimeter(const FlatSet<Attribute>& attributes) const override
      {
        return getHost().getPerimeter(attributes);
      }

    private:
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

  /// @brief Deduces the distributed resident mesh of a host mesh.
  template <class Backend>
  Mesh(Rodin::Context::Kokkos<Rodin::Context::MPI, Backend>,
    const Mesh<Rodin::Context::MPI>&)
    -> Mesh<Rodin::Context::Kokkos<Rodin::Context::MPI, Backend>>;
}

#endif
