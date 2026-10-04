/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or copy at https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_LOCATION_AABB_H
#define RODIN_MPI_LOCATION_AABB_H

#include <type_traits>
#include "Rodin/Location/AABB.h"
#include "Rodin/MPI/Geometry/Mesh.h"

namespace Rodin::Location
{
  /**
   * @brief Noncollective point location in owned entities of an MPI mesh.
   *
   * @section MPIAABBArchitecture Architecture
   * The MPI layer snapshots owned shard-local indices in every dimension.
   * A composed local AABB builds bounds and performs inversion on that subset.
   * Successful points are attached to the MPI mesh with unchanged dimension,
   * local index and coordinates. The MPI mesh delegates their transformations
   * to the same shard. Construction, setters and queries perform no communication.
   * This specialization also supports derived MPI meshes, including SubMesh.
   *
   * An empty result is a rank-local miss. Different ranks may each return an
   * owned incident entity at a shared spatial boundary; there is no global
   * arbitration. Nonowned support vertices remain available to owned entities.
   * Physical tolerance is relative to the full shard's vertex-box diagonal,
   * so its scale can vary with the partition and overlap. Reconstruct the
   * locator after geometry, topology, ownership or reconciliation changes.
   * The mesh must outlive the locator and its returned points.
   */
  template <class MeshType>
    requires std::is_base_of_v<Geometry::Mesh<Context::MPI>, MeshType>
  class AABB<MeshType> final
  {
    public:
      explicit AABB(const MeshType& mesh)
        : m_mesh(mesh),
          m_shard(mesh.getShard(), ownedCandidates(mesh))
      {}

      Real getTolerance() const
      {
        return m_shard.getTolerance();
      }
      Real getReferenceTolerance() const
      {
        return m_shard.getReferenceTolerance();
      }

      AABB& setTolerance(Real tolerance)
      {
        m_shard.setTolerance(tolerance);
        return *this;
      }

      AABB& setReferenceTolerance(Real tolerance)
      {
        m_shard.setReferenceTolerance(tolerance);
        return *this;
      }

      AABB& setExhaustiveFallback(bool enabled)
      {
        m_shard.setExhaustiveFallback(enabled);
        return *this;
      }

      AABB& setProjectionPruning(bool enabled)
      {
        m_shard.setProjectionPruning(enabled);
        return *this;
      }

      /// Searches owned entities of dimension d and lifts the local result.
      Optional<Geometry::Point> locate(
        size_t dimension, const Math::SpatialPoint& x) const
      {
        auto hit = m_shard.locate(dimension, x);
        if (!hit)
          return {};
        const auto& p = hit->getPolytope();
        return Geometry::Point(
          Geometry::Polytope(p.getDimension(), p.getIndex(), m_mesh.get()),
          hit->getReferenceCoordinates(), hit->getPhysicalCoordinates());
      }

      /// Uses the MPI mesh's logical dimension, including on empty submesh ranks.
      Optional<Geometry::Point> locate(const Math::SpatialPoint& x) const
      {
        return locate(m_mesh.get().getDimension(), x);
      }

    private:
      static typename AABB<Geometry::LocalMesh>::Candidates ownedCandidates(
        const MeshType& mesh)
      {
        const auto& shard = mesh.getShard();
        typename AABB<Geometry::LocalMesh>::Candidates candidates(
          mesh.getDimension() + 1);
        for (size_t d = 0; d <= shard.getDimension(); ++d)
        {
          for (Index i = 0; i < shard.getPolytopeCount(d); ++i)
            if (shard.isOwned(d, i))
              candidates[d].push_back(i);
        }
        return candidates;
      }

      std::reference_wrapper<const MeshType> m_mesh;
      AABB<Geometry::LocalMesh> m_shard;
  };
}
#endif
