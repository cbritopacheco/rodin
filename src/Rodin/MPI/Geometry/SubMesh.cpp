/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <vector>
#include <utility>

#include <boost/mpi/collectives.hpp>

#include "Rodin/Geometry/Polytope.h"

#include "SubMesh.h"

namespace Rodin::Geometry
{
  // ---------------------------------------------------------------------------
  // Builder
  // ---------------------------------------------------------------------------

  SubMesh<Context::MPI>::Builder&
  SubMesh<Context::MPI>::Builder::initialize(const Mesh<Context::MPI>& parent)
  {
    const auto& shard = parent.getShard();
    const size_t dim = shard.getDimension();

    m_parent = parent;

    m_shardBuilder.initialize(shard);

    m_s2ps.clear();
    m_s2ps.resize(dim + 1);

    m_sidx.assign(dim + 1, 0);
    m_dimension = 0;

    return *this;
  }

  SubMesh<Context::MPI>::Builder&
  SubMesh<Context::MPI>::Builder::include(size_t d, Index parentLocalIdx)
  {
    assert(m_parent.has_value());
    const auto& parentMesh = m_parent.value().get();
    const auto& parentShard = parentMesh.getShard();

    // Fast path: entity already included.
    if (m_s2ps[d].right.count(parentLocalIdx))
      return *this;

    // For d > 0: pre-include all vertices of this polytope.
    // Shard::Builder::include({d, idx}) for d > 0 requires that every vertex
    // of the polytope is already in the shard builder's vertex map.
    if (d > 0)
    {
      const auto& conn = parentShard.getConnectivity();
      const auto& verts = conn.getPolytope(d, parentLocalIdx);
      for (Index j = 0; j < static_cast<Index>(verts.size()); ++j)
      {
        const Index v = verts(j);
        if (m_s2ps[0].right.count(v))
          continue;

        const Shard::State vState = parentShard.getState(0)[v];
        const auto [subV, vInserted] = m_shardBuilder.include({0, v}, vState);

        if (vInserted)
        {
          m_s2ps[0].left.push_back(v);
          m_s2ps[0].right.emplace(v, subV);
          ++m_sidx[0];

          // Preserve ownership metadata from the parent shard.
          if (vState == Shard::State::Shared || vState == Shard::State::Ghost)
          {
            const auto& parentOwner = parentShard.getOwner(0);
            auto oit = parentOwner.find(v);
            if (oit != parentOwner.end())
              m_shardBuilder.setOwner(0, subV, oit->second);
          }
          else
          {
            const auto& parentHalo = parentShard.getHalo(0);
            auto hit = parentHalo.find(v);
            if (hit != parentHalo.end())
            {
              for (const Index neighborRank : hit->second)
                m_shardBuilder.halo(0, subV, neighborRank);
            }
          }
        }
      }
    }

    // Include the entity with its parent-shard state.
    const Shard::State state = parentShard.getState(d)[parentLocalIdx];
    const auto [subIdx, inserted] = m_shardBuilder.include({d, parentLocalIdx}, state);

    if (inserted)
    {
      m_s2ps[d].left.push_back(parentLocalIdx);
      m_s2ps[d].right.emplace(parentLocalIdx, subIdx);
      ++m_sidx[d];

      // Preserve ownership metadata.
      if (state == Shard::State::Shared || state == Shard::State::Ghost)
      {
        const auto& parentOwner = parentShard.getOwner(d);
        auto oit = parentOwner.find(parentLocalIdx);
        if (oit != parentOwner.end())
          m_shardBuilder.setOwner(d, subIdx, oit->second);
      }
      else
      {
        const auto& parentHalo = parentShard.getHalo(d);
        auto hit = parentHalo.find(parentLocalIdx);
        if (hit != parentHalo.end())
        {
          for (const Index neighborRank : hit->second)
            m_shardBuilder.halo(d, subIdx, neighborRank);
        }
      }

      m_dimension = std::max(m_dimension, d);
    }

    return *this;
  }

  SubMesh<Context::MPI> SubMesh<Context::MPI>::Builder::finalize()
  {
    assert(m_parent.has_value());
    const auto& parentMesh = m_parent.value().get();
    const auto& parentShard = parentMesh.getShard();

    Shard shard = m_shardBuilder.finalize();

    // Remap the shard PolytopeMap so that left[subLocal] holds the distributed
    // global index rather than the parent-shard-local index.
    //
    // After Shard::Builder::finalize() in Parent mode, left[subLocal] contains
    // the parent-shard-local index for each dimension. We need distributed
    // indices so that global operations (e.g. reconcile, MPI reductions) work
    // correctly.
    const size_t dim = shard.getDimension();
    for (size_t d = 0; d <= dim; ++d)
    {
      auto& pm = shard.getPolytopeMap(d);
      const auto& parentPM = parentShard.getPolytopeMap(d);

      UnorderedMap<Index, Index> newRight;
      newRight.reserve(pm.left.size());

      for (size_t i = 0; i < pm.left.size(); ++i)
      {
        const Index parentShardLocal = pm.left[i];
        const Index distributed = parentPM.left.at(parentShardLocal);
        pm.left[i] = distributed;
        newRight.emplace(distributed, static_cast<Index>(i));
      }

      pm.right = std::move(newRight);
    }

    // -------------------------------------------------------------------------
    // Orphan-ownership resolution.
    //
    // A Ghost/Shared entity in this shard can be "orphaned" when its stated
    // parent-shard owner did not include it in the SubMesh (e.g. a boundary
    // vertex whose parent owner has no boundary faces on its shard).  In that
    // case no rank will ever send a DOF for it, crashing or silently corrupting
    // FES constructors (P1, etc.) that rely on owner→ghost DOF exchange.
    //
    // Fix: for each dimension, do one all-gather of (distributedGlobalID, rank)
    // pairs for Owned entities.  Any Ghost/Shared entity whose stated owner is
    // absent from the gathered set is re-owned by the minimum rank that has it,
    // determined from a second all-gather of all (not just Owned) entity GIDs.
    // -------------------------------------------------------------------------
    {
      const auto& comm = parentMesh.getContext().getCommunicator();
      const int   rank = comm.rank();
      const int   P    = comm.size();

      for (size_t d = 0; d <= dim; ++d)
      {
        const auto& pm    = shard.getPolytopeMap(d);
        const auto& state = shard.getState(d);
        const size_t n    = pm.left.size();

        // Collect global IDs of Owned entities on this rank.
        std::vector<Index> myOwnedGids;
        myOwnedGids.reserve(n);
        for (size_t i = 0; i < n; ++i)
          if (state[i] == Shard::State::Owned)
            myOwnedGids.push_back(pm.left[i]);

        // All-gather: find out which global IDs are owned by at least one rank.
        std::vector<std::vector<Index>> allOwnedGids;
        boost::mpi::all_gather(comm, myOwnedGids, allOwnedGids);

        UnorderedSet<Index> globallyOwned;
        for (const auto& v : allOwnedGids)
          for (Index gid : v)
            globallyOwned.insert(gid);

        // Check whether any Ghost/Shared entity on this rank is orphaned.
        bool hasOrphan = false;
        for (size_t i = 0; i < n && !hasOrphan; ++i)
          if (state[i] != Shard::State::Owned && !globallyOwned.count(pm.left[i]))
            hasOrphan = true;

        bool anyOrphan = false;
        boost::mpi::all_reduce(comm, hasOrphan, anyOrphan, std::logical_or<bool>());
        if (!anyOrphan)
          continue;

        // Second all-gather: for each rank, the list of ALL entity GIDs
        // present in its SubMesh shard (regardless of state).
        std::vector<Index> myAllGids;
        myAllGids.reserve(n);
        for (size_t i = 0; i < n; ++i)
          myAllGids.push_back(pm.left[i]);

        std::vector<std::vector<Index>> allAllGids;
        boost::mpi::all_gather(comm, myAllGids, allAllGids);

        // For each orphaned GID, the new owner is the minimum rank (in
        // communicator order) that has it in its SubMesh.
        UnorderedMap<Index, int> orphanNewOwner;
        for (int r = 0; r < P; ++r)
        {
          for (Index gid : allAllGids[r])
          {
            if (!globallyOwned.count(gid) && !orphanNewOwner.count(gid))
              orphanNewOwner.emplace(gid, r); // first (= minimum) rank wins
          }
        }

        if (orphanNewOwner.empty())
          continue;

        // Apply ownership changes to this rank's shard.
        auto& mutableState = shard.getState(d);
        auto& ownerMap     = shard.getOwner(d);
        auto& haloMap      = shard.getHalo(d);

        for (size_t i = 0; i < n; ++i)
        {
          const Index gid = pm.left[i];
          auto it = orphanNewOwner.find(gid);
          if (it == orphanNewOwner.end())
            continue;

          const int newOwner = it->second;

          if (newOwner == rank)
          {
            // This rank becomes the Owned entity.
            mutableState[i] = Shard::State::Owned;
            ownerMap.erase(i);
            // Register all other ranks that have this entity as halo members.
            for (int r = 0; r < P; ++r)
            {
              if (r == rank)
                continue;
              for (Index otherGid : allAllGids[r])
              {
                if (otherGid == gid)
                {
                  haloMap[i].insert(static_cast<Index>(r));
                  break;
                }
              }
            }
          }
          else
          {
            // Update the stated owner to the new owner rank.
            ownerMap[i] = static_cast<Index>(newOwner);
          }
        }
      }
    }

    Mesh<Context::MPI>::Builder meshBuilder(parentMesh.getContext());
    meshBuilder.initialize(std::move(shard));

    SubMesh result(parentMesh);
    result.Parent::operator=(meshBuilder.finalize());
    result.m_s2ps = std::move(m_s2ps);
    return result;
  }

  // ---------------------------------------------------------------------------
  // SubMesh<Context::MPI>
  // ---------------------------------------------------------------------------

  SubMesh<Context::MPI>::SubMesh(
      std::reference_wrapper<const Mesh<Context::MPI>> parent)
    : Parent(parent.get().getContext()),
      m_parent(parent)
  {
    if (parent.get().isSubMesh())
    {
      const auto& submesh = parent.get().asSubMesh();
      m_ancestors = submesh.getAncestors();
    }
    m_ancestors.push_front(parent.get());
  }

  SubMesh<Context::MPI>::SubMesh(const SubMesh& other)
    : Parent(other),
      m_parent(other.m_parent),
      m_s2ps(other.m_s2ps),
      m_ancestors(other.m_ancestors)
  {}

  SubMesh<Context::MPI>::SubMesh(SubMesh&& other)
    : Parent(std::move(other)),
      m_parent(std::move(other.m_parent)),
      m_s2ps(std::move(other.m_s2ps)),
      m_ancestors(std::move(other.m_ancestors))
  {}

  const Mesh<Context::MPI>& SubMesh<Context::MPI>::getParent() const
  {
    return m_parent.get();
  }

  bool SubMesh<Context::MPI>::isLocalPoint(const Point& p) const
  {
    const auto& mesh = p.getPolytope().getMesh();
    return mesh == static_cast<const MeshBase&>(*this) || mesh == this->getShard();
  }

  Optional<Point> SubMesh<Context::MPI>::restriction(const Point& p) const
  {
    const auto& polytope = p.getPolytope();
    const size_t d = polytope.getDimension();
    Index i = polytope.getIndex();

    const auto& ancestors = getAncestors();

    Deque<std::reference_wrapper<const SubMeshBase>> descendants;
    descendants.push_front(*this);

    for (auto it = ancestors.begin(); it != ancestors.end(); ++it)
    {
      if (it->get() == polytope.getMesh())
      {
        break;
      }
      else if (!it->get().isSubMesh())
      {
        // Invalid restriction: the SubMesh is not a descendant of the mesh
        // to which the Point belongs.
        return {};
      }
      else
      {
        descendants.push_front(it->get().asSubMesh());
      }
    }

    for (auto it = descendants.begin(); it != descendants.end(); ++it)
    {
      const auto& polytopeMap = it->get().getPolytopeMap(d);
      auto find = polytopeMap.right.find(i);
      if (find == polytopeMap.right.end())
      {
        // Polytope (d, i) is not present in this submesh.
        return {};
      }
      i = find->second;
    }

    return Point(*getPolytope(d, i), p.getReferenceCoordinates(), p.getPhysicalCoordinates());
  }
}
