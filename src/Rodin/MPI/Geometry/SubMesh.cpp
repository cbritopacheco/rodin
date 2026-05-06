/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <vector>
#include <utility>
#include <algorithm>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
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
    // Shared-vertex ownership resolution for dimension 0 (vertices).
    //
    // The SubMesh builder copies vertex ownership verbatim from the parent
    // volume-partition shard.  This is wrong for SubMeshes built from a
    // strict subset of polytopes (e.g. boundary faces only): a vertex v with
    // state Shared(owner=A) in the volume partition may have NO owned boundary
    // face on A's shard, so A never includes v in its SubMesh.  No rank then
    // owns v → the DOF exchange in P1 (and similar FE spaces) breaks.
    //
    // Fix — 2-round neighbor exchange:
    //
    //   Round 1  Each rank sends the global IDs of its SubMesh-Shared vertices
    //            to their stated owner, using the parent halo/owner maps for
    //            symmetric communication (both owner rank and halo ranks appear
    //            in each other's neighbor sets).
    //
    //   Round 2  The stated owner replies per GID:
    //            - empty list     → "I do include it as Owned; keep Shared."
    //            - sorted querier list → "I do not include it; min querier
    //              becomes the new Owned, others redirect their owner pointer."
    //
    // Communication cost: O(shared-boundary-vertices / P) per rank, which is
    // the same as the subsequent P1 DOF exchange.  No all-gather is needed.
    // -------------------------------------------------------------------------
    {
      const auto& comm = parentMesh.getContext().getCommunicator();
      const int   rank = comm.rank();

      // Only vertices (d == 0) can be Shared in an owned-face-driven SubMesh.
      // Faces are always Owned (getBoundary() returns only owned faces).
      const size_t d = 0;

      if (dim == 0)
        goto after_ownership_fix; // 0-D mesh: no edges/faces, nothing to do

      {
        const auto& pm    = shard.getPolytopeMap(d);
        const size_t n    = pm.left.size();

        auto& subState = shard.getState(d);
        auto& subOwner = shard.getOwner(d);
        auto& subHalo  = shard.getHalo(d);

        // Build the neighbor set from the PARENT shard so that the send/recv
        // pattern is symmetric across all ranks without global communication.
        // Rationale: if v is Shared(owner=A) on B in the parent, A has B in its
        // parent halo — both A and B therefore appear in each other's set.
        const auto& pOwner = parentShard.getOwner(d);
        const auto& pHalo  = parentShard.getHalo(d);

        UnorderedSet<int> neighborSet;
        for (const auto& [li, r] : pOwner)
          neighborSet.insert(static_cast<int>(r));
        for (const auto& [li, rs] : pHalo)
          for (Index r : rs)
            neighborSet.insert(static_cast<int>(r));

        if (neighborSet.empty())
          goto after_ownership_fix;

        // queryMap[ownerRank] = GIDs of SubMesh-Shared vertices with that owner.
        // Initialise all neighbors to empty so we always send (even if nothing
        // to ask), ensuring the matching irecv on the other side is satisfied.
        UnorderedMap<int, std::vector<Index>> queryMap;
        for (int r : neighborSet)
          queryMap[r]; // default-construct empty vector

        for (size_t i = 0; i < n; ++i)
        {
          if (subState[i] != Shard::State::Shared)
            continue;
          auto oit = subOwner.find(i);
          if (oit != subOwner.end())
            queryMap[static_cast<int>(oit->second)].push_back(pm.left[i]);
        }

        // ── Round 1: send queries, receive queries ───────────────────────────
        std::vector<boost::mpi::request> reqs;
        UnorderedMap<int, std::vector<Index>> recvQuery;
        for (int r : neighborSet)
        {
          recvQuery[r]; // default-construct
          reqs.push_back(comm.irecv(r, 10, recvQuery[r]));
          reqs.push_back(comm.isend(r, 10, queryMap[r]));
        }
        boost::mpi::wait_all(reqs.begin(), reqs.end());

        // ── Build replies ────────────────────────────────────────────────────
        // For each GID queried of us: if it is Owned in our SubMesh → reply
        // with empty querier list.  Otherwise collect all queriers and reply
        // with the sorted list (each querier computes the min independently).

        // Build owned-GID → sub-local-index lookup.
        UnorderedMap<Index, Index> ownedGidIdx;
        for (size_t i = 0; i < n; ++i)
          if (subState[i] == Shard::State::Owned)
            ownedGidIdx.emplace(pm.left[i], static_cast<Index>(i));

        // Aggregate queriers per GID.
        UnorderedMap<Index, std::vector<int>> gidQueriers;
        for (auto& [r, gids] : recvQuery)
          for (Index gid : gids)
            gidQueriers[gid].push_back(r);

        // ── Prune stale halo entries for Owned vertices ───────────────────────
        // The SubMesh builder copies the parent shard's halo map verbatim.
        // Some entries may reference ranks that hold the vertex in the volume
        // partition but did NOT include it in this SubMesh (e.g. they have no
        // boundary face adjacent to it).  Those ranks will not participate in
        // the subsequent P1 / P1-vector DOF exchange, so if an Owned vertex's
        // halo still names them the P1 constructor posts an irecv/isend for a
        // rank that never posts the matching operation — causing a deadlock.
        //
        // After Round 1 we know exactly which ranks queried about each GID:
        // only ranks that actually included the vertex as Shared in their
        // SubMesh send a query.  We therefore prune every halo entry that
        // corresponds to a rank that sent no query for that vertex.
        for (size_t i = 0; i < n; ++i)
        {
          if (subState[i] != Shard::State::Owned)
            continue;
          auto haloIt = subHalo.find(i);
          if (haloIt == subHalo.end())
            continue;

          const Index gid = pm.left[i];
          IndexSet& haloSet = haloIt->second;

          auto qIt = gidQueriers.find(gid);
          if (qIt == gidQueriers.end())
          {
            // No rank queried this vertex → no rank has it as Shared in SubMesh.
            subHalo.erase(haloIt);
          }
          else
          {
            const UnorderedSet<int> querierSet(qIt->second.begin(), qIt->second.end());
            IndexSet newHalo;
            for (const Index r : haloSet)
              if (querierSet.count(static_cast<int>(r)))
                newHalo.insert(r);
            if (newHalo.empty())
              subHalo.erase(haloIt);
            else
              haloSet = std::move(newHalo);
          }
        }

        // Reply type: { gid, querier_list }.
        // Empty querier_list = "I own it."
        // Non-empty sorted querier_list = "I don't; min rank should own it."
        using GidInfo = std::pair<Index, std::vector<int>>;
        UnorderedMap<int, std::vector<GidInfo>> replyMap;
        for (int r : neighborSet)
          replyMap[r]; // default-construct

        for (auto& [gid, queriers] : gidQueriers)
        {
          if (ownedGidIdx.count(gid))
          {
            for (int q : queriers)
              replyMap[q].push_back({gid, {}});
          }
          else
          {
            std::sort(queriers.begin(), queriers.end());
            for (int q : queriers)
              replyMap[q].push_back({gid, queriers});
          }
        }

        // ── Round 2: send replies, receive replies ───────────────────────────
        reqs.clear();
        UnorderedMap<int, std::vector<GidInfo>> recvReply;
        for (int r : neighborSet)
        {
          recvReply[r]; // default-construct
          reqs.push_back(comm.irecv(r, 11, recvReply[r]));
          reqs.push_back(comm.isend(r, 11, replyMap[r]));
        }
        boost::mpi::wait_all(reqs.begin(), reqs.end());

        // ── Apply ownership corrections ──────────────────────────────────────
        for (auto& [r, replies] : recvReply)
        {
          for (auto& [gid, queriers] : replies)
          {
            if (queriers.empty())
              continue; // stated owner does include it; keep Shared(owner=r)

            // Find the sub-local index for this GID.
            auto pit = pm.right.find(gid);
            if (pit == pm.right.end())
              continue; // should not happen

            const Index localIdx = static_cast<Index>(pit->second);
            const int   newOwner = queriers.front(); // already sorted, front = min

            if (newOwner == rank)
            {
              // This rank becomes the SubMesh owner of the vertex.
              subState[localIdx] = Shard::State::Owned;
              subOwner.erase(localIdx);
              // Halo: all other queriers need the DOF from us in P1.
              for (int q : queriers)
                if (q != rank)
                  subHalo[localIdx].insert(static_cast<Index>(q));
            }
            else
            {
              // Redirect the owner pointer from r (old owner) to newOwner.
              subOwner[localIdx] = static_cast<Index>(newOwner);
            }
          }
        }
      }
      after_ownership_fix:;
    } // ownership resolution block

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
