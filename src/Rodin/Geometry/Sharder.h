/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SHARDER_H
#define RODIN_GEOMETRY_SHARDER_H

/**
 * @file
 * @brief Sharder interface for mesh decomposition.
 */

#include "Rodin/Geometry/MeshPartitioner.h"

#include "Shard.h"

namespace Rodin::Geometry
{
  template <class Context>
  class SharderBase
  {
    public:
      using ContextType = Context;

      /**
       * @brief Construct a SharderBase with the given context.
       * @param context The context.
       */
      SharderBase(const Context& context)
        : m_context(context)
      {};

      SharderBase& shard(Partitioner& partitioner)
      {
        m_shards.clear();

        const auto& mesh = partitioner.getMesh();
        const auto& conn = mesh.getConnectivity();
        const size_t cellDim = mesh.getDimension();
        const size_t numShards = partitioner.getCount();

        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, cellDim, cellDim);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, cellDim, 0);

        assert(numShards > 0);

        std::vector<Shard::Builder> sbs(numShards);
        for (auto& sb : sbs)
          sb.initialize(mesh);

        std::vector<std::vector<Index>> owner(cellDim + 1);     // Owner of each polytope
        std::vector<std::vector<Index>> local(cellDim + 1);     // Owning shard index for each polytope
        std::vector<std::vector<IndexSet>> halo(cellDim + 1);   // Halo information for each polytope
        std::vector<std::vector<uint8_t>> visited(cellDim + 1); // Visited polytopes

        for (size_t d = 0; d < cellDim + 1; d++)
        {
          owner[d].resize(mesh.getPolytopeCount(d));
          local[d].resize(mesh.getPolytopeCount(d));
          halo[d].resize(mesh.getPolytopeCount(d));
          visited[d].resize(mesh.getPolytopeCount(d), false);
        }

        for (Index cellIdx = 0; cellIdx < mesh.getCellCount(); cellIdx++)
        {
          const size_t partition = partitioner.getPartition(cellIdx);
          for (size_t d = 0; d <= cellDim - 1; d++)
          {
            const auto& inc = conn.getIncidence(cellDim, d);
            if (inc.empty())
              continue;
            for (const Index& idx : inc.at(cellIdx))
            {
              if (visited[d][idx])
              {
                if (owner[d][idx] != partition)
                  halo[d][idx].insert(partition);
                const auto [localIdx, inserted] =
                  sbs[partition].include({ d, idx }, Shard::State::Shared);
                if (inserted)
                  sbs[partition].getOwner(d)[localIdx] = owner[d][idx];
              }
              else
              {
                const auto [localIdx, inserted] =
                  sbs[partition].include({ d, idx }, Shard::State::Owned);
                assert(inserted);
                owner[d][idx] = partition;
                local[d][idx] = localIdx;
                visited[d][idx] = true;
              }
            }
          }
          assert(!visited[cellDim][cellIdx]);
          const auto [localIdx, inserted] =
            sbs[partition].include({ cellDim, cellIdx }, Shard::State::Owned);
          assert(inserted);
          owner[cellDim][cellIdx] = partition;
          local[cellDim][cellIdx] = localIdx;
          visited[cellDim][cellIdx] = true;
        }

        for (Index cellIdx = 0; cellIdx < mesh.getCellCount(); cellIdx++)
        {
          const size_t partition = partitioner.getPartition(cellIdx);
          for (const Index& nbr : conn.getIncidence({ cellDim, cellDim }, cellIdx))
          {
            if (partition == partitioner.getPartition(nbr))
              continue;
            for (size_t d = 0; d <= cellDim - 1; d++)
            {
              const auto& inc = conn.getIncidence(cellDim, d);
              if (inc.empty())
                continue;
              for (const Index idx : inc.at(nbr))
              {
                auto find = sbs[partition].getPolytopeMap(d).right.find(idx);
                if (find == sbs[partition].getPolytopeMap(d).right.end())
                {
                  halo[d][idx].insert(partition);
                  const auto [localIdx, inserted] =
                    sbs[partition].include({ d, idx }, Shard::State::Ghost);
                  if (inserted)
                    sbs[partition].getOwner(d)[localIdx] = owner[d][idx];
                }
              }
            }
            halo[cellDim][nbr].insert(partition);
            const auto [localIdx, inserted] =
              sbs[partition].include({ cellDim, nbr }, Shard::State::Ghost);
            if (inserted)
              sbs[partition].getOwner(cellDim)[localIdx] = owner[cellDim][nbr];
          }
        }

        // Now move the halo information into the shards
        for (size_t d = 0; d < cellDim + 1; d++)
        {
          for (size_t i = 0; i < mesh.getPolytopeCount(d); i++)
          {
            if (!visited[d][i])
              continue;
            auto& hm = sbs[owner[d][i]].getHalo(d);
            const size_t k = local[d][i];
            auto [it, inserted] = hm.try_emplace(k, std::move(halo[d][i]));
            if (!inserted)
              it->second.insert(halo[d][i].begin(), halo[d][i].end());
          }
        }

        m_shards.resize(numShards);

#pragma omp parallel for schedule(dynamic, 1)
        for (size_t i = 0; i < numShards; i++)
          m_shards[i] = sbs[i].finalize();

        return *this;
      }

      auto& getShards()
      {
        return m_shards;
      }

      const auto& getShards() const
      {
        return m_shards;
      }

      const ContextType& getContext() const
      {
        return m_context;
      }

    private:
      std::vector<Shard> m_shards;
      ContextType m_context;
  };

  /**
   * @brief Interface for mesh sharding algorithms.
   *
   * A Sharder decomposes a mesh into smaller pieces (shards) that can be
   * processed independently. This is useful for distributed computing and
   * parallel finite element analysis.
   *
   * @tparam Context The context type (e.g., Context::Local, Context::MPI)
   *
   * @see Shard, MeshShard
   */
  template <class Context>
  class Sharder;
}

#endif
