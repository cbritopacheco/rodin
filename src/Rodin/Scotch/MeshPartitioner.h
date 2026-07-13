/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MeshPartitioner.h
 * @brief Scotch-based mesh partitioner.
 *
 * Implements the @ref Rodin::Geometry "Geometry" partitioner interface
 * using the SCOTCH graph partitioning library, for distributing meshes
 * across processes.
 */
#ifndef RODIN_SCOTCH_MESHPARTITIONER_H
#define RODIN_SCOTCH_MESHPARTITIONER_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/MeshPartitioner.h"

#include <scotch.h>

namespace Rodin::Scotch
{
  class Partitioner : public Geometry::Partitioner
  {
    public:
      /// @brief Mesh type.
      using MeshType = Geometry::Mesh<Context::Local>;

      Partitioner(const MeshType& mesh);

      ~Partitioner() override;

      const MeshType& getMesh() const override;

      void partition(size_t count)
      {
        partition(count, getMesh().getDimension());
      }

      void partition(size_t numPartitions, size_t d) override;

      size_t getPartition(Index index) const override;

      Partitioner& setStrategy(const SCOTCH_Strat& strat);

      size_t getCount() const override
      {
        return m_numPartitions;
      }

    private:
      size_t m_numPartitions;
      std::reference_wrapper<const MeshType> m_mesh;
      std::vector<SCOTCH_Num> m_partition;

      SCOTCH_Graph m_graph;
      Optional<SCOTCH_Strat> m_strat;
  };
}

#endif
