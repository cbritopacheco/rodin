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

#include <string>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/MeshPartitioner.h"

#include <scotch.h>

namespace Rodin::Scotch
{
  /**
   * @brief Mesh partitioner implemented with the SCOTCH library.
   *
   * Partitions a local mesh by building the graph representation expected by
   * SCOTCH and storing the resulting partition index for each mesh entity.
   */
  class Partitioner : public Geometry::Partitioner
  {
    public:
      /// @brief Mesh type.
      using MeshType = Geometry::Mesh<Context::Local>;

      /**
       * @brief Constructs a partitioner for a local mesh.
       * @param mesh Mesh to partition.
       */
      Partitioner(const MeshType& mesh);

      ~Partitioner() override;

      const MeshType& getMesh() const override;

      /**
       * @brief Partitions the mesh into a given number of parts.
       * @param count Number of partitions to create.
       */
      void partition(size_t count)
      {
        partition(count, getMesh().getDimension());
      }

      void partition(size_t numPartitions, size_t d) override;

      /**
       * @brief Gets a stored partition assignment.
       * @param index Entity index in the most recently partitioned dimension.
       * @return Partition index for the selected entity.
       * @throws Alert::Exception If no assignment exists for @p index.
       */
      size_t getPartition(Index index) const override;

      /**
       * @brief Sets the SCOTCH partitioning strategy.
       * Strategy contents are copied; the caller retains ownership of @p strat
       * and may release it immediately after this call.
       * @param strat Graph mapping strategy object used by SCOTCH.
       * @return Reference to this partitioner.
       */
      Partitioner& setStrategy(const SCOTCH_Strat& strat);

      size_t getCount() const override
      {
        return m_numPartitions;
      }

    private:
      /// @brief Requested count of the most recent successful partition operation.
      size_t m_numPartitions;

      /// @brief Mesh whose incidence graph is partitioned.
      std::reference_wrapper<const MeshType> m_mesh;

      /// @brief Most recent successful assignment, indexed by mesh entity.
      std::vector<SCOTCH_Num> m_partition;

      /// @brief Independent serialized mapping strategy, absent for the Scotch default.
      Optional<std::string> m_strategy;
  };
}

#endif
