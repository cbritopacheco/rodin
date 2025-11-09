/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_GREEDYPARTITIONER_H
#define RODIN_GEOMETRY_GREEDYPARTITIONER_H

/**
 * @file
 * @brief Greedy partitioning strategy for mesh decomposition.
 */

#include "MeshPartitioner.h"

namespace Rodin::Geometry
{
  /**
   * @brief Partitioner using a greedy algorithm.
   *
   * This partitioning strategy uses a greedy approach to divide the mesh
   * into a specified number of partitions. The algorithm processes polytopes
   * sequentially and assigns them to partitions to balance the load.
   *
   * @note This is a simpler and faster alternative to BalancedCompactPartitioner,
   * but may produce less spatially compact partitions.
   *
   * @see Partitioner, BalancedCompactPartitioner
   */
  class GreedyPartitioner : public Partitioner
  {
    public:
      /**
       * @brief Type of mesh used by this partitioner.
       */
      using MeshType = Geometry::Mesh<Context::Local>;

      /**
       * @brief Constructs a greedy partitioner for the given mesh.
       * @param[in] mesh Mesh to partition
       */
      GreedyPartitioner(const MeshType& mesh);

      /**
       * @brief Virtual destructor.
       */
      virtual ~GreedyPartitioner() = default;

      /**
       * @brief Gets the mesh being partitioned.
       * @returns Reference to the mesh
       */
      virtual const MeshType& getMesh() const override;

      /**
       * @brief Partitions the mesh cells.
       * @param[in] count Number of partitions to create
       *
       * Convenience method that partitions the highest-dimensional polytopes
       * (cells) of the mesh into @p count partitions.
       */
      void partition(size_t count)
      {
        partition(count, getMesh().getDimension());
      }

      /**
       * @brief Partitions the mesh polytopes of dimension @p d.
       * @param[in] numPartitions Number of partitions to create
       * @param[in] d Dimension of polytopes to partition
       *
       * Uses a greedy algorithm to divide polytopes of dimension @p d into
       * @p numPartitions approximately equal-sized partitions.
       */
      virtual void partition(size_t numPartitions, size_t d) override;

      /**
       * @brief Gets the partition index for a given polytope.
       * @param[in] index Index of the polytope
       * @returns Partition number assigned to the polytope
       */
      virtual size_t getPartition(Index index) const override;

      /**
       * @brief Gets the total number of partitions.
       * @returns Number of partitions created
       */
      size_t getCount() const override;

    private:
      size_t m_count;
      std::reference_wrapper<const MeshType> m_mesh;
      std::vector<size_t> m_partition;
  };
}

#endif
