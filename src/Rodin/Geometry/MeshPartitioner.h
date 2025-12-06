/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MESHPARTITIONER_H
#define RODIN_GEOMETRY_MESHPARTITIONER_H

/**
 * @file
 * @brief Base interface for mesh partitioning algorithms.
 */

#include "Mesh.h"

namespace Rodin::Geometry
{
  /**
   * @brief Abstract base class for mesh partitioning strategies.
   *
   * A Partitioner divides mesh polytopes into disjoint partitions, which is
   * useful for load balancing, parallel processing, and domain decomposition
   * methods. Each polytope is assigned to exactly one partition.
   *
   * @see BalancedCompactPartitioner, GreedyPartitioner
   */
  class Partitioner
  {
    public:
      /**
       * @brief Default constructor.
       */
      Partitioner() = default;

      /**
       * @brief Virtual destructor.
       */
      virtual ~Partitioner() = default;

      /**
       * @brief Gets the mesh being partitioned.
       * @returns Reference to the mesh
       */
      virtual const Mesh<Context::Local>& getMesh() const = 0;

      /**
       * @brief Performs the partitioning of mesh polytopes.
       * @param[in] numPartitions Number of partitions to create
       * @param[in] d Dimension of polytopes to partition
       *
       * Divides the polytopes of dimension @p d into @p numPartitions
       * disjoint partitions according to the partitioning strategy.
       */
      virtual void partition(size_t numPartitions, size_t d) = 0;

      /**
       * @brief Gets the partition index for a given polytope.
       * @param[in] index Index of the polytope
       * @returns Partition number assigned to the polytope
       */
      virtual size_t getPartition(Index index) const = 0;

      /**
       * @brief Operator overload for getting partition index.
       * @param[in] index Index of the polytope
       * @returns Partition number assigned to the polytope
       */
      size_t operator[](Index index) const
      {
        return getPartition(index);
      }

      /**
       * @brief Gets the total number of partitions.
       * @returns Number of partitions created
       */
      virtual size_t getCount() const = 0;
  };
}

#endif
