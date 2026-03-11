/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_BALANCEDCOMPACTPARTITIONER_H
#define RODIN_GEOMETRY_BALANCEDCOMPACTPARTITIONER_H

/**
 * @file
 * @brief Balanced compact partitioning strategy for mesh decomposition.
 */

#include "MeshPartitioner.h"

namespace Rodin::Geometry
{
  /**
   * @brief Partitioner that creates balanced, spatially compact partitions.
   *
   * This partitioning strategy divides the mesh into partitions of
   * approximately equal size while maintaining spatial locality. The
   * algorithm aims to minimize the interface between partitions.
   *
   * @note This is particularly useful for parallel computing where
   * minimizing communication between processors is important.
   *
   * @see Partitioner, GreedyPartitioner
   */
  class BalancedCompactPartitioner : public Partitioner
  {
    public:
      /**
       * @brief Type of mesh used by this partitioner.
       */
      using MeshType = Mesh<Context::Local>;

      /**
       * @brief Constructs a partitioner for the given mesh.
       * @param[in] mesh Mesh to partition
       */
      BalancedCompactPartitioner(const MeshType& mesh);

      /**
       * @brief Virtual destructor.
       */
      virtual ~BalancedCompactPartitioner() = default;

      /**
       * @brief Gets the mesh being partitioned.
       * @returns Reference to the mesh
       */
      virtual const MeshType& getMesh() const override;

      /**
       * @brief Partitions the mesh polytopes of dimension @p d.
       * @param[in] maxPartitionSize Maximum number of polytopes per partition
       * @param[in] d Dimension of polytopes to partition
       *
       * Creates partitions such that no partition contains more than
       * @p maxPartitionSize polytopes, while maintaining spatial compactness.
       */
      virtual void partition(size_t maxPartitionSize, size_t d) override;

      /**
       * @brief Partitions the mesh cells.
       * @param[in] maxPartitionSize Maximum number of cells per partition
       *
       * Convenience method that partitions the highest-dimensional polytopes
       * (cells) of the mesh.
       */
      void partition(size_t maxPartitionSize)
      {
        this->partition(maxPartitionSize, getMesh().getDimension());
      }

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
