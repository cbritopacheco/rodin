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
   * The partitions are made by recursive bisection. A set of polytopes for k
   * partitions splits into a region for the first k / 2 partitions and the
   * rest for the others. The region starts from a polytope at the periphery
   * of the set and grows through the adjacency, lowest first along the axis
   * to the polytope furthest from it, until it has its exact share of the set.
   *
   * The sizes of the partitions differ by at most one polytope, and the
   * region of a bisection is connected on a connected set. A partition is
   * thus usually connected, though the rest of a bisection is not always.
   *
   * @see <a href="class_rodin_1_1_geometry_1_1_partitioner.html">Partitioner</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_greedy_partitioner.html">GreedyPartitioner</a>
   */
  class BalancedCompactPartitioner : public Partitioner
  {
    public:
      /// @brief Type of mesh used by this partitioner.
      using MeshType = Mesh<Context::Local>;

      /**
       * @brief Constructs a partitioner for the given mesh.
       * @param[in] mesh Mesh to partition
       */
      BalancedCompactPartitioner(const MeshType& mesh);

      /// @brief Virtual destructor.
      virtual ~BalancedCompactPartitioner() = default;

      /**
       * @brief Gets the mesh being partitioned.
       * @returns Reference to the mesh
       */
      virtual const MeshType& getMesh() const override;

      /**
       * @brief Partitions the mesh polytopes of dimension @p d.
       * @param[in] count Number of partitions
       * @param[in] d Dimension of polytopes to partition
       *
       * Requires the incidence (d, d). When @p count exceeds the number of
       * polytopes, some partitions are empty.
       */
      virtual void partition(size_t count, size_t d) override;

      /**
       * @brief Partitions the mesh cells.
       * @param[in] count Number of partitions
       */
      void partition(size_t count)
      {
        this->partition(count, getMesh().getDimension());
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
