/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SHARD_H
#define RODIN_GEOMETRY_SHARD_H

/**
 * @file
 * @brief Simple mesh shard for mesh decomposition (deprecated/legacy).
 *
 * @note This file defines a simplified Shard class. For the full-featured
 * implementation with ghost/owned polytope tracking, see Shard.h.
 */

#include "Mesh.h"

namespace Rodin::Geometry
{
  /**
   * @brief Simplified mesh shard representing a subset of a mesh.
   *
   * A Shard is a mesh that represents a portion of a parent mesh, typically
   * used for domain decomposition or distributed computing. This is a
   * simplified version that maintains bidirectional mappings between shard
   * and parent polytope indices.
   *
   * # Relationship to Parent Mesh
   *
   * A Shard maintains a reference relationship to its parent mesh through
   * index mappings. Each polytope dimension has a separate bidirectional
   * map (boost::bimap) that allows efficient lookup in both directions.
   *
   * # Use Cases
   *
   * - Domain decomposition for parallel computing
   * - Mesh partitioning for load balancing
   * - Extracting mesh regions for specialized processing
   *
   * @note For distributed MPI applications with ghost polytope tracking,
   * consider using the full Shard class defined in Shard.h instead.
   *
   * @see Sharder, Mesh, SubMesh
   */
  class Shard final : public Mesh<Context::Local>
  {
    public:
      /**
       * @brief Parent mesh type.
       */
      using Parent = Mesh<Rodin::Context::Local>;
      
      /**
       * @brief Context type.
       */
      using Context = typename Parent::Context;

      /**
       * @brief Builder for constructing Shard instances.
       *
       * Provides a fluent interface for incrementally building a shard from
       * a parent mesh by selecting specific polytopes.
       */
      class Builder
      {
        public:
          /**
           * @brief Default constructor.
           */
          Builder() = default;

          /**
           * @brief Initializes the builder with a parent mesh.
           * @param[in] parent Parent mesh from which to extract the shard
           * @returns Reference to this builder for method chaining
           */
          Builder& initialize(const Mesh<Context::Local>& parent);

          /**
           * @brief Includes a single polytope in the shard.
           * @param[in] d Dimension of the polytope
           * @param[in] parentIdx Index of the polytope in the parent mesh
           * @returns Reference to this builder for method chaining
           */
          Builder& include(size_t d, Index parentIdx);

          /**
           * @brief Includes multiple polytopes in the shard.
           * @param[in] d Dimension of the polytopes
           * @param[in] indices Set of polytope indices in the parent mesh
           * @returns Reference to this builder for method chaining
           */
          Builder& include(size_t d, const IndexSet& indices);

          /**
           * @brief Finalizes construction and returns the built shard.
           * @returns Newly constructed Shard object
           */
          Shard finalize();

        private:
          Optional<std::reference_wrapper<const Mesh<Context>>> m_parent;
          Mesh<Context>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<boost::bimap<Index, Index>> m_s2ps;
          size_t m_dimension;
      };

      /**
       * @brief Default constructor.
       */
      Shard();

      /**
       * @brief Copy constructor.
       * @param[in] other Shard to copy from
       */
      Shard(const Shard& other);

      /**
       * @brief Move constructor.
       * @param[in] other Shard to move from
       */
      Shard(Shard&& other);

      /**
       * @brief Copy assignment operator (deleted).
       */
      Shard& operator=(const Shard&) = delete;

      /**
       * @brief Move assignment operator.
       * @param[in] other Shard to move from
       * @returns Reference to this object
       */
      Shard& operator=(Shard&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_s2ps = std::move(other.m_s2ps);
        }
        return *this;
      }

      /**
       * @brief Gets the bidirectional polytope index map for a dimension.
       * @param[in] d Dimension of polytopes
       * @returns Bidirectional map between shard and parent indices
       *
       * The returned bimap allows efficient lookup in both directions:
       * - left: shard index -> parent index
       * - right: parent index -> shard index
       */
      const boost::bimap<Index, Index>& getPolytopeMap(size_t d) const
      {
        return m_s2ps.at(d);
      }

    private:
      std::vector<boost::bimap<Index, Index>> m_s2ps; ///< Shard-to-parent index mappings per dimension
  };
}

#endif
