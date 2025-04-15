/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MPI_MESH_H
#define RODIN_GEOMETRY_MPI_MESH_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_MPI

#include "Rodin/Context/MPI.h"

#include "MPIConnectivity.h"

#include "Mesh.h"
#include "Shard.h"


namespace Rodin::Geometry
{
  using MPIMesh = Mesh<Context::MPI>;

  template <>
  class Mesh<Context::MPI> : public MeshBase
  {
    public:
      class Builder
      {
        public:
          Builder();

          Builder& initialize(const Context::MPI& context, Shard&& shard);

          Mesh finalize();

        private:
          Shard m_shard;
          Context::MPI m_context;
      };

      Mesh(const Context::MPI& context)
        : m_context(context)
      {}

      Shard& getShard();

      const Shard& getShard() const;

      Index getGlobalIndex(const std::pair<size_t, Index>& p, Index fragmentId);

      Mesh& scale(Real c) override;

      void flush() override;

      bool isSubMesh() const override;

      size_t getDimension() const override;

      size_t getSpaceDimension() const override;

      // const Connectivity<Context::MPI>& getConnectivity() const override
      // {
      //   return m_connectivity;
      // }

      const Context::MPI& getContext() const override;

      size_t getPolytopeCount(size_t d) const override;

      size_t getPolytopeCount(Polytope::Type g) const override;

      /**
       * @brief Returns a PolytopeIterator for the given global polytope index.
       * @param dimension The topological dimension of the polytopes to iterate over.
       * @param globalIdx The global index (consistent across processes) of the
       * polytope to start iteration from.
       *
       * This function maps a global polytope index into the local index space
       * using an exclusive scan on the local polytope count for the given
       * dimension. Each MPI process computes its own global offset (via an
       * exclusive scan) so that its local data is numbered within the global
       * index space.
       *
       * - If the provided global index is within the range that belongs to
       *   this process, the function constructs a PolytopeIterator that starts
       *   at the corresponding local index.
       * - Otherwise, the function returns an "empty" PolytopeIterator,
       *   indicating that the requested global polytope is not available on
       *   the calling process.
       *
       * @return A PolytopeIterator starting at the appropriate local index if
       * this process owns that global index; otherwise, an empty iterator.
       */
      PolytopeIterator getPolytope(size_t dimension, Index globalIdx) const override;

    private:
      Context::MPI m_context;
      Shard m_shard;
      Connectivity<Context::MPI> m_connectivity;
  };
}

#endif
#endif
