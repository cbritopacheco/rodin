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

          // Mesh finalize();

        private:
          Shard m_shard;
          Context::MPI m_context;
      };

      Mesh(const Context::MPI& context)
        : m_context(context)
      {}

      Shard& local();

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

    private:
      Context::MPI m_context;
      Shard m_shard;
      Connectivity<Context::MPI> m_connectivity;
  };
}

#endif
#endif
