/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_MPI_CONVERGENCE_H
#define RODIN_TESTS_CONVERGENCE_MPI_CONVERGENCE_H

#include "Convergence.h"
#include "Rodin/Geometry/BalancedCompactPartitioner.h"
#include "Rodin/MPI/Geometry/Sharder.h"

namespace Rodin::Tests::Convergence
{
  /** Unit-box grids with explicit incidence requirements, partitioned on rank zero. */
  class DistributedUniformGrid
  {
    public:
      DistributedUniformGrid(
        const Context::MPI& context, Geometry::Polytope::Type geometry)
        : m_context(context),
          m_grid(geometry)
      {}

      Geometry::Mesh<Context::MPI> makeMesh(size_t pointsPerAxis) const
      {
        Geometry::Sharder<Context::MPI> sharder(m_context);
        const auto& comm = m_context.getCommunicator();
        if (comm.rank() == 0)
        {
          auto mesh = m_grid.makeMesh(pointsPerAxis);
          const size_t dim = mesh.getDimension();
          mesh.getConnectivity().compute(dim, dim);
          mesh.getConnectivity().compute(dim, 0);
          mesh.getConnectivity().compute(dim, dim - 1);
          mesh.getConnectivity().compute(dim - 1, dim);
          mesh.getConnectivity().compute(dim - 1, 0);
          Geometry::BalancedCompactPartitioner partitioner(mesh);
          partitioner.partition(static_cast<size_t>(comm.size()));
          sharder.shard(partitioner);
          sharder.scatter(0);
        }
        return sharder.gather(0);
      }

    private:
      Context::MPI m_context;
      UniformGrid m_grid;
  };
}

#endif
