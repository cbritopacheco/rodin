/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_MPICONNECTIVITY_H
#define RODIN_GEOMETRY_MPICONNECTIVITY_H

#include "Rodin/MPI/MPIContext.h"

#include "Rodin/Geometry/Connectivity.h"

namespace Rodin::Geometry
{
  using MPIConnectivity = Connectivity<Context::MPI>;

  template <>
  class Connectivity<Context::MPI> : public ConnectivityBase
  {
    public:
      using Context = Context::MPI;

      Connectivity(const Context& context)
      {}

  };
}

#endif

