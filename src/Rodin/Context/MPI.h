/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CONTEXT_MPI_H
#define RODIN_CONTEXT_MPI_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_MPI

#include <mpi.h>
#include <boost/mpi.hpp>

namespace Rodin::Context
{
  class MPI
  {
    boost::mpi::communicator comm;
  };
}

#endif
#endif

