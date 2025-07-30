/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_MPI_H
#define RODIN_MPI_MPI_H

#include <mpi.h>
#include <boost/mpi.hpp>

#include "Rodin/Context/Base.h"

namespace Rodin::Context
{
  class MPI : public Context::Base
  {
    public:
      MPI(const boost::mpi::environment& env, const boost::mpi::communicator& world)
          : m_env(env),
            m_world(world)
      {}

      const boost::mpi::communicator& getCommunicator() const
      {
        return m_world.get();
      }

      const boost::mpi::environment& getEnvironment() const
      {
        return m_env.get();
      }

    private:
      std::reference_wrapper<const boost::mpi::environment> m_env;
      std::reference_wrapper<const boost::mpi::communicator> m_world;
  };
}

#endif

