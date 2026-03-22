/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_MPI_H
#define RODIN_MPI_MPI_H

/**
 * @file
 * @brief Definition of the MPI execution context wrapper.
 */

#include <mpi.h>
#include <boost/mpi.hpp>

#include "Rodin/Context/Base.h"

namespace Rodin::Context
{
  /**
   * @brief MPI context for parallel computing with Message Passing Interface.
   *
   * This class provides an MPI context that wraps the Boost.MPI environment
   * and communicator objects. It serves as the interface for parallel 
   * computations across multiple processes in distributed memory systems.
   *
   * The context is passed to distributed mesh, assembly, and finite-element
   * components so they can access rank-local and collective communication
   * services in a uniform way.
   */
  class MPI : public Context::Base
  {
    public:
      /**
       * @brief Constructs an MPI context.
       * @param env MPI environment managing MPI initialization/finalization
       * @param world MPI communicator for process communication
       */
      MPI(const boost::mpi::environment& env, const boost::mpi::communicator& world)
          : m_env(env),
            m_world(world)
      {}

      /**
       * @brief Gets the MPI communicator.
       * @return Const reference to the MPI communicator
       */
      const boost::mpi::communicator& getCommunicator() const
      {
        return m_world.get();
      }

      /**
       * @brief Gets the MPI environment.
       * @return Const reference to the MPI environment
       */
      const boost::mpi::environment& getEnvironment() const
      {
        return m_env.get();
      }

    private:
      /// Reference to the MPI environment
      std::reference_wrapper<const boost::mpi::environment> m_env;

      /// Reference to the MPI communicator
      std::reference_wrapper<const boost::mpi::communicator> m_world;
  };
}

#endif
