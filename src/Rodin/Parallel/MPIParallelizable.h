#ifndef RODIN_MPIPROCESS_H
#define RODIN_MPIPROCESS_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_MPI

#include <optional>
#include <boost/mpi.hpp>

namespace Rodin::Parallel
{
   template <class ParallelClass>
   class MPIParallelizable
   {
      public:
         virtual ParallelClass setMPIComm(boost::mpi::communicator comm) = 0;
         virtual const boost::mpi::communicator& getMPIComm() const = 0;
   };
}

#else

namespace Rodin::Parallel
{
   struct MPIParallelizable
   {};
}

#endif
#endif
