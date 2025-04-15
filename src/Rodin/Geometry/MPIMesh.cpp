#include "MPIMesh.h"

#ifdef RODIN_USE_MPI

namespace Rodin::Geometry
{
  MPIMesh::Builder& MPIMesh::Builder::initialize(const Context::MPI& context, Shard&& shard)
  {
    m_context = std::move(context);
    m_shard = std::move(shard);
    return *this;
  }

  // MPIMesh MPIMesh::Builder::finalize()
  // {
  //   MPIMesh mesh(m_context);
  //   mesh.m_shard = std::move(m_shard);
  //   return mesh;
  // }

  Index MPIMesh::getGlobalIndex(const std::pair<size_t, Index>& p, Index fragmentId)
  {
    return 0;
  }

  MPIMesh& MPIMesh::scale(Real c)
  {
    m_shard.scale(c);
    return *this;
  }

  void MPIMesh::flush()
  {
    m_shard.flush();
  }

  bool MPIMesh::isSubMesh() const
  {
    return false;
  }

  size_t MPIMesh::getDimension() const
  {
    return m_shard.getDimension();
  }

  size_t MPIMesh::getSpaceDimension() const
  {
    return m_shard.getSpaceDimension();
  }

  const Context::MPI& MPIMesh::getContext() const
  {
    return m_context;
  }

  size_t MPIMesh::getPolytopeCount(size_t d) const
  {
    const auto& comm = m_context.getCommunicator();
    size_t count = m_shard.getPolytopeCount(d);
    boost::mpi::all_reduce(comm, count, std::plus<size_t>());
    return count;
  }

  size_t MPIMesh::getPolytopeCount(Polytope::Type g) const
  {
    const auto& comm = m_context.getCommunicator();
    size_t count = m_shard.getPolytopeCount(g);
    boost::mpi::all_reduce(comm, count, std::plus<size_t>());
    return count;
  }

  Shard& MPIMesh::getShard()
  {
    return m_shard;
  }

  const Shard& MPIMesh::getShard() const
  {
    return m_shard;
  }

  PolytopeIterator MPIMesh::getPolytope(size_t dimension, Index globalIdx) const
  {
    // Get local count from the shard.
    size_t localCount = m_shard.getPolytopeCount(dimension);

    // Compute inclusive scan of localCount across processes.
    size_t inclusiveSum = 0;
    boost::mpi::communicator comm = m_context.getCommunicator();
    boost::mpi::scan(comm, localCount, inclusiveSum, std::plus<size_t>());

    // Compute exclusive scan by subtracting the local count.
    size_t globalOffset = inclusiveSum - localCount;

    // Check if the global index belongs to this rank.
    if (globalIdx < globalOffset || globalIdx >= globalOffset + localCount)
    {
      // Not owned by this rank: return an empty iterator.
      return PolytopeIterator(dimension, m_shard, BoundedIndexGenerator(0, 0));
    }
    else
    {
      // Map global index to local index.
      size_t localIdx = globalIdx - globalOffset;
      // Create an iterator that starts at this local index up to localCount.
      return PolytopeIterator(dimension, m_shard, BoundedIndexGenerator(localIdx, localCount));
    }
  }
}

#endif
