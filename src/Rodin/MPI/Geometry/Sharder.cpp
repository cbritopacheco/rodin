#include <cassert>
#include <boost/multi_array.hpp>

#include "Rodin/Geometry/Mesh.h"

#include "Rodin/Serialization/BitSet.h"
#include "Rodin/Serialization/FlatMap.h"
#include "Rodin/Serialization/FlatSet.h"
#include "Rodin/Serialization/Optional.h"

#include "Sharder.h"

namespace Rodin::Geometry
{
  Sharder<Context::MPI>::Sharder(const Context::MPI& context)
    : Parent(context)
  {}

  Mesh<Context::MPI> Sharder<Context::MPI>::distribute(Partitioner& p, int root)
  {
    this->shard(p);
    return this->scatter(root).gather(root);
  }

  Sharder<Context::MPI>& Sharder<Context::MPI>::scatter(int root)
  {
    const auto& ctx = this->getContext();
    const auto& comm = ctx.getCommunicator();
    const int tag = ctx.getEnvironment().collectives_tag();
    auto& shards = this->getShards();
    std::vector<boost::mpi::request> reqs(shards.size());
    for (size_t i = 0; i < shards.size(); i++)
    {
      assert(root >= 0);
      if (i == static_cast<size_t>(root))
        continue;
      reqs[i] = comm.isend(i, tag, shards[i]);
    }
    boost::mpi::wait_all(reqs.begin(), reqs.end());
    return *this;
  }

  MPIMesh Sharder<Context::MPI>::gather(int root)
  {
    const auto& ctx = this->getContext();
    const auto& comm = ctx.getCommunicator();
    const int tag = ctx.getEnvironment().collectives_tag();
    auto& shards = this->getShards();
    if (comm.rank() == root)
    {
      return Mesh<Context::MPI>::Builder(ctx).initialize(std::move(shards[root]))
                                             .finalize();
    }
    else
    {
      Shard s;
      comm.recv(root, tag, s);
      Mesh<Context::MPI>::Builder build(ctx);
      build.initialize(std::move(s));
      return build.finalize();
    }
  }
}
