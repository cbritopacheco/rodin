/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <sstream>
#include <string>
#include <boost/multi_array.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "Rodin/Geometry/Mesh.h"

#include "Rodin/Serialization/BitSet.h"
#include "Rodin/Serialization/FlatMap.h"
#include "Rodin/Serialization/FlatSet.h"
#include "Rodin/Serialization/Optional.h"

#include "Sharder.h"

// Instantiate exports after the packed MPI and binary archive declarations.
// This single translation unit owns the built-in polymorphic geometry
// registrations.
#include "Rodin/Serialization/Export.h"

namespace
{
  /**
   * @brief Largest message of a shard transfer.
   *
   * MPI counts are int, and boost::mpi packs with int positions, so a shard
   * goes as a binary archive in chunks of this size.
   */
  constexpr size_t ChunkBytes = size_t(1) << 30;
}

namespace Rodin::Geometry
{
  /// @cond
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
    assert(root >= 0);
    std::vector<std::string> buffers(shards.size());
    std::vector<std::uint64_t> sizes(shards.size(), 0);
    std::vector<boost::mpi::request> reqs;
    for (size_t i = 0; i < shards.size(); i++)
    {
      if (i == static_cast<size_t>(root))
        continue;
      {
        std::ostringstream os(std::ios::binary);
        boost::archive::binary_oarchive oa(os);
        oa << shards[i];
        buffers[i] = std::move(os).str();
      }
      sizes[i] = buffers[i].size();
      reqs.push_back(comm.isend(static_cast<int>(i), tag, sizes[i]));
      for (size_t offset = 0; offset < buffers[i].size(); offset += ChunkBytes)
      {
        const size_t n = std::min(ChunkBytes, buffers[i].size() - offset);
        reqs.push_back(comm.isend(static_cast<int>(i), tag,
          buffers[i].data() + offset, static_cast<int>(n)));
      }
    }
    boost::mpi::wait_all(reqs.begin(), reqs.end());
    return *this;
  }

  Mesh<Context::MPI> Sharder<Context::MPI>::gather(int root)
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
      std::uint64_t size = 0;
      comm.recv(root, tag, size);
      std::string buffer(size, '\0');
      for (size_t offset = 0; offset < buffer.size(); offset += ChunkBytes)
      {
        const size_t n = std::min<size_t>(ChunkBytes, buffer.size() - offset);
        comm.recv(root, tag, buffer.data() + offset, static_cast<int>(n));
      }
      Shard s;
      {
        std::istringstream is(std::move(buffer), std::ios::binary);
        boost::archive::binary_iarchive ia(is);
        ia >> s;
      }
      Mesh<Context::MPI>::Builder build(ctx);
      build.initialize(std::move(s));
      return build.finalize();
    }
  }
  /// @endcond
}
