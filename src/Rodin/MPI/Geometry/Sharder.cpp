/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
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
   * goes as a binary archive in chunks of this size. 1 GiB is the largest
   * power of two below INT_MAX bytes, which leaves room for the size header
   * that boost::mpi packs with each chunk.
   */
  constexpr size_t ChunkBytes = size_t(1) << 30;
  static_assert(ChunkBytes < size_t(std::numeric_limits<int>::max()));
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
    // One serialized shard at a time, released once sent, so that the root
    // holds at most one archive besides the shards left to send.
    for (size_t i = 0; i < shards.size(); i++)
    {
      if (i == static_cast<size_t>(root))
        continue;
      std::string buffer;
      {
        std::ostringstream os(std::ios::binary);
        {
          boost::archive::binary_oarchive oa(os);
          oa << shards[i];
        }
        shards[i] = Shard();
        buffer = std::move(os).str();
      }
      const std::uint64_t size = buffer.size();
      comm.send(static_cast<int>(i), tag, size);
      for (size_t offset = 0; offset < buffer.size(); offset += ChunkBytes)
      {
        const size_t n = std::min(ChunkBytes, buffer.size() - offset);
        comm.send(static_cast<int>(i), tag, buffer.data() + offset, static_cast<int>(n));
      }
    }
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
