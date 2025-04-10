#ifndef RODIN_GEOMETRY_MESHSHARDER_H
#define RODIN_GEOMETRY_MESHSHARDER_H

#include "MeshPartitioner.h"

namespace Rodin::Geometry
{
  template <class TargetContext>
  class MeshSharder;

  template <>
  class MeshSharder<Context::MPI>
  {
    public:
      using ContextType = Context::MPI;

      MeshSharder() = default;

      virtual ~MeshSharder() = default;

      const Partitioner& getPartitioner() const;

      const Mesh<Context::Local>& getMesh() const;

      Mesh<ContextType> shard(size_t numShards);
  };
}

#endif
