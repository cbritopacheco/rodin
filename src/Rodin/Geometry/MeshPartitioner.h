#ifndef RODIN_GEOMETRY_MESHPARTITIONER_H
#define RODIN_GEOMETRY_MESHPARTITIONER_H

#include "Mesh.h"

namespace Rodin::Geometry
{
  class Partitioner
  {
    public:
      Partitioner() = default;

      virtual ~Partitioner() = default;

      virtual const Mesh<Context::Local>& getMesh() const = 0;

      virtual void partition(size_t numPartitions, size_t d) = 0;

      virtual size_t getPartition(Index index) const = 0;

      size_t operator[](Index index) const
      {
        return getPartition(index);
      }

      virtual size_t getCount() const;
  };
}

#endif
