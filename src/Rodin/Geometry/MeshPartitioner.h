#ifndef RODIN_GEOMETRY_MESHPARTITIONER_H
#define RODIN_GEOMETRY_MESHPARTITIONER_H

#include "Mesh.h"

namespace Rodin::Geometry
{
    class MeshPartitioner
    {
      public:
        MeshPartitioner() = default;

        virtual ~MeshPartitioner() = default;

        virtual const MeshBase& getMesh() const = 0;

        virtual void partition(size_t numPartitions, size_t d) = 0;

        virtual size_t getPartition(Index index) const = 0;
    };
}

#endif
