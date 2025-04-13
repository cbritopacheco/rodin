#ifndef RODIN_GEOMETRY_GREEDYPARTITIONER_H
#define RODIN_GEOMETRY_GREEDYPARTITIONER_H

#include "MeshPartitioner.h"

namespace Rodin::Geometry
{
  class GreedyPartitioner : public Partitioner
  {
    public:
      using MeshType = Geometry::Mesh<Context::Local>;

      GreedyPartitioner(const MeshType& mesh);

      virtual ~GreedyPartitioner() = default;

      virtual const MeshType& getMesh() const override;

      void partition(size_t count)
      {
        partition(count, getMesh().getDimension());
      }

      virtual void partition(size_t numPartitions, size_t d) override;

      virtual size_t getPartition(Index index) const override;

      size_t getCount() const override;

    private:
      size_t m_count;
      std::reference_wrapper<const MeshType> m_mesh;
      std::vector<size_t> m_partition;
  };
}

#endif
