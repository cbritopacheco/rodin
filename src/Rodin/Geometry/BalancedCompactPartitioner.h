#ifndef RODIN_GEOMETRY_BALANCEDCOMPACTPARTITIONER_H
#define RODIN_GEOMETRY_BALANCEDCOMPACTPARTITIONER_H

#include "MeshPartitioner.h"

namespace Rodin::Geometry
{
  class BalancedCompactPartitioner : public Partitioner
  {
    public:
      using MeshType = Mesh<Context::Local>;

      BalancedCompactPartitioner(const MeshType& mesh);

      virtual ~BalancedCompactPartitioner() = default;

      virtual const MeshType& getMesh() const override;

      virtual void partition(size_t maxPartitionSize, size_t d) override;

      void partition(size_t maxPartitionSize)
      {
        partition(maxPartitionSize, getMesh().getDimension());
      }

      virtual size_t getPartition(Index index) const override;

      size_t getCount() const override;

    private:
      size_t m_count;
      std::reference_wrapper<const MeshType> m_mesh;
      std::vector<size_t> m_partition;
  };
}

#endif
