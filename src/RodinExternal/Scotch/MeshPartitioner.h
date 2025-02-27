#ifndef RODIN_RODINEXTERNAL_SCOTCH_MESH_PARTITIONER_H
#define RODIN_RODINEXTERNAL_SCOTCH_MESH_PARTITIONER_H

#include "Rodin/Geometry/MeshPartitioner.h"

namespace Rodin::External::Scotch
{
  class MeshPartitioner : public Geometry::MeshPartitioner
  {
    public:
      using MeshType = Geometry::Mesh<Context::Local>;

      MeshPartitioner(const MeshType& mesh);

      const MeshType& getMesh() const override;

      void partition(size_t numPartitions, size_t d) override;

      size_t getPartition(Index index) const override;

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      std::vector<Index> m_partition;
  };
}

#endif
