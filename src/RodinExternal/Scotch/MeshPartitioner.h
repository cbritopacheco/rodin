#ifndef RODIN_RODINEXTERNAL_SCOTCH_MESH_PARTITIONER_H
#define RODIN_RODINEXTERNAL_SCOTCH_MESH_PARTITIONER_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/MeshPartitioner.h"

#include <scotch.h>

namespace Rodin::External::Scotch
{
  class MeshPartitioner : public Geometry::MeshPartitioner
  {
    public:
      using MeshType = Geometry::Mesh<Context::Local>;

      MeshPartitioner(const MeshType& mesh);

      ~MeshPartitioner() override;

      const MeshType& getMesh() const override;

      void partition(size_t numPartitions)
      {
        partition(numPartitions, getMesh().getDimension());
      }

      void partition(size_t numPartitions, size_t d) override;

      size_t getPartition(Index index) const override;

      MeshPartitioner& setStrategy(const SCOTCH_Strat& strat);

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      std::vector<SCOTCH_Num> m_partition;

      SCOTCH_Graph m_graph;
      std::optional<SCOTCH_Strat> m_strat;
  };
}

#endif
