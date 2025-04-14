#ifndef RODIN_RODINEXTERNAL_SCOTCH_MESH_PARTITIONER_H
#define RODIN_RODINEXTERNAL_SCOTCH_MESH_PARTITIONER_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/MeshPartitioner.h"

#include <scotch.h>

namespace Rodin::External::Scotch
{
  class Partitioner : public Geometry::Partitioner
  {
    public:
      using MeshType = Geometry::Mesh<Context::Local>;

      Partitioner(const MeshType& mesh);

      ~Partitioner() override;

      const MeshType& getMesh() const override;

      void partition(size_t count)
      {
        partition(count, getMesh().getDimension());
      }

      void partition(size_t numPartitions, size_t d) override;

      size_t getPartition(Index index) const override;

      Partitioner& setStrategy(const SCOTCH_Strat& strat);

      size_t getCount() const override
      {
        return m_numPartitions;
      }

    private:
      size_t m_numPartitions;
      std::reference_wrapper<const MeshType> m_mesh;
      std::vector<SCOTCH_Num> m_partition;

      SCOTCH_Graph m_graph;
      std::optional<SCOTCH_Strat> m_strat;
  };
}

#endif
