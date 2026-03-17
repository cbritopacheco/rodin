#include <Rodin/Scotch/MeshPartitioner.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main()
{
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, {32, 32});

  mesh.getConnectivity().compute(2, 2);

  Scotch::Partitioner partitioner(mesh);

  partitioner.partition(16);

  for (auto it = mesh.getCell(); it; ++it)
  {
    const auto idx = it->getIndex();
    const auto part = partitioner.getPartition(idx);
    mesh.setAttribute({2, idx}, part + 1);
  }

  mesh.save("Partitioned.mfem", IO::FileFormat::MFEM);

  return 0;
}
