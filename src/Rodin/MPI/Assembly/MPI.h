#ifndef RODIN_MPI_ASSEMBLY_MPI_H
#define RODIN_MPI_ASSEMBLY_MPI_H

#include "Rodin/Variational/Integrator.h"

#include "Rodin/MPI/Geometry/Mesh.h"

namespace Rodin::Assembly
{
  class MPIIteration
  {
    public:
      using MeshType = Geometry::Mesh<Context::MPI>;

      MPIIteration(const MeshType& mesh, Variational::Integrator::Region);

      Geometry::PolytopeIterator getIterator() const;

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      Variational::Integrator::Region m_region;
  };
}

namespace Rodin::Assembly
{
  template <class LinearAlgebraType, class Operand>
  class MPI;
}

#endif
