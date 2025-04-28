#include "PETSc.h"

namespace Rodin::PETSc
{
  void copyVector(const Math::Vector<PetscScalar>& source, ::Vec& destination)
  {
    throw std::runtime_error("Not implemented");
  }
}
