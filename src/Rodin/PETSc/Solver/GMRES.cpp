#include "Rodin/PETSc/Math/LinearSystem.h"

#include "GMRES.h"

namespace Rodin::Solver
{
  GMRES<PETSc::Math::LinearSystem>::GMRES(ProblemBaseType& pb)
    : Parent(pb)
  {
    this->setType(KSPGMRES);
  }

  GMRES<PETSc::Math::LinearSystem>::GMRES(const GMRES& other)
    : Parent(other)
  {}

  GMRES<PETSc::Math::LinearSystem>::GMRES(GMRES&& other)
    : Parent(std::move(other))
  {}
}


