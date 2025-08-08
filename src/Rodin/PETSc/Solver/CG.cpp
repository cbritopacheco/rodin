#include "Rodin/PETSc/Math/LinearSystem.h"

#include "CG.h"

namespace Rodin::Solver
{
  CG<PETSc::Math::LinearSystem>::CG(ProblemBaseType& pb)
    : Parent(pb)
  {
    this->setType(KSPCG);
  }

  CG<PETSc::Math::LinearSystem>::CG(const CG& other)
    : Parent(other)
  {}

  CG<PETSc::Math::LinearSystem>::CG(CG&& other)
    : Parent(std::move(other))
  {}
}

