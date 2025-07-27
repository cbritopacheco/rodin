#include "CG.h"
#include "Rodin/PETSc/Math/LinearSystem.h"

namespace Rodin::Solver
{
  CG<PETSc::LinearSystem>::CG(ProblemBaseType& pb)
    : Parent(pb)
  {
    this->setType(KSPCG);
  }

  CG<PETSc::LinearSystem>::CG(const CG& other)
    : Parent(other)
  {}

  CG<PETSc::LinearSystem>::CG(CG&& other)
    : Parent(std::move(other))
  {}
}

