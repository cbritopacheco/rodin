#include "CG.h"

namespace Rodin::Solver
{
  CG<PETSc::Matrix, PETSc::Vector>::CG(ProblemType& pb)
    : Parent(pb)
  {
    setType(KSPCG);
  }

  CG<PETSc::Matrix, PETSc::Vector>::CG(const CG& other)
    : Parent(other)
  {}

  CG<PETSc::Matrix, PETSc::Vector>::CG(CG&& other)
    : Parent(std::move(other))
  {}
}

