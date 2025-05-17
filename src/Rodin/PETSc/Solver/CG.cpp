#include "CG.h"

namespace Rodin::Solver
{
  CG<::Mat, ::Vec>::CG(ProblemType& pb)
    : Parent(pb)
  {
    setType(KSPCG);
  }

  CG<::Mat, ::Vec>::CG(const CG& other)
    : Parent(other)
  {}

  CG<::Mat, ::Vec>::CG(CG&& other)
    : Parent(std::move(other))
  {}

  CG<::Mat, ::Vec>::~CG() = default;

} // namespace Rodin::Solver

