/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
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

