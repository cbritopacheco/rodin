/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
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


