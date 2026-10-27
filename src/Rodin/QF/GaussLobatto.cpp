/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include "GaussLobatto.h"

namespace Rodin::QF
{
  const QuadratureFormula& GaussLobatto::get(Geometry::Polytope::Type g, size_t count)
  {
    return QuadratureFormula::get<GaussLobatto>(g, count);
  }
}
