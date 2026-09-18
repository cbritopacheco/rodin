/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "RotatedNitscheIntegrator.h"

namespace KelvinBall
{
  RotatedNitscheIntegrator::RotatedNitscheIntegrator(const Mesh& mesh,
    const FlatSet<Attribute>& masterCuts, Real physicalTolerance,
    Real referenceTolerance)
    : m_locator(mesh, masterCuts, physicalTolerance, referenceTolerance)
  {}

  const RotatedNitscheIntegrator::Locator& RotatedNitscheIntegrator::getLocator() const
  {
    return m_locator;
  }
}
