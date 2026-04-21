/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the Fluid mechanics module.
 */
#ifndef RODIN_FLUID_FORWARDDECLS_H
#define RODIN_FLUID_FORWARDDECLS_H

namespace Rodin::Fluid
{
  class FlowPoint;

  namespace Tags
  {
    struct Density;
    struct DynamicViscosity;
    struct WallDistance;
    struct YieldStress;
    struct ConsistencyIndex;
    struct FlowIndex;
    struct Regularization;
    struct TurbulentKineticEnergy;
    struct DissipationRate;
    struct SpecificDissipationRate;
  }

  template <class Derived>
  class Input;

  class StrainRate;
  class ShearRate;
  class Vorticity;

  template <class LawDerived>
  class CauchyStress;

  template <class Derived>
  class RheologyLaw;

  class Newtonian;
  class PowerLaw;
  class CarreauYasuda;
  class Bingham;
  class HerschelBulkley;
}

#endif
