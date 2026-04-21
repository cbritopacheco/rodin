/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Tags.h
 * @brief Standard typed tags for Fluid::FlowPoint auxiliary data.
 */
#ifndef RODIN_FLUID_LOCAL_TAGS_H
#define RODIN_FLUID_LOCAL_TAGS_H

#include "Rodin/Types.h"

namespace Rodin::Fluid::Tags
{
  struct Density
  {
    using Type = Real;
  };

  struct DynamicViscosity
  {
    using Type = Real;
  };

  struct WallDistance
  {
    using Type = Real;
  };

  struct YieldStress
  {
    using Type = Real;
  };

  struct ConsistencyIndex
  {
    using Type = Real;
  };

  struct FlowIndex
  {
    using Type = Real;
  };

  struct Regularization
  {
    using Type = Real;
  };

  struct TurbulentKineticEnergy
  {
    using Type = Real;
  };

  struct DissipationRate
  {
    using Type = Real;
  };

  struct SpecificDissipationRate
  {
    using Type = Real;
  };
}

#endif
