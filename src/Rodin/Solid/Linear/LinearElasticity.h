/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LinearElasticity.h
 * @brief Linear elasticity utilities within the Solid module.
 *
 * Re-exports the existing Rodin linear elasticity integrator through the
 * Solid module namespace, and provides helper utilities.
 *
 * This header bridges the existing Rodin::Variational::LinearElasticityIntegral
 * into the Rodin::Solid namespace for consistency.
 */
#ifndef RODIN_SOLID_LINEAR_LINEARELASTICITY_H
#define RODIN_SOLID_LINEAR_LINEARELASTICITY_H

#include "Rodin/LinearElasticity.h"

#include "Hooke.h"

namespace Rodin::Solid
{
  /**
   * @brief Brings the existing linear elasticity integrator into the Solid namespace.
   *
   * This allows users to access the linear elasticity bilinear form integrator
   * through both `Rodin::Variational::LinearElasticityIntegral` and
   * `Rodin::Solid::LinearElasticityIntegral`.
   */
  using Variational::LinearElasticityIntegral;
  using Variational::LinearElasticityIntegrator;
}

#endif
