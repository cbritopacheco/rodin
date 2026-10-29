/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_SOLID_PARAMETERS_H
#define RODIN_ATRIUM_SOLID_PARAMETERS_H

#include <Eigen/Core>
#include <Rodin/Geometry.h>
#include <Rodin/Solid/Constitutive/DemirayHolzapfelGasserOgden.h>

namespace Rodin::Examples::Heart::Atrium
{
  using Vec = Eigen::Matrix<Real, 3, 1>;

  /// @brief Passive, viscous, active, and prescribed-field parameters of the wall.
  struct SolidParameters
  {
    Geometry::Attribute baseLabel = 6;
    Vec appendageApex = Vec(46.39638367, 187.44307667, -205.03317767);
    Real apexRadius = 3;
    bool fibersOnly = false;

    /// @brief Reuses active constitutive blocks within a synchronized Newton state.
    bool cacheConstitutive = false;
    Real activation = 20; // Positive electrical peak, in 1/s.
    Real relaxation = 10; // Magnitude of the resting electrical drive, in 1/s.
    Real speed = 200;
    Real width = 6;

    /// @brief Literature-based atrial passive coefficients; see the solid README.
    Rodin::Solid::DemirayHolzapfelGasserOgden::Parameters passive{.matrixScale = 2.92,
      .matrixExponent = 5.6,
      .fiberScale = 11.84,
      .fiberExponent = 17.95,
      .volumePenalty = 500};
    Real viscosity = 0.2; // Tissue viscosity, in kPa s.
    Real activeStiffness = 20; // Serial elastic stiffness, in kPa.
    Real activeDamping = 3; // Contractile dashpot, in kPa s.
    Real crossBridgeStiffness = 100; // Cross-bridge stiffness source, in kPa.
    Real contractility = 20; // Cross-bridge stress source, in kPa.
    Real destructionRate = 0.4; // Dimensionless strain-rate destruction factor.
    Real density = 1060; // Reference mass density in kg/m^3.
  };

}

#endif
