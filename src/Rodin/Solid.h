/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Solid.h
 * @brief Aggregator header for the Rodin::Solid module.
 *
 * The Solid module provides tools for linear and nonlinear (hyperelastic)
 * solid mechanics:
 *
 * - **Kinematics**: Deformation gradient, Cauchy-Green tensors, invariants
 * - **Constitutive laws**: Neo-Hookean, Saint-Venant-Kirchhoff, Mooney-Rivlin
 * - **Integrators**: Internal force, material tangent, body/traction forces
 * - **Linear**: Hooke's law, linear elasticity bridge
 * - **Post-processing**: Green-Lagrange strain, Piola-Kirchhoff stress, Cauchy stress
 *
 * ## Usage
 * @code
 * #include <Rodin/Solid.h>
 * using namespace Rodin::Solid;
 * @endcode
 */
#ifndef RODIN_SOLID_H
#define RODIN_SOLID_H

// Forward declarations
#include "Solid/ForwardDecls.h"

// Kinematics
#include "Solid/Kinematics/KinematicState.h"
#include "Solid/Kinematics/Invariants.h"

// Constitutive laws
#include "Solid/Constitutive/HyperElasticLaw.h"
#include "Solid/Constitutive/NeoHookean.h"
#include "Solid/Constitutive/SaintVenantKirchhoff.h"
#include "Solid/Constitutive/MooneyRivlin.h"

// Linear elasticity
#include "Solid/Linear/Hooke.h"
#include "Solid/Linear/LinearElasticity.h"

// Integrators
#include "Solid/Integrators/InternalForce.h"
#include "Solid/Integrators/MaterialTangent.h"
#include "Solid/Integrators/BodyForce.h"
#include "Solid/Integrators/TractionForce.h"

// Post-processing
#include "Solid/PostProcessing/GreenLagrangeStrain.h"
#include "Solid/PostProcessing/FirstPiolaKirchhoffStress.h"
#include "Solid/PostProcessing/CauchyStress.h"

#endif
