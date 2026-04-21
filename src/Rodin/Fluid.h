/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Fluid.h
 * @brief Aggregator header for the Rodin::Fluid module.
 *
 * The Fluid module provides quadrature-point closure infrastructure for
 * rheology and future turbulence/flow closures:
 *
 * - **Local**: FlowPoint state bundle, auxiliary typed tags, input injection
 * - **Fields**: Strain-rate, shear-rate, vorticity and Cauchy stress helpers
 * - **Constitutive**: CRTP rheology base and concrete generalized Newtonian laws
 *
 * This module intentionally focuses on local constitutive evaluation and does
 * not provide solver/problem wrappers or weak-form integrator catalogs.
 */
#ifndef RODIN_FLUID_H
#define RODIN_FLUID_H

// Forward declarations
#include "Fluid/ForwardDecls.h"

// Local
#include "Fluid/Local/Tags.h"
#include "Fluid/Local/Input.h"
#include "Fluid/Local/FlowPoint.h"

// Fields
#include "Fluid/Fields/StrainRate.h"
#include "Fluid/Fields/ShearRate.h"
#include "Fluid/Fields/Vorticity.h"
#include "Fluid/Fields/CauchyStress.h"

// Constitutive
#include "Fluid/Constitutive/RheologyLaw.h"
#include "Fluid/Constitutive/Newtonian.h"
#include "Fluid/Constitutive/PowerLaw.h"
#include "Fluid/Constitutive/CarreauYasuda.h"
#include "Fluid/Constitutive/Bingham.h"
#include "Fluid/Constitutive/HerschelBulkley.h"

#endif
