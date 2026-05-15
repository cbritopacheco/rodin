/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TMOP_H
#define RODIN_ADAPTATION_TMOP_H

/**
 * @file
 * @brief Rodin-native TMOP mesh adaptation components.
 *
 * The production direction is fixed-topology TMOP assembled as ordinary Rodin
 * residual and tangent terms into LinearProblem objects, then solved by
 * NewtonSolver and the existing local/PETSc/MPI/OpenMP backends.
 *
 * TMOP exposes metrics, target Jacobians, and Rodin-native residual/tangent
 * terms. The residual/tangent terms are the production path: they assemble into
 * LinearProblem objects and use the existing solver backends. The detached
 * high-order geometry/objective/optimizer helpers remain diagnostic bridges and
 * are not the production nonlinear solve architecture.
 */

#include "TMOP/Metrics.h"
#include "TMOP/Functions.h"
#include "TMOP/Terms.h"

// Fixed-topology geometry and adaptation helpers.
#include "TMOP/Geometry.h"
#include "TMOP/Objective.h"
#include "TMOP/Optimizer.h"
#include "TMOP/Problem.h"

#endif
