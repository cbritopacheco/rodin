/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_H

/**
 * @file
 * @brief Rodin-native TMOP mesh adaptation components.
 *
 * The production direction is fixed-topology TMOP assembled as ordinary Rodin
 * residual and tangent terms into LinearProblem objects, then solved by
 * NewtonSolver and the existing local/PETSc/MPI/OpenMP backends.
 *
 * TMOP exposes metrics, target Jacobians, deformation helpers, and Rodin-native
 * residual/tangent terms. The residual/tangent terms are the production path:
 * they assemble into LinearProblem objects and use the existing solver
 * backends.
 */

#include "TargetMatrixOptimization/Metrics.h"
#include "TargetMatrixOptimization/Terms.h"
#include "TargetMatrixOptimization/IsoparametricGeometry.h"
#include "TargetMatrixOptimization/IsoparametricTMOPSolver.h"

#endif
