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
 * TMOP exposes metrics, target Jacobians, and Rodin-native residual/tangent
 * terms. The residual/tangent terms are the production path: they assemble into
 * LinearProblem objects and use the existing solver backends. The detached
 * high-order geometry/objective/optimizer helpers remain diagnostic bridges and
 * are not the production nonlinear solve architecture.
 */

#include "TargetMatrixOptimization/Metrics.h"
#include "TargetMatrixOptimization/Functions.h"
#include "TargetMatrixOptimization/Terms.h"

// Fixed-topology geometry and adaptation helpers.
#include "TargetMatrixOptimization/Geometry.h"
#include "TargetMatrixOptimization/Objective.h"
#include "TargetMatrixOptimization/Optimizer.h"
#include "TargetMatrixOptimization/Problem.h"

// Interface-vertex tangential redistribution (Laplacian-with-projection),
// applied between the topological stage and the geometric (TMOP) stage.
#include "TargetMatrixOptimization/InterfaceTangentialRelaxation.h"

// Isoparametric (curved) H1<2> TMOP: shared geometry + Newton solver.
#include "TargetMatrixOptimization/IsoparametricGeometry.h"
#include "TargetMatrixOptimization/IsoparametricTMOPSolver.h"

#endif
