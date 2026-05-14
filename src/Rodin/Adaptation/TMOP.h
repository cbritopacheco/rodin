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
 * @brief Initial TMOP-style mesh quality infrastructure.
 *
 * These utilities provide the reusable pieces needed after level-set cutting:
 * quality metrics on target-normalized Jacobians, P2 geometry nodes on linear
 * triangles, objective evaluation, and a small damped optimizer skeleton. They
 * keep topology fixed; they do not remesh, coarsen, or alter interface
 * constraints.
 */

#include "TMOP/Metrics.h"
#include "TMOP/Functions.h"
#include "TMOP/Geometry.h"
#include "TMOP/Objective.h"
#include "TMOP/Optimizer.h"
#include "TMOP/Problem.h"

#endif
