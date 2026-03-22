/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MMG.h
 * @brief Public umbrella header for the @ref Rodin::MMG module.
 *
 * Include this header to access Rodin's C++ wrappers around the
 * [Mmg Platform](https://www.mmgtools.org/) for:
 * - mesh optimization,
 * - metric-based mesh adaptation,
 * - implicit-domain and level-set discretization workflows.
 *
 * The module extends @ref Rodin::Geometry::Mesh "Geometry::Mesh" with MMG-specific
 * boundary tags and provides high-level operators that internally dispatch to
 * MMG2D, MMG3D, or MMGS depending on mesh dimension and manifold type.
 */
#ifndef RODIN_EXTERNAL_MMG_H
#define RODIN_EXTERNAL_MMG_H

#include "MMG/MMG5.h"

#include "MMG/Mesh.h"
#include "MMG/Adapt.h"
#include "MMG/GridFunction.h"
#include "MMG/MeshOptimizer.h"
#include "MMG/LevelSetDiscretizer.h"

#endif
