/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_ORDERING_H
#define RODIN_GEOMETRY_ORDERING_H

/**
 * @file
 * @brief Concrete local mesh vertex ordering policies.
 *
 * Morton and Hilbert expose bounding-box quantization bit depth. RCM uses
 * deterministic pseudo-peripheral roots and degree/index tie breaking.
 * Each policy is callable and returns a new-to-old vertex permutation.
 * Mesh::renumber() accepts each concrete policy through an explicit overload
 * and defaults to Morton when no policy is supplied.
 */

/**
 * @namespace Rodin::Geometry::Ordering
 * @brief Concrete algorithms for ordering local mesh vertices.
 *
 * Each algorithm computes a new-to-old vertex permutation and is accepted
 * by an explicit Mesh::renumber() overload. Morton and Hilbert provide
 * configurable quantization precision; RCM uses deterministic graph traversal.
 */

#include "Ordering/Morton.h"
#include "Ordering/Hilbert.h"
#include "Ordering/RCM.h"

#endif
