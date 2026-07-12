/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_VARIATIONAL_H
#define RODIN_MPI_VARIATIONAL_H

/**
 * @file
 * @brief Aggregated include for MPI variational finite element components.
 *
 * This header provides distributed finite element-space abstractions for MPI
 * meshes, including:
 * - @ref Rodin::Variational::FiniteElementSpace specializations on
 *   @ref Rodin::Geometry::Mesh<Rodin::Context::MPI>
 * - Distributed @ref Rodin::Variational::P0 spaces.
 * - Distributed @ref Rodin::Variational::P0g spaces.
 * - Distributed @ref Rodin::Variational::P1 spaces.
 * - Distributed @ref Rodin::Variational::H1 spaces.
 */

#include "Variational/P0.h"
#include "Variational/P0g.h"
#include "Variational/P1.h"
#include "Variational/H1.h"
#include "Variational/FiniteElementSpace.h"

#endif
