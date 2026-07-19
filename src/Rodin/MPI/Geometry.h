/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_GEOMETRY_H
#define RODIN_MPI_GEOMETRY_H

/**
 * @file
 * @brief Aggregated include for distributed geometry primitives.
 *
 * This header groups the MPI-enabled geometry components:
 * - @ref Rodin::Geometry::Mesh<Rodin::Context::MPI> for distributed mesh access
 * - @ref Rodin::Geometry::SubMesh<Rodin::Context::MPI> for distributed submeshes
 * - @ref Rodin::Geometry::Sharder<Rodin::Context::MPI> for mesh partitioning
 *   and distribution across ranks.
 */

#include "Geometry/Mesh.h"
#include "Geometry/SubMesh.h"
#include "Geometry/Sharder.h"

#endif
