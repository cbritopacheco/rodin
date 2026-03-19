#ifndef RODIN_MPI_GEOMETRY_H
#define RODIN_MPI_GEOMETRY_H

/**
 * @file
 * @brief Aggregated include for distributed geometry primitives.
 *
 * This header groups the MPI-enabled geometry components:
 * - @ref Rodin::Geometry::Mesh<Rodin::Context::MPI> for distributed mesh access
 * - @ref Rodin::Geometry::Sharder<Rodin::Context::MPI> for mesh partitioning
 *   and distribution across ranks.
 */

#include "Geometry/Mesh.h"
#include "Geometry/Sharder.h"

#endif
