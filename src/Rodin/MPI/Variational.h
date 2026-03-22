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
 * - Distributed @ref Rodin::Variational::P1 spaces.
 */

#include "Variational/P1.h"
#include "Variational/FiniteElementSpace.h"

#endif
