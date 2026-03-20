#ifndef RODIN_MPI_ASSEMBLY_H
#define RODIN_MPI_ASSEMBLY_H

/**
 * @file
 * @brief Aggregated include for distributed assembly backends.
 *
 * This header exposes MPI specializations of Rodin assembly policies and
 * executors, including:
 * - @ref Rodin::Assembly::Default<Rodin::Context::MPI>
 * - @ref Rodin::Assembly::MPI
 */

#include "Assembly/Default.h"
#include "Assembly/MPI.h"

#endif
