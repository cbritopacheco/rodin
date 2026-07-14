/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ASSEMBLY_H
#define RODIN_PETSC_ASSEMBLY_H

/**
 * @file
 * @brief Top level include for PETSc assembly backends.
 *
 * This header exposes PETSc-compatible assembly kernels used by variational
 * forms and problems. It aggregates sequential and MPI assembly, and includes
 * OpenMP assembly when enabled.
 *
 * Included components:
 * - @ref Rodin::PETSc::Assembly::Default "Default" selector
 * - Sequential assembly kernels
 * - MPI assembly kernels
 * - OpenMP assembly kernels (optional)
 */

#include "Assembly/Sequential.h"
#include "Assembly/MPI.h"

#ifdef RODIN_USE_OPENMP
#include "Assembly/OpenMP.h"
#endif

#endif
