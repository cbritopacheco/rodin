/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_KOKKOS_H
#define RODIN_PETSC_KOKKOS_H

/**
 * @file
 * @brief PETSc algebra assembled using local or distributed Kokkos meshes.
 */

#include "Rodin/PETSc.h"
#include "Rodin/Kokkos.h"
#ifdef RODIN_USE_MPI
#include "Rodin/Kokkos/MPI.h"
#endif
#include "Rodin/PETSc/Kokkos/Assembly/Kokkos.h"

#endif
