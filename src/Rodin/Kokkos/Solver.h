/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_SOLVER_H
#define RODIN_KOKKOS_SOLVER_H

/**
 * @file
 * @brief Kokkos Krylov solvers and optional vendor direct solvers.
 */

#include "Solver/CG.h"
#include "Solver/GMRES.h"

#ifdef RODIN_USE_ROCSOLVER
#include "Solver/RocSolverRF.h"
#endif

#endif
