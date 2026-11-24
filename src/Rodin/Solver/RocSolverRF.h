/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Compatibility include for the Kokkos rocSOLVER integration.
 *
 * Link Rodin::Kokkos and initialize Kokkos before using Solver::RocSolverRF.
 * The implementation is provided by Rodin/Kokkos/Solver/RocSolverRF.h.
 */
#ifndef RODIN_SOLVER_ROCSOLVERRF_H
#define RODIN_SOLVER_ROCSOLVERRF_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_ROCSOLVER
#include "Rodin/Kokkos/Solver/RocSolverRF.h"
#endif

#endif
