/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Compatibility include for the Kokkos cuDSS integration.
 *
 * Link Rodin::Kokkos and initialize Kokkos before using Solver::CuDSS.
 * The implementation is provided by Rodin/Kokkos/Solver/CuDSS.h.
 */
#ifndef RODIN_SOLVER_CUDSS_H
#define RODIN_SOLVER_CUDSS_H

#include "Rodin/Configure.h"

#ifdef RODIN_USE_CUDSS
#include "Rodin/Kokkos/Solver/CuDSS.h"
#endif

#endif
