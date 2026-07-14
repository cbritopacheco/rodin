/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_SOLVER_H
#define RODIN_PETSC_SOLVER_H

/**
 * @file
 * @brief Top level include for PETSc solver wrappers.
 *
 * Provides PETSc-backed linear and nonlinear solvers integrated with Rodin
 * variational problems. This includes generic KSP/SNES wrappers and convenient
 * aliases for common Krylov methods.
 */

#include "Solver/KSP.h"
#include "Solver/SNES.h"

#include "Solver/CG.h"
#include "Solver/GMRES.h"

#endif
