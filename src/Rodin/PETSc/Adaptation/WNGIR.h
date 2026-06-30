/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_ADAPTATION_WNGIR_H
#define RODIN_PETSC_ADAPTATION_WNGIR_H

/**
 * @file
 * @brief PETSc-backed WNGIR public include.
 *
 * Provides the PETSc specialization of @ref Rodin::Adaptation::WNGIRBackend, so
 * that the backend-independent @ref Rodin::Adaptation::WNGIR solver, when
 * constructed from a PETSc displacement GridFunction, runs its step problems on
 * PETSc objects with no per-backend solver code.
 */

#include "WNGIRBackend.h"

#endif
