/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_H
#define RODIN_PETSC_H

/**
 * @file
 * @brief Top level include for the Rodin::PETSc module.
 *
 * The PETSc module provides integration with the Portable, Extensible Toolkit
 * for Scientific Computation (PETSc). This enables the use of PETSc's powerful
 * parallel linear and nonlinear solvers with Rodin's finite element framework.
 *
 * Components include:
 * - **Math**: PETSc vector and matrix types
 * - **Solver**: PETSc-based linear and nonlinear solvers
 * - **Assembly**: PETSc-compatible matrix assembly
 * - **Variational**: Variational forms with PETSc backend
 * - **IO**: PETSc data input/output
 *
 * @note This module requires PETSc to be available and properly configured.
 *
 * @see Rodin::PETSc
 */

#include "PETSc/ForwardDecls.h"

#include "PETSc/IO.h"
#include "PETSc/Math.h"
#include "PETSc/Solver.h"
#include "PETSc/Assembly.h"
#include "PETSc/Variational.h"

#endif
