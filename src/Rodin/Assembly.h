/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ASSEMBLY_H
#define RODIN_ASSEMBLY_H

/**
 * @file
 * @brief Top level include for the Rodin::Assembly module.
 *
 * The Assembly module provides various strategies for assembling finite element
 * matrices and vectors. Assembly is the process of computing the global system
 * @f$ Ax = b @f$ from local element contributions.
 *
 * Available assembly strategies:
 * - **Sequential**: Single-threaded assembly (always available)
 * - **OpenMP**: Multi-threaded parallel assembly (requires RODIN_USE_OPENMP)
 *
 * @see Rodin::Assembly
 */

#include "Assembly/Generic.h"
#include "Assembly/Sequential.h"
#include "Assembly/AssemblyBase.h"

#ifdef RODIN_USE_OPENMP
#include "Assembly/OpenMP.h"
#endif

#include "Assembly/Default.h"

#endif