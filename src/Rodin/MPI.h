/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_H
#define RODIN_MPI_H

/**
 * @file
 * @brief Top level include for the Rodin::MPI module.
 *
 * The MPI module provides distributed computing support using the Message
 * Passing Interface (MPI). This enables parallel finite element computations
 * across multiple processes and nodes.
 *
 * Components include:
 * - **Context**: MPI initialization and environment management
 * - **Assembly**: Distributed matrix and vector assembly
 * - **Variational**: Distributed variational formulations
 * - **Geometry**: Distributed mesh handling
 * - **IO**: Parallel input/output operations
 *
 * @note This module requires MPI to be available and properly configured.
 *
 * @see Rodin::MPI
 */

#include "MPI/Context.h"
#include "MPI/Assembly.h"
#include "MPI/Variational.h"
#include "MPI/Geometry.h"
#include "MPI/IO.h"

#endif
