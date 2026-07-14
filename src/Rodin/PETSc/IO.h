/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_IO_H
#define RODIN_PETSC_IO_H

/**
 * @file
 * @brief Top level include for PETSc-specific IO support.
 *
 * This header aggregates PETSc-backed grid function printers and loaders for
 * Rodin IO formats.
 *
 * Supported formats:
 * - @ref Rodin::IO::FileFormat::HDF5 "HDF5"
 * - @ref Rodin::IO::FileFormat::MFEM "MFEM"
 * - @ref Rodin::IO::FileFormat::MEDIT "MEDIT"
 */

#include "IO/GridFunctionPrinter.h"
#include "IO/HDF5.h"
#include "IO/MFEM.h"
#include "IO/MEDIT.h"

#endif
