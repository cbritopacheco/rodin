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
