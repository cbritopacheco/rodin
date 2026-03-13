/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_H
#define RODIN_IO_H

/**
 * @file
 * @brief Top level include for the Rodin::IO module.
 *
 * The IO module provides input/output functionality for meshes, grid functions,
 * and other finite element data. It supports various file formats commonly used
 * in computational mechanics and scientific visualization.
 *
 * Supported formats include:
 * - MEDIT (.mesh)
 * - MFEM (.mfem)
 *
 * @see Rodin::IO
 */

#include "IO/ForwardDecls.h"
#include "IO/Loader.h"
#include "IO/Printer.h"
#include "IO/HDF5.h"
#include "IO/XDMF.h"

#endif
