/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_IO_MEDIT_H
#define RODIN_MPI_IO_MEDIT_H

/**
 * @file
 * @brief Placeholder for MEDIT IO specializations on MPI meshes.
 *
 * This header reserves the extension point for
 * `MeshLoader<FileFormat::MEDIT, Context::MPI>` and
 * `MeshPrinter<FileFormat::MEDIT, Context::MPI>`.
 *
 * MPI workflows currently rely on HDF5 support for distributed shard
 * persistence; MEDIT MPI specializations can be added here without changing
 * public include structure.
 */

namespace Rodin::IO
{}

#endif
