/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_IO_HDF5_H
#define RODIN_MPI_IO_HDF5_H

/**
 * @file HDF5.h
 * @brief HDF5 mesh persistence specializations for distributed MPI meshes.
 *
 * Provides MeshPrinter and MeshLoader specializations for
 * `Mesh<Context::MPI>` using the HDF5 file format.
 *
 * Each MPI rank persists its own local shard independently: the printer
 * delegates to `MeshPrinter<FileFormat::HDF5, Context::Local>` operating on
 * the mesh shard, and the loader delegates to
 * `MeshLoader<FileFormat::HDF5, Context::Local>` operating on the shard.
 *
 * @see MeshPrinter
 * @see MeshLoader
 * @see Geometry::Mesh<Context::MPI>
 */

#include "Rodin/IO/HDF5.h"
#include "Rodin/IO/MeshPrinter.h"
#include "Rodin/IO/MeshLoader.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/MPI/Geometry/Mesh.h"

namespace Rodin::IO
{
  /**
   * @brief Loads a distributed MPI mesh from an HDF5 file.
   *
   * Each MPI rank independently loads its local shard from the given HDF5
   * file path using the canonical local mesh HDF5 layout. The file is
   * expected to contain the standard datasets produced by
   * `MeshPrinter<FileFormat::HDF5, Context::Local>`.
   *
   * @note Each rank must be given a rank-specific file path (e.g. via
   *       the callable filename overload on Mesh<Context::MPI>::load).
   *
   * @see MeshLoader<FileFormat::HDF5, Context::Local>
   */
  template <>
  class MeshLoader<FileFormat::HDF5, Context::MPI>
    : public MeshLoaderBase<Context::MPI>
  {
    public:
      using ContextType = Context::MPI;
      using ObjectType = Geometry::Mesh<ContextType>;
      using Parent = MeshLoaderBase<ContextType>;

      /**
       * @brief Constructs a loader for the given distributed mesh.
       * @param[in] mesh  Distributed mesh whose local shard will be loaded.
       */
      MeshLoader(ObjectType& mesh)
        : Parent(mesh)
      {}

      /**
       * @brief Stream-based loading is not supported for HDF5 format.
       * @throws Alert::MemberFunctionException always.
       */
      void load(std::istream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 MPI mesh loading is file-path based. "
          << "Please use the "
          << Alert::Identifier::Function("load(const boost::filesystem::path&)")
          << " overload."
          << Alert::Raise;
      }

      /**
       * @brief Loads the local mesh shard from the given HDF5 file.
       *
       * Delegates to `MeshLoader<FileFormat::HDF5, Context::Local>`
       * operating on the shard of the distributed mesh.
       *
       * @param[in] filename  Path to the HDF5 file for this rank's shard.
       */
      void load(const boost::filesystem::path& filename) override
      {
        auto& mesh = this->getObject();
        auto& shard = mesh.getShard();
        MeshLoader<FileFormat::HDF5, Context::Local> localLoader(shard);
        localLoader.load(filename);
      }
  };

  /**
   * @brief Prints a distributed MPI mesh to an HDF5 file.
   *
   * Each MPI rank independently saves its local shard to the given HDF5
   * file path using the canonical local mesh HDF5 layout.
   *
   * @note Each rank must be given a rank-specific file path (e.g. via
   *       the callable filename overload on Mesh<Context::MPI>::save).
   *
   * @see MeshPrinter<FileFormat::HDF5, Context::Local>
   */
  template <>
  class MeshPrinter<FileFormat::HDF5, Context::MPI>
    : public MeshPrinterBase<Context::MPI>
  {
    public:
      using ContextType = Context::MPI;
      using ObjectType = Geometry::Mesh<ContextType>;
      using Parent = MeshPrinterBase<ContextType>;

      /**
       * @brief Constructs a printer for the given distributed mesh.
       * @param[in] mesh  Distributed mesh whose local shard will be saved.
       */
      MeshPrinter(const ObjectType& mesh)
        : Parent(mesh)
      {}

      /**
       * @brief Stream-based printing is not supported for HDF5 format.
       * @throws Alert::MemberFunctionException always.
       */
      void print(std::ostream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 MPI mesh printing is file-path based. "
          << "Please use the "
          << Alert::Identifier::Function("print(const boost::filesystem::path&)")
          << " overload."
          << Alert::Raise;
      }

      /**
       * @brief Saves the local mesh shard to the given HDF5 file.
       *
       * Delegates to `MeshPrinter<FileFormat::HDF5, Context::Local>`
       * operating on the shard of the distributed mesh.
       *
       * @param[in] filename  Path to the HDF5 file for this rank's shard.
       */
      void print(const boost::filesystem::path& filename) override
      {
        const auto& mesh = this->getObject();
        const auto& shard = mesh.getShard();
        MeshPrinter<FileFormat::HDF5, Context::Local> localPrinter(shard);
        localPrinter.print(filename);
      }
  };
}

#endif
