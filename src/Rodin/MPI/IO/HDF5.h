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
 * delegates to `MeshPrinter<FileFormat::HDF5, Context::Local>` to write the
 * base mesh topology and then appends shard-specific metadata (ownership
 * flags, polytope index maps, ghost-to-owner maps, and halo CSR) under the
 * `/Shard/...` HDF5 group.
 *
 * The loader delegates to `MeshLoader<FileFormat::HDF5, Context::Local>` to
 * restore the base mesh and then reconstructs the MPI shard metadata needed
 * for distributed assembly from the `/Shard/...` HDF5 datasets written by the
 * MPI printer.
 *
 * @see MeshPrinter
 * @see MeshLoader
 * @see Geometry::Mesh<Context::MPI>
 * @see Geometry::Shard
 */

#include "Rodin/IO/HDF5.h"
#include "Rodin/IO/MeshPrinter.h"
#include "Rodin/IO/MeshLoader.h"
#include "Rodin/IO/GridFunctionPrinter.h"
#include "Rodin/IO/GridFunctionLoader.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/MPI/Geometry/Mesh.h"

namespace Rodin::IO
{
  /**
   * @brief Loads a distributed MPI mesh from an HDF5 file.
   *
   * Each MPI rank independently loads its local shard from the given HDF5 file
   * path using the canonical local mesh HDF5 layout and the vertex/cell shard
   * metadata stored under `/Shard/...`.  The file is expected to have been
   * produced by `MeshPrinter<FileFormat::HDF5, Context::MPI>`.
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
      /**
       * @brief Distributed context type handled by this loader.
       */
      using ContextType = Context::MPI;
      /**
       * @brief Distributed mesh object type loaded by this specialization.
       */
      using ObjectType = Geometry::Mesh<ContextType>;
      /**
       * @brief Base loader interface specialization.
       */
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
       * @param[in] filename  Path to the HDF5 file for this rank's shard.
       */
      void load(const boost::filesystem::path& filename) override
      {
        auto& mesh = this->getObject();
        mesh.getShard() = loadShard(filename);
      }

    private:
      struct ShardDimensionMetadata
      {
        std::vector<HDF5::U8> state;
        std::vector<HDF5::U64> left;
        std::vector<HDF5::U64> ownerKeys;
        std::vector<HDF5::U64> ownerValues;
        std::vector<HDF5::U64> haloKeys;
        std::vector<HDF5::U64> haloOffsets;
        std::vector<HDF5::U64> haloIndices;
      };

      static Geometry::Shard loadShard(const boost::filesystem::path& filename)
      {
        Geometry::Mesh<Context::Local> baseMesh;
        baseMesh.load(filename, FileFormat::HDF5);

        const auto file =
          HDF5::File(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!file)
        {
          Alert::Exception()
            << "Failed to open HDF5 file for shard metadata: " << filename
            << Alert::Raise;
        }

        const size_t Dmax = baseMesh.getDimension();
        const auto vertexMetadata =
          readShardDimensionMetadata(file.get(), baseMesh, 0);
        const auto cellMetadata =
          readShardDimensionMetadata(file.get(), baseMesh, Dmax);

        Geometry::Shard::Builder builder;
        builder.initialize(baseMesh);

        includeDimension(builder, 0, vertexMetadata);
        if (Dmax > 0)
          includeDimension(builder, Dmax, cellMetadata);

        applyOwnership(builder, 0, vertexMetadata);
        if (Dmax > 0)
          applyOwnership(builder, Dmax, cellMetadata);

        Geometry::Shard shard = builder.finalize();
        restorePolytopeMap(shard, 0, vertexMetadata.left);
        if (Dmax > 0)
          restorePolytopeMap(shard, Dmax, cellMetadata.left);

        return shard;
      }

      static ShardDimensionMetadata readShardDimensionMetadata(
          hid_t file,
          const Geometry::Mesh<Context::Local>& baseMesh,
          size_t d)
      {
        ShardDimensionMetadata metadata;
        metadata.state =
          HDF5::readVectorDataset<HDF5::U8>(file, HDF5::shardStatePath(d));
        metadata.left =
          HDF5::readVectorDataset<HDF5::U64>(file, HDF5::shardPolytopeMapLeftPath(d));

        const size_t count = baseMesh.getPolytopeCount(d);
        if (metadata.state.size() != count || metadata.left.size() != count)
        {
          Alert::Exception()
            << "Invalid shard metadata size for dimension " << d
            << ": expected " << count << " entries, got "
            << metadata.state.size() << " state entries and "
            << metadata.left.size() << " map entries."
            << Alert::Raise;
        }

        const std::string ownerPath = HDF5::shardOwnerGroupPath(d);
        metadata.ownerKeys =
          HDF5::readVectorDataset<HDF5::U64>(file, ownerPath + "/Keys");
        metadata.ownerValues =
          HDF5::readVectorDataset<HDF5::U64>(file, ownerPath + "/Values");
        if (metadata.ownerKeys.size() != metadata.ownerValues.size())
        {
          Alert::Exception()
            << "Invalid shard owner metadata for dimension " << d
            << ": key/value sizes differ."
            << Alert::Raise;
        }

        const std::string haloPath = HDF5::shardHaloGroupPath(d);
        metadata.haloKeys =
          HDF5::readVectorDataset<HDF5::U64>(file, haloPath + "/Keys");
        metadata.haloOffsets =
          HDF5::readVectorDataset<HDF5::U64>(file, haloPath + "/Offsets");
        metadata.haloIndices =
          HDF5::readVectorDataset<HDF5::U64>(file, haloPath + "/Indices");
        if (metadata.haloOffsets.size() != metadata.haloKeys.size() + 1)
        {
          Alert::Exception()
            << "Invalid shard halo metadata for dimension " << d
            << ": offsets size must equal keys size plus one."
            << Alert::Raise;
        }

        return metadata;
      }

      static Geometry::Shard::State stateOf(HDF5::U8 flag)
      {
        if (flag == 1)
          return Geometry::Shard::State::Owned;
        if (flag == 2)
          return Geometry::Shard::State::Ghost;
        return Geometry::Shard::State::Shared;
      }

      static void includeDimension(
          Geometry::Shard::Builder& builder,
          size_t d,
          const ShardDimensionMetadata& metadata)
      {
        for (Index i = 0; i < static_cast<Index>(metadata.state.size()); ++i)
          builder.include({ d, i }, stateOf(metadata.state[i]));
      }

      static void applyOwnership(
          Geometry::Shard::Builder& builder,
          size_t d,
          const ShardDimensionMetadata& metadata)
      {
        for (size_t i = 0; i < metadata.ownerKeys.size(); ++i)
        {
          builder.setOwner(
              d,
              static_cast<Index>(metadata.ownerKeys[i]),
              static_cast<Index>(metadata.ownerValues[i]));
        }

        for (size_t i = 0; i < metadata.haloKeys.size(); ++i)
        {
          const Index key = static_cast<Index>(metadata.haloKeys[i]);
          const size_t from = static_cast<size_t>(metadata.haloOffsets[i]);
          const size_t to = static_cast<size_t>(metadata.haloOffsets[i + 1]);
          if (from > to || to > metadata.haloIndices.size())
          {
            Alert::Exception()
              << "Invalid shard halo offsets for dimension " << d << "."
              << Alert::Raise;
          }
          for (size_t j = from; j < to; ++j)
            builder.halo(d, key, static_cast<Index>(metadata.haloIndices[j]));
        }
      }

      static void restorePolytopeMap(
          Geometry::Shard& shard,
          size_t d,
          const std::vector<HDF5::U64>& left)
      {
        auto& pmap = shard.getPolytopeMap(d);
        if (pmap.left.size() != left.size())
        {
          Alert::Exception()
            << "Invalid reconstructed shard map size for dimension " << d
            << ": expected " << pmap.left.size() << " entries, got "
            << left.size() << "."
            << Alert::Raise;
        }
        pmap.right.clear();
        for (Index i = 0; i < static_cast<Index>(left.size()); ++i)
        {
          pmap.left[i] = static_cast<Index>(left[i]);
          pmap.right[static_cast<Index>(left[i])] = i;
        }
      }
  };

  /**
   * @brief Prints a distributed MPI mesh to an HDF5 file.
   *
   * Each MPI rank independently saves its local shard to the given HDF5
   * file path.  The process consists of two phases:
   *
   * 1. **Base mesh** — delegated to
   *    `MeshPrinter<FileFormat::HDF5, Context::Local>`, which writes the
   *    canonical `/Mesh/...` datasets (vertices, connectivity, attributes,
   *    transformations).
   *
   * 2. **Shard metadata** — the file is reopened in read-write mode and
   *    the following datasets are appended under `/Shard/`:
   *
   *    | Dataset | Type | Description |
   *    |---------|------|-------------|
   *    | `/Shard/Flags/{d}` | U8[] | Per-entity ownership (0=None, 1=Owned, 2=Ghost) |
   *    | `/Shard/PolytopeMap/{d}/Left` | U64[] | Shard→distributed index |
   *    | `/Shard/PolytopeMap/{d}/Right/Keys` | U64[] | Distributed indices (keys) |
   *    | `/Shard/PolytopeMap/{d}/Right/Values` | U64[] | Shard indices (values) |
   *    | `/Shard/Owner/{d}/Keys` | U64[] | Ghost local indices (keys) |
   *    | `/Shard/Owner/{d}/Values` | U64[] | Owner ranks (values) |
   *    | `/Shard/Halo/{d}/Keys` | U64[] | Owned local indices (keys) |
   *    | `/Shard/Halo/{d}/Offsets` | U64[] | CSR offsets into Indices |
   *    | `/Shard/Halo/{d}/Indices` | U64[] | Flattened neighbor ranks |
   *
   * @note Each rank must be given a rank-specific file path (e.g. via
   *       the callable filename overload on Mesh<Context::MPI>::save).
   *
   * @see MeshPrinter<FileFormat::HDF5, Context::Local>
   * @see Geometry::Shard
   */
  template <>
  class MeshPrinter<FileFormat::HDF5, Context::MPI>
    : public MeshPrinterBase<Context::MPI>
  {
    public:
      /**
       * @brief Distributed context type handled by this printer.
       */
      using ContextType = Context::MPI;
      /**
       * @brief Distributed mesh object type written by this specialization.
       */
      using ObjectType = Geometry::Mesh<ContextType>;
      /**
       * @brief Base printer interface specialization.
       */
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
       * First writes the base mesh via the local printer, then reopens the
       * file and appends shard metadata datasets.
       *
       * @param[in] filename  Path to the HDF5 file for this rank's shard.
       */
      void print(const boost::filesystem::path& filename) override
      {
        const auto& mesh = this->getObject();
        const auto& shard = mesh.getShard();

        // Phase 1: write base mesh topology via local printer
        MeshPrinter<FileFormat::HDF5, Context::Local> localPrinter(shard);
        localPrinter.print(filename);

        // Phase 2: reopen and append shard metadata
        const auto file =
          HDF5::File(H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
        if (!file)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to reopen HDF5 file for shard metadata: " << filename
            << Alert::Raise;
        }

        writeShardMetadata(file.get(), shard);
      }

    private:
      /**
       * @brief Appends all shard metadata to an open HDF5 file.
       * @param[in] file   Open HDF5 file identifier (read-write).
       * @param[in] shard  Shard whose metadata will be written.
       */
      static void writeShardMetadata(hid_t file, const Geometry::Shard& shard)
      {
        const size_t Dmax = shard.getDimension();

        // Create top-level groups: /Shard, /Shard/Flags, etc.
        {
          const auto g = HDF5::Group(
            H5Gcreate2(file, HDF5::Path::Shard, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::Exception()
              << "Failed to create " << HDF5::Path::Shard << " group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(
            file, HDF5::Path::ShardState, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::Exception()
              << "Failed to create " << HDF5::Path::ShardState << " group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(
            file, HDF5::Path::ShardPolytopeMap, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::Exception()
              << "Failed to create " << HDF5::Path::ShardPolytopeMap << " group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(
            file, HDF5::Path::ShardOwner, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::Exception()
              << "Failed to create " << HDF5::Path::ShardOwner << " group."
              << Alert::Raise;
          }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(
            file, HDF5::Path::ShardHalo, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::Exception()
              << "Failed to create " << HDF5::Path::ShardHalo << " group."
              << Alert::Raise;
          }
        }

        for (size_t d = 0; d <= Dmax; ++d)
        {
          writeFlags(file, shard, d);
          writePolytopeMap(file, shard, d);
          writeOwner(file, shard, d);
          writeHalo(file, shard, d);
        }
      }

      /**
       * @brief Writes ownership flags for dimension `d`.
       *
       * Encodes each Flags value as U8: 0=None, 1=Owned, 2=Ghost.
       */
      static void writeFlags(hid_t file, const Geometry::Shard& shard, size_t d)
      {
        const auto& state = shard.getState(d);
        std::vector<HDF5::U8> buf(state.size());
        for (size_t i = 0; i < state.size(); ++i)
        {
          if (state[i] == Geometry::Shard::State::Owned)
            buf[i] = 1;
          else if (state[i] == Geometry::Shard::State::Ghost)
            buf[i] = 2;
          else
            buf[i] = 0;
        }
        HDF5::writeVectorDataset(file, HDF5::shardStatePath(d), buf);
      }

      /**
       * @brief Writes bidirectional polytope index map for dimension `d`.
       */
      static void writePolytopeMap(hid_t file, const Geometry::Shard& shard, size_t d)
      {
        const auto& pmap = shard.getPolytopeMap(d);

        // Per-dimension group
        {
          const auto g =
            HDF5::Group(H5Gcreate2(file, HDF5::shardPolytopeMapGroupPath(d).c_str(),
              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::Exception()
              << "Failed to create " << HDF5::shardPolytopeMapGroupPath(d) << " group."
              << Alert::Raise;
          }
        }

        // Left: shard-to-distributed (vector<Index>)
        {
          std::vector<HDF5::U64> left(pmap.left.size());
          for (size_t i = 0; i < pmap.left.size(); ++i)
            left[i] = static_cast<HDF5::U64>(pmap.left[i]);
          HDF5::writeVectorDataset(file, HDF5::shardPolytopeMapLeftPath(d), left);
        }

        // Right: distributed-to-shard (UnorderedMap)
        {
          const auto rightGroup = HDF5::shardPolytopeMapRightGroupPath(d);
          {
            const auto g = HDF5::Group(H5Gcreate2(
              file, rightGroup.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
            if (!g)
            {
              Alert::Exception()
                << "Failed to create " << rightGroup << " group."
                << Alert::Raise;
            }
          }

          std::vector<HDF5::U64> keys;
          std::vector<HDF5::U64> vals;
          keys.reserve(pmap.right.size());
          vals.reserve(pmap.right.size());
          for (const auto& [k, v] : pmap.right)
          {
            keys.push_back(static_cast<HDF5::U64>(k));
            vals.push_back(static_cast<HDF5::U64>(v));
          }
          HDF5::writeVectorDataset(file, rightGroup + "/Keys", keys);
          HDF5::writeVectorDataset(file, rightGroup + "/Values", vals);
        }
      }

      /**
       * @brief Writes ghost-to-owner rank map for dimension `d`.
       */
      static void writeOwner(hid_t file, const Geometry::Shard& shard, size_t d)
      {
        const auto& owner = shard.getOwner(d);

        {
          const auto g = HDF5::Group(H5Gcreate2(file,
            HDF5::shardOwnerGroupPath(d).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::Exception()
              << "Failed to create " << HDF5::shardOwnerGroupPath(d) << " group."
              << Alert::Raise;
          }
        }

        std::vector<HDF5::U64> keys;
        std::vector<HDF5::U64> vals;
        keys.reserve(owner.size());
        vals.reserve(owner.size());
        for (const auto& [k, v] : owner)
        {
          keys.push_back(static_cast<HDF5::U64>(k));
          vals.push_back(static_cast<HDF5::U64>(v));
        }
        HDF5::writeVectorDataset(file, HDF5::shardOwnerGroupPath(d) + "/Keys", keys);
        HDF5::writeVectorDataset(file, HDF5::shardOwnerGroupPath(d) + "/Values", vals);
      }

      /**
       * @brief Writes owned-to-halo CSR for dimension `d`.
       *
       * The halo is serialized as three arrays:
       * - Keys: owned local indices that appear in other ranks
       * - Offsets: CSR offsets into Indices (length = Keys.size() + 1)
       * - Indices: flattened neighbor rank sets
       */
      static void writeHalo(hid_t file, const Geometry::Shard& shard, size_t d)
      {
        const auto& halo = shard.getHalo(d);

        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::shardHaloGroupPath(d).c_str(),
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g)
          {
            Alert::Exception()
              << "Failed to create " << HDF5::shardHaloGroupPath(d) << " group."
              << Alert::Raise;
          }
        }

        std::vector<HDF5::U64> keys;
        std::vector<HDF5::U64> offsets;
        std::vector<HDF5::U64> indices;

        keys.reserve(halo.size());
        offsets.reserve(halo.size() + 1);
        offsets.push_back(0);

        for (const auto& [k, ranks] : halo)
        {
          keys.push_back(static_cast<HDF5::U64>(k));
          for (const auto& r : ranks)
            indices.push_back(static_cast<HDF5::U64>(r));
          offsets.push_back(static_cast<HDF5::U64>(indices.size()));
        }

        const auto base = HDF5::shardHaloGroupPath(d);
        HDF5::writeVectorDataset(file, base + "/Keys", keys);
        HDF5::writeVectorDataset(file, base + "/Offsets", offsets);
        HDF5::writeVectorDataset(file, base + "/Indices", indices);
      }
  };

  // ---- MPI GridFunction IO --------------------------------------------------
  //
  // Grid function load/print for MPI meshes reuses the local specializations
  // defined in Rodin/IO/HDF5.h (included above).  The local
  // GridFunctionLoader<FileFormat::HDF5, FES, Math::Vector<Scalar>> and
  // GridFunctionPrinter<FileFormat::HDF5, FES, Math::Vector<Scalar>>
  // already handle per-rank DOF vector persistence, so no MPI-specific
  // overrides are needed here.
}

#endif
