/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file HDF5.h
 * @brief HDF5-based mesh and grid function persistence for Rodin.
 *
 * This file provides the HDF5 storage backend for Rodin meshes and grid
 * functions using the official HDF5 C API. It implements:
 *
 * - Full mesh persistence and direct reconstruction from disk, including
 *   vertex geometry, full connectivity (entities, incidences, state),
 *   attributes, and polytope transformations.
 * - Standalone grid function field-file serialization (one `.h5` per grid
 *   function), without mesh or finite element space embedding.
 * - Utility functions for XDMF mixed-topology stream generation and
 *   minimal visualization mesh export, used by the IO::XDMF writer.
 *   These helpers (`writeXDMFMesh`, `writeXDMFTopology`, `writeXDMFVertices`,
 *   `writeXDMFRegionAttribute`, `writeXDMFNodeAttribute`,
 *   `writeXDMFCellAttribute`) write only the minimal HDF5 datasets needed
 *   by ParaView/XDMF, and are completely separate from the canonical
 *   persistence path.
 *
 * ## HDF5 Mesh Layout
 *
 * ```
 * /Mesh
 *   /Meta
 *     SpaceDimension                  scalar (U64)
 *   /Geometry
 *     Vertices                        [nv × sdim] (F64)
 *   /Connectivity
 *     /Meta
 *       MaximalDimension              scalar (U64)
 *       Dimension                     scalar (U64)
 *     /Counts
 *       ByDimension                   [Dmax+1] (U64)
 *       ByGeometry                    [num_types] (U64)
 *     /Entities/{d}
 *       Types                         [n_d] (I32)
 *       Offsets                        [n_d+1] (U64)
 *       Indices                        [nnz] (U64)
 *     /State
 *       Present                       [(Dmax+1)×(Dmax+1)] (U8)
 *       Dirty                         [(Dmax+1)×(Dmax+1)] (U8)
 *     /Incidence/{d}_{dp}
 *       Offsets                        [n_d+1] (U64)
 *       Indices                        [nnz] (U64)
 *   /Attributes/{d}                   [n_d] (U64)
 *   /Transformations/{d}
 *     Kind                            [n_d] (I32)
 * ```
 *
 * ## XDMF Visualization Datasets (written by IO::XDMF via writeXDMFMesh)
 *
 * ```
 * /Mesh                               (visualization-only file, no canonical persistence)
 *   /Geometry
 *     Vertices                          [nv × sdim] (F64)
 *   /XDMF
 *     Topology                          [mixed_stream_size] (U64)
 *     TopologySize                      scalar (U64)
 *   /Attributes/{d}                     [n_d] (U64)
 * ```
 *
 * ## HDF5 GridFunction Layout (standalone field file)
 *
 * ```
 * /GridFunction
 *   /Meta
 *     Size                            scalar (U64)
 *     Dimension                       scalar (U64)
 *   /Values
 *     Data                            [N] (F64)
 * ```
 *
 * The IO classes follow the MeshLoader/MeshPrinter and
 * GridFunctionLoader/GridFunctionPrinter API patterns established by the
 * MFEM and MEDIT format specializations.
 *
 * @see IO::MeshLoader, IO::MeshPrinter, IO::GridFunctionLoader,
 *      IO::GridFunctionPrinter, IO::XDMF
 */
#ifndef RODIN_IO_HDF5_H
#define RODIN_IO_HDF5_H

#include <hdf5.h>

#include <array>
#include <limits>
#include <string>
#include <vector>
#include <type_traits>
#include <boost/filesystem/path.hpp>

#include "ForwardDecls.h"

#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Context/Local.h"
#include "Rodin/Geometry/PointCloud.h"
#include "Rodin/Geometry/Connectivity.h"
#include "Rodin/Geometry/AttributeIndex.h"
#include "Rodin/Geometry/PolytopeTransformationIndex.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/IO/MeshLoader.h"
#include "Rodin/IO/MeshPrinter.h"
#include "Rodin/IO/GridFunctionLoader.h"
#include "Rodin/IO/GridFunctionPrinter.h"
#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::IO
{
  /**
   * @brief HDF5 utility namespace providing type aliases, path constants,
   *        RAII handle wrappers, and low-level dataset read/write helpers.
   *
   * This namespace contains all HDF5-specific helpers used by the MeshLoader,
   * MeshPrinter, GridFunctionLoader, and GridFunctionPrinter specializations
   * for the HDF5 file format. It wraps the HDF5 C API in a type-safe,
   * RAII-managed interface.
   */
  namespace HDF5
  {
    /// @brief Unsigned 64-bit integer type used for HDF5 dataset storage.
    using U64 = unsigned long long;

    /// @brief Signed 32-bit integer type used for HDF5 dataset storage.
    using I32 = int;

    /// @brief 64-bit floating point type used for HDF5 dataset storage.
    using F64 = double;

    /// @brief Unsigned 8-bit integer type used for boolean flag datasets.
    using U8  = unsigned char;

    /**
     * @brief Sentinel value used to represent the absence of an attribute
     *        in HDF5 attribute datasets.
     *
     * Polytopes without an assigned attribute are stored as this value
     * in the `/Mesh/Attributes/{d}` arrays. On load, entries equal to
     * this marker are left unset.
     */
    static constexpr U64 NullAttributeMarker = std::numeric_limits<U64>::max();

    /**
     * @brief HDF5 dataset path constants for the mesh and grid function
     *        storage layouts.
     *
     * All paths are absolute within their respective HDF5 file. Mesh
     * datasets live under `/Mesh/...` and grid function datasets under
     * `/GridFunction/...`. The XDMF-specific derived topology lives
     * under `/Mesh/XDMF/...`.
     */
    namespace Path
    {
      static constexpr const char* Mesh = "/Mesh";  ///< Root mesh group.

      static constexpr const char* MeshMeta = "/Mesh/Meta";  ///< Mesh metadata group.
      static constexpr const char* MeshMetaSpaceDimension = "/Mesh/Meta/SpaceDimension";  ///< Spatial embedding dimension (scalar).

      static constexpr const char* MeshGeometry = "/Mesh/Geometry";  ///< Geometry group.
      static constexpr const char* MeshGeometryVertices = "/Mesh/Geometry/Vertices";  ///< Vertex coordinate matrix [nv × sdim].

      static constexpr const char* MeshConnectivity = "/Mesh/Connectivity";  ///< Connectivity root group.
      static constexpr const char* MeshConnectivityMeta = "/Mesh/Connectivity/Meta";  ///< Connectivity metadata group.
      static constexpr const char* MeshConnectivityMetaMaximalDimension = "/Mesh/Connectivity/Meta/MaximalDimension";  ///< Maximal topological dimension.
      static constexpr const char* MeshConnectivityMetaDimension = "/Mesh/Connectivity/Meta/Dimension";  ///< Top-dimensional cell dimension.

      static constexpr const char* MeshConnectivityCounts = "/Mesh/Connectivity/Counts";  ///< Count arrays group.
      static constexpr const char* MeshConnectivityCountsByDimension = "/Mesh/Connectivity/Counts/ByDimension";  ///< Entity counts per dimension [Dmax+1].
      static constexpr const char* MeshConnectivityCountsByGeometry = "/Mesh/Connectivity/Counts/ByGeometry";  ///< Entity counts per geometry type.

      static constexpr const char* MeshConnectivityEntities = "/Mesh/Connectivity/Entities";  ///< CSR entity groups per dimension.
      static constexpr const char* MeshConnectivityState = "/Mesh/Connectivity/State";  ///< Incidence state group.
      static constexpr const char* MeshConnectivityStatePresent = "/Mesh/Connectivity/State/Present";  ///< Incidence present flags [(Dmax+1)²].
      static constexpr const char* MeshConnectivityStateDirty = "/Mesh/Connectivity/State/Dirty";  ///< Incidence dirty flags [(Dmax+1)²].
      static constexpr const char* MeshConnectivityIncidence = "/Mesh/Connectivity/Incidence";  ///< CSR incidence groups.

      static constexpr const char* MeshAttributes = "/Mesh/Attributes";  ///< Polytope attribute arrays group.
      static constexpr const char* MeshTransformations = "/Mesh/Transformations";  ///< Polytope transformation groups.

      static constexpr const char* MeshXDMF = "/Mesh/XDMF";  ///< XDMF-specific derived data group.
      static constexpr const char* MeshXDMFTopology = "/Mesh/XDMF/Topology";  ///< XDMF mixed topology stream.
      static constexpr const char* MeshXDMFTopologySize = "/Mesh/XDMF/TopologySize";  ///< Length of the mixed topology stream.

      static constexpr const char* Shard = "/Shard";  ///< Root shard metadata group.
      static constexpr const char* ShardFlags = "/Shard/Flags";  ///< Ownership flags per dimension.
      static constexpr const char* ShardPolytopeMap = "/Shard/PolytopeMap";  ///< Index maps per dimension.
      static constexpr const char* ShardOwner = "/Shard/Owner";  ///< Ghost-to-owner map per dimension.
      static constexpr const char* ShardHalo = "/Shard/Halo";  ///< Owned-to-halo map per dimension.

      static constexpr const char* GridFunction = "/GridFunction";  ///< Root grid function group.
      static constexpr const char* GridFunctionMeta = "/GridFunction/Meta";  ///< Grid function metadata group.
      static constexpr const char* GridFunctionMetaName = "/GridFunction/Meta/Name";  ///< Optional name string.
      static constexpr const char* GridFunctionMetaSize = "/GridFunction/Meta/Size";  ///< Number of DOFs (scalar).
      static constexpr const char* GridFunctionMetaDimension = "/GridFunction/Meta/Dimension";  ///< Vector dimension (scalar).
      static constexpr const char* GridFunctionValues = "/GridFunction/Values";  ///< Values group.
      static constexpr const char* GridFunctionValuesData = "/GridFunction/Values/Data";  ///< Raw DOF data vector.
    }

    /**
     * @brief RAII wrapper for HDF5 identifier handles.
     *
     * Manages the lifetime of an HDF5 `hid_t` identifier by calling the
     * appropriate close function (e.g. `H5Fclose`, `H5Gclose`, `H5Dclose`,
     * `H5Sclose`) on destruction or reset. Move-only; copying is disabled.
     *
     * Use the factory functions File(), Group(), DataSet(), and Space()
     * to construct handles with the correct close function.
     *
     * ## Usage Example
     * ```cpp
     * auto file = HDF5::File(H5Fcreate("out.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
     * if (!file)
     *   // handle error ...
     *
     * auto group = HDF5::Group(H5Gcreate2(file.get(), "/Data", ...));
     * // group is automatically closed when it goes out of scope.
     * ```
     */
    class Handle
    {
      public:
        /// @brief Default constructor; creates an invalid handle.
        Handle()
          : m_id(-1),
            m_close(nullptr)
        {}

        /**
         * @brief Constructs a handle wrapping the given HDF5 identifier.
         * @param[in] id       HDF5 identifier (`hid_t`) to manage.
         * @param[in] closeFn  HDF5 close function matching the identifier type.
         */
        Handle(hid_t id, herr_t (*closeFn)(hid_t))
          : m_id(id),
            m_close(closeFn)
        {}

        Handle(const Handle&) = delete;
        Handle& operator=(const Handle&) = delete;

        Handle(Handle&& other) noexcept
          : m_id(other.m_id),
            m_close(other.m_close)
        {
          other.m_id = -1;
          other.m_close = nullptr;
        }

        Handle& operator=(Handle&& other) noexcept
        {
          if (this != &other)
          {
            reset();
            m_id = other.m_id;
            m_close = other.m_close;
            other.m_id = -1;
            other.m_close = nullptr;
          }
          return *this;
        }

        ~Handle()
        {
          reset();
        }

        /// @brief Releases the managed identifier by calling its close function.
        void reset()
        {
          if (m_id >= 0 && m_close)
            m_close(m_id);
          m_id = -1;
          m_close = nullptr;
        }

        /**
         * @brief Returns the raw HDF5 identifier.
         * @returns The managed `hid_t` value, or a negative value if invalid.
         */
        hid_t get() const
        {
          return m_id;
        }

        /**
         * @brief Tests whether this handle holds a valid HDF5 identifier.
         * @returns `true` if the identifier is non-negative (valid).
         */
        explicit operator bool() const
        {
          return m_id >= 0;
        }

      private:
        hid_t m_id;
        herr_t (*m_close)(hid_t);
    };

    /**
     * @brief Creates a Handle for an HDF5 file identifier.
     * @param[in] id  File identifier returned by `H5Fcreate` or `H5Fopen`.
     * @returns Handle that calls `H5Fclose` on destruction.
     */
    inline
    Handle File(hid_t id)
    {
      return Handle(id, H5Fclose);
    }

    /**
     * @brief Creates a Handle for an HDF5 group identifier.
     * @param[in] id  Group identifier returned by `H5Gcreate2` or `H5Gopen2`.
     * @returns Handle that calls `H5Gclose` on destruction.
     */
    inline
    Handle Group(hid_t id)
    {
      return Handle(id, H5Gclose);
    }

    /**
     * @brief Creates a Handle for an HDF5 dataset identifier.
     * @param[in] id  Dataset identifier returned by `H5Dcreate2` or `H5Dopen2`.
     * @returns Handle that calls `H5Dclose` on destruction.
     */
    inline
    Handle DataSet(hid_t id)
    {
      return Handle(id, H5Dclose);
    }

    /**
     * @brief Creates a Handle for an HDF5 dataspace identifier.
     * @param[in] id  Dataspace identifier returned by `H5Screate_simple` etc.
     * @returns Handle that calls `H5Sclose` on destruction.
     */
    inline
    Handle Space(hid_t id)
    {
      return Handle(id, H5Sclose);
    }

    /**
     * @brief Returns the HDF5 native memory type corresponding to `T`.
     * @tparam T  One of U64, I32, F64, or U8.
     * @returns The HDF5 native type identifier for use with `H5Dread`/`H5Dwrite`.
     */
    template <class T>
    hid_t getNativeType();

    template <>
    inline
    hid_t getNativeType<U64>()
    {
      return H5T_NATIVE_ULLONG;
    }

    template <>
    inline
    hid_t getNativeType<I32>()
    {
      return H5T_NATIVE_INT;
    }

    template <>
    inline
    hid_t getNativeType<F64>()
    {
      return H5T_NATIVE_DOUBLE;
    }

    template <>
    inline
    hid_t getNativeType<U8>()
    {
      return H5T_NATIVE_UCHAR;
    }

    /**
     * @brief Tests whether an HDF5 link exists at the given path.
     * @param[in] loc   HDF5 location identifier (file or group).
     * @param[in] path  Absolute or relative path to test.
     * @returns `true` if the link exists.
     */
    inline
    bool exists(hid_t loc, const std::string& path)
    {
      return H5Lexists(loc, path.c_str(), H5P_DEFAULT) > 0;
    }

    /**
     * @brief Returns the HDF5 group path for entities of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Mesh/Connectivity/Entities/2"`.
     */
    inline
    std::string entityGroupPath(size_t d)
    {
      return std::string(Path::MeshConnectivityEntities) + "/" + std::to_string(d);
    }

    /**
     * @brief Returns the dataset path for entity type codes of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Mesh/Connectivity/Entities/2/Types"`.
     */
    inline
    std::string entityTypesPath(size_t d)
    {
      return entityGroupPath(d) + "/Types";
    }

    /**
     * @brief Returns the dataset path for entity CSR offsets of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Mesh/Connectivity/Entities/2/Offsets"`.
     */
    inline
    std::string entityOffsetsPath(size_t d)
    {
      return entityGroupPath(d) + "/Offsets";
    }

    /**
     * @brief Returns the dataset path for entity CSR vertex indices of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Mesh/Connectivity/Entities/2/Indices"`.
     */
    inline
    std::string entityIndicesPath(size_t d)
    {
      return entityGroupPath(d) + "/Indices";
    }

    /**
     * @brief Returns the HDF5 group path for the incidence relation `d → dp`.
     * @param[in] d   Source topological dimension.
     * @param[in] dp  Target topological dimension.
     * @returns Path string, e.g. `"/Mesh/Connectivity/Incidence/2_0"`.
     */
    inline
    std::string incidenceGroupPath(size_t d, size_t dp)
    {
      return std::string(Path::MeshConnectivityIncidence) + "/" + std::to_string(d) + "_" + std::to_string(dp);
    }

    /**
     * @brief Returns the dataset path for incidence CSR offsets.
     * @param[in] d   Source topological dimension.
     * @param[in] dp  Target topological dimension.
     * @returns Path string, e.g. `"/Mesh/Connectivity/Incidence/2_0/Offsets"`.
     */
    inline
    std::string incidenceOffsetsPath(size_t d, size_t dp)
    {
      return incidenceGroupPath(d, dp) + "/Offsets";
    }

    /**
     * @brief Returns the dataset path for incidence CSR indices.
     * @param[in] d   Source topological dimension.
     * @param[in] dp  Target topological dimension.
     * @returns Path string, e.g. `"/Mesh/Connectivity/Incidence/2_0/Indices"`.
     */
    inline
    std::string incidenceIndicesPath(size_t d, size_t dp)
    {
      return incidenceGroupPath(d, dp) + "/Indices";
    }

    /**
     * @brief Returns the dataset path for polytope attributes of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Mesh/Attributes/2"`.
     */
    inline
    std::string attributePath(size_t d)
    {
      return std::string(Path::MeshAttributes) + "/" + std::to_string(d);
    }

    /**
     * @brief Returns the HDF5 group path for polytope transformations of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Mesh/Transformations/2"`.
     */
    inline
    std::string transformationGroupPath(size_t d)
    {
      return std::string(Path::MeshTransformations) + "/" + std::to_string(d);
    }

    /**
     * @brief Returns the dataset path for transformation kind values of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Mesh/Transformations/2/Kind"`.
     */
    inline
    std::string transformationKindPath(size_t d)
    {
      return transformationGroupPath(d) + "/Kind";
    }

    // ---- Shard metadata path helpers ----------------------------------------

    /**
     * @brief Returns the dataset path for shard ownership flags of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Shard/Flags/2"`.
     */
    inline
    std::string shardFlagsPath(size_t d)
    {
      return std::string(Path::ShardFlags) + "/" + std::to_string(d);
    }

    /**
     * @brief Returns the group path for shard polytope map of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Shard/PolytopeMap/2"`.
     */
    inline
    std::string shardPolytopeMapGroupPath(size_t d)
    {
      return std::string(Path::ShardPolytopeMap) + "/" + std::to_string(d);
    }

    /**
     * @brief Returns the dataset path for shard-to-distributed index mapping.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Shard/PolytopeMap/2/Left"`.
     */
    inline
    std::string shardPolytopeMapLeftPath(size_t d)
    {
      return shardPolytopeMapGroupPath(d) + "/Left";
    }

    /**
     * @brief Returns the group path for distributed-to-shard index mapping.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Shard/PolytopeMap/2/Right"`.
     */
    inline
    std::string shardPolytopeMapRightGroupPath(size_t d)
    {
      return shardPolytopeMapGroupPath(d) + "/Right";
    }

    /**
     * @brief Returns the group path for ghost-to-owner map of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Shard/Owner/2"`.
     */
    inline
    std::string shardOwnerGroupPath(size_t d)
    {
      return std::string(Path::ShardOwner) + "/" + std::to_string(d);
    }

    /**
     * @brief Returns the group path for owned-to-halo map of dimension `d`.
     * @param[in] d  Topological dimension.
     * @returns Path string, e.g. `"/Shard/Halo/2"`.
     */
    inline
    std::string shardHaloGroupPath(size_t d)
    {
      return std::string(Path::ShardHalo) + "/" + std::to_string(d);
    }

    // ---- end shard metadata path helpers ------------------------------------

    /**
     * @brief Returns the size of the geometry-count array.
     * @returns Number of slots needed to index by `Polytope::Type`.
     */
    inline
    size_t getGeometryCountArraySize()
    {
      return Geometry::Polytope::Types.size();
    }

    /**
     * @brief Maps a Rodin polytope type to its XDMF mixed-topology type id.
     * @param[in] t  Polytope geometry type.
     * @returns XDMF type id (e.g. Triangle=4, Quadrilateral=5, Tetrahedron=6).
     *
     * @see <a href="https://www.xdmf.org/index.php/XDMF_Model_and_Format">
     *      XDMF Model and Format Specification</a>
     */
    inline
    U64 getXDMFMixedTopologyId(Geometry::Polytope::Type t)
    {
      using PT = Geometry::Polytope::Type;
      switch (t)
      {
        case PT::Point:         return 1;
        case PT::Segment:       return 2;
        case PT::Triangle:      return 4;
        case PT::Quadrilateral: return 5;
        case PT::Tetrahedron:   return 6;
        case PT::Wedge:         return 8;
        case PT::Hexahedron:    return 9;
      }

      Alert::Exception()
        << "Unsupported polytope type for XDMF mixed topology."
        << Alert::Raise;
    }

    /**
     * @brief Computes the total length of the XDMF mixed-topology stream
     *        for the given mesh.
     * @param[in] mesh  Mesh whose cells contribute to the stream.
     * @returns Total number of entries in the mixed topology array.
     *
     * Each cell contributes: 1 (type id) + optional count + number of vertices.
     */
    inline
    size_t getXDMFMixedTopologySize(const Geometry::MeshBase& mesh)
    {
      size_t size = 0;
      for (auto it = mesh.getCell(); !it.end(); ++it)
      {
        const auto geometry = it->getGeometry();
        const size_t nv = it->getVertices().size();
        size += 1 + nv;
        if (geometry == Geometry::Polytope::Type::Segment)
          size += 1;
      }
      return size;
    }

    /**
     * @brief Reads a 1D HDF5 dataset into a `std::vector`.
     * @tparam T    Element type (must be one of U64, I32, F64, U8).
     * @param[in] file  Open HDF5 file identifier.
     * @param[in] path  Absolute dataset path within the file.
     * @returns Vector containing all elements of the dataset.
     */
    template <class T>
    std::vector<T> readVectorDataset(hid_t file, const std::string& path)
    {
      const auto dset = DataSet(H5Dopen2(file, path.c_str(), H5P_DEFAULT));
      if (!dset)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataset: " << path
          << Alert::Raise;
      }

      const auto space = Space(H5Dget_space(dset.get()));
      if (!space)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataspace: " << path
          << Alert::Raise;
      }

      const auto count = static_cast<size_t>(H5Sget_simple_extent_npoints(space.get()));
      std::vector<T> values(count);
      if (count > 0)
      {
        const auto status = H5Dread(
            dset.get(),
            getNativeType<T>(),
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            values.data());
        if (status < 0)
        {
          Alert::Exception()
            << "Failed to read HDF5 dataset: " << path
            << Alert::Raise;
        }
      }
      return values;
    }

    /**
     * @brief Reads a scalar (single-element) HDF5 dataset.
     * @tparam T    Element type (must be one of U64, I32, F64, U8).
     * @param[in] file  Open HDF5 file identifier.
     * @param[in] path  Absolute dataset path within the file.
     * @returns The single scalar value stored in the dataset.
     */
    template <class T>
    T readScalarDataset(hid_t file, const std::string& path)
    {
      const auto values = readVectorDataset<T>(file, path);
      if (values.size() != 1)
      {
        Alert::Exception()
          << "Expected scalar HDF5 dataset at path: " << path
          << Alert::Raise;
      }
      return values[0];
    }

    /**
     * @brief Writes a `std::vector` as a 1D HDF5 dataset.
     * @tparam T        Element type (must be one of U64, I32, F64, U8).
     * @param[in] file    Open HDF5 file identifier (writable).
     * @param[in] path    Absolute dataset path to create.
     * @param[in] values  Data to write.
     */
    template <class T>
    void writeVectorDataset(hid_t file, const std::string& path, const std::vector<T>& values)
    {
      const hsize_t dims[1] = { static_cast<hsize_t>(values.size()) };
      const auto space = Space(H5Screate_simple(1, dims, nullptr));
      if (!space)
      {
        Alert::Exception()
          << "Failed to create HDF5 dataspace for dataset: " << path
          << Alert::Raise;
      }

      const auto dset = DataSet(H5Dcreate2(
          file,
          path.c_str(),
          getNativeType<T>(),
          space.get(),
          H5P_DEFAULT,
          H5P_DEFAULT,
          H5P_DEFAULT));
      if (!dset)
      {
        Alert::Exception()
          << "Failed to create HDF5 dataset: " << path
          << Alert::Raise;
      }

      if (!values.empty())
      {
        const auto status = H5Dwrite(
            dset.get(),
            getNativeType<T>(),
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            values.data());
        if (status < 0)
        {
          Alert::Exception()
            << "Failed to write HDF5 dataset: " << path
            << Alert::Raise;
        }
      }
    }

    /**
     * @brief Writes a single scalar value as a 1D dataset of length 1.
     * @tparam T       Element type.
     * @param[in] file   Open HDF5 file identifier (writable).
     * @param[in] path   Absolute dataset path to create.
     * @param[in] value  Scalar value to write.
     */
    template <class T>
    void writeScalarDataset(hid_t file, const std::string& path, const T& value)
    {
      writeVectorDataset(file, path, std::vector<T>{ value });
    }

    /**
     * @brief Writes a row-major packed vector as a 2D HDF5 dataset.
     * @tparam T        Element type.
     * @param[in] file    Open HDF5 file identifier (writable).
     * @param[in] path    Absolute dataset path to create.
     * @param[in] values  Row-major packed matrix data (`rows × cols` elements).
     * @param[in] rows    Number of rows.
     * @param[in] cols    Number of columns.
     */
    template <class T>
    void writeMatrixDataset(
        hid_t file,
        const std::string& path,
        const std::vector<T>& values,
        hsize_t rows,
        hsize_t cols)
    {
      if (values.size() != static_cast<size_t>(rows * cols))
      {
        Alert::Exception()
          << "Invalid HDF5 matrix payload size for dataset: " << path
          << Alert::Raise;
      }

      const hsize_t dims[2] = { rows, cols };
      const auto space = Space(H5Screate_simple(2, dims, nullptr));
      if (!space)
      {
        Alert::Exception()
          << "Failed to create HDF5 dataspace for dataset: " << path
          << Alert::Raise;
      }

      const auto dset = DataSet(H5Dcreate2(
          file,
          path.c_str(),
          getNativeType<T>(),
          space.get(),
          H5P_DEFAULT,
          H5P_DEFAULT,
          H5P_DEFAULT));
      if (!dset)
      {
        Alert::Exception()
          << "Failed to create HDF5 dataset: " << path
          << Alert::Raise;
      }

      if (!values.empty())
      {
        const auto status = H5Dwrite(
            dset.get(),
            getNativeType<T>(),
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            values.data());
        if (status < 0)
        {
          Alert::Exception()
            << "Failed to write HDF5 dataset: " << path
            << Alert::Raise;
        }
      }
    }

    /**
     * @brief Reads the shape (rows, cols) of a rank-2 HDF5 dataset.
     * @param[in] file  Open HDF5 file identifier.
     * @param[in] path  Absolute dataset path.
     * @returns Pair of `(rows, cols)`.
     */
    inline
    std::pair<hsize_t, hsize_t> readMatrixShape(hid_t file, const std::string& path)
    {
      const auto dset = DataSet(H5Dopen2(file, path.c_str(), H5P_DEFAULT));
      if (!dset)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataset: " << path
          << Alert::Raise;
      }

      const auto space = Space(H5Dget_space(dset.get()));
      if (!space)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataspace: " << path
          << Alert::Raise;
      }

      hsize_t dims[2] = { 0, 0 };
      const int rank = H5Sget_simple_extent_ndims(space.get());
      if (rank != 2)
      {
        Alert::Exception()
          << "Expected rank-2 HDF5 dataset at path: " << path
          << Alert::Raise;
      }

      const auto status = H5Sget_simple_extent_dims(space.get(), dims, nullptr);
      if (status != 2)
      {
        Alert::Exception()
          << "Failed to read matrix shape for HDF5 dataset: " << path
          << Alert::Raise;
      }

      return { dims[0], dims[1] };
    }

    /**
     * @brief Reads the full shape of an HDF5 dataset as a vector of dimensions.
     * @param[in] file  Open HDF5 file identifier.
     * @param[in] path  Absolute dataset path.
     * @returns Vector of dimension sizes (length equals the dataset rank).
     */
    inline
    std::vector<hsize_t> readDatasetShape(hid_t file, const std::string& path)
    {
      const auto dset = DataSet(H5Dopen2(file, path.c_str(), H5P_DEFAULT));
      if (!dset)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataset: " << path
          << Alert::Raise;
      }

      const auto space = Space(H5Dget_space(dset.get()));
      if (!space)
      {
        Alert::Exception()
          << "Failed to open HDF5 dataspace: " << path
          << Alert::Raise;
      }

      const int rank = H5Sget_simple_extent_ndims(space.get());
      if (rank < 0)
      {
        Alert::Exception()
          << "Failed to read rank of HDF5 dataset: " << path
          << Alert::Raise;
      }

      std::vector<hsize_t> dims(static_cast<size_t>(rank), 0);
      if (rank > 0)
      {
        const auto status = H5Sget_simple_extent_dims(space.get(), dims.data(), nullptr);
        if (status < 0)
        {
          Alert::Exception()
            << "Failed to read shape of HDF5 dataset: " << path
            << Alert::Raise;
        }
      }

      return dims;
    }

    /**
     * @brief Writes the XDMF mixed-topology stream to an open HDF5 file.
     *
     * Creates the `/Mesh/XDMF` group and writes the mixed-topology dataset
     * at `/Mesh/XDMF/Topology` and its size at `/Mesh/XDMF/TopologySize`.
     *
     * @note Only local (sequential) meshes are supported. The mesh parameter
     *       is internally cast to `Geometry::Mesh<Context::Local>`.
     *
     * @param[in] file  Open HDF5 file identifier with write access.
     * @param[in] mesh  Local mesh whose cells provide the topology data.
     */
    inline
    void writeXDMFTopology(hid_t file, const Geometry::MeshBase& mesh)
    {
      {
        const auto g = Group(H5Gcreate2(file, Path::MeshXDMF, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!g)
        {
          Alert::Exception()
            << "Failed to create /Mesh/XDMF group."
            << Alert::Raise;
        }
      }

      const auto& connectivity =
          static_cast<const Geometry::Mesh<Context::Local>&>(mesh).getConnectivity();
      const size_t D = connectivity.getDimension();

      std::vector<U64> topology;
      topology.reserve(getXDMFMixedTopologySize(mesh));

      for (Index i = 0; i < static_cast<Index>(connectivity.getCount(D)); ++i)
      {
        const auto geometry = connectivity.getGeometry(D, i);
        const auto& key = connectivity.getPolytope(D, i);

        topology.push_back(getXDMFMixedTopologyId(geometry));

        if (geometry == Geometry::Polytope::Type::Segment)
          topology.push_back(static_cast<U64>(key.size()));

        for (size_t k = 0; k < key.size(); ++k)
          topology.push_back(static_cast<U64>(key[k]));
      }

      writeVectorDataset(file, Path::MeshXDMFTopology, topology);
      writeScalarDataset(file, Path::MeshXDMFTopologySize, static_cast<U64>(topology.size()));
    }

    /**
     * @brief Writes the XDMF mixed-topology stream to an existing HDF5 file.
     *
     * Opens the given HDF5 file for read/write access, creates the
     * `/Mesh/XDMF` group, and writes the mixed-topology dataset at
     * `/Mesh/XDMF/Topology` and its size at `/Mesh/XDMF/TopologySize`.
     *
     * @note Only local (sequential) meshes are supported. The mesh parameter
     *       is internally cast to `Geometry::Mesh<Context::Local>`.
     *
     * @param[in] filename  Path to an existing HDF5 mesh file.
     * @param[in] mesh      Local mesh whose cells provide the topology data.
     */
    inline
    void writeXDMFTopology(
        const boost::filesystem::path& filename,
        const Geometry::MeshBase& mesh)
    {
      const auto file = File(H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
      if (!file)
      {
        Alert::Exception()
          << "Failed to open HDF5 file for XDMF topology: " << filename
          << Alert::Raise;
      }

      writeXDMFTopology(file.get(), mesh);
    }

    /**
     * @brief Writes vertex coordinate data to an open HDF5 file for XDMF
     *        visualization.
     *
     * Creates the `/Mesh/Geometry` group and writes the vertex coordinate
     * matrix at `/Mesh/Geometry/Vertices` as a `[nv × sdim]` dataset.
     *
     * @param[in] file  Open HDF5 file identifier with write access.
     * @param[in] mesh  Local mesh whose vertex coordinates are exported.
     */
    inline
    void writeXDMFVertices(hid_t file, const Geometry::MeshBase& mesh)
    {
      {
        const auto g = Group(H5Gcreate2(file, Path::MeshGeometry, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!g)
        {
          Alert::Exception()
            << "Failed to create /Mesh/Geometry group."
            << Alert::Raise;
        }
      }

      const auto& localMesh = static_cast<const Geometry::Mesh<Context::Local>&>(mesh);
      const size_t nv = localMesh.getVertexCount();
      const size_t sdim = localMesh.getSpaceDimension();

      std::vector<F64> packed(nv * sdim);
      for (Index i = 0; i < static_cast<Index>(nv); ++i)
      {
        const auto x = localMesh.getVertexCoordinates(i);
        for (size_t d = 0; d < sdim; ++d)
          packed[static_cast<size_t>(i) * sdim + d] = static_cast<F64>(x(d));
      }

      writeMatrixDataset(
          file,
          Path::MeshGeometryVertices,
          packed,
          static_cast<hsize_t>(nv),
          static_cast<hsize_t>(sdim));
    }

    /**
     * @brief Writes polytope region attributes to an open HDF5 file for
     *        XDMF visualization.
     *
     * Creates the `/Mesh/Attributes` group and writes an attribute array
     * for each topological dimension from 0 to the mesh dimension.
     *
     * @param[in] file  Open HDF5 file identifier with write access.
     * @param[in] mesh  Local mesh whose attributes are exported.
     */
    inline
    void writeXDMFRegionAttribute(hid_t file, const Geometry::MeshBase& mesh)
    {
      {
        const auto g = Group(H5Gcreate2(file, Path::MeshAttributes, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!g)
        {
          Alert::Exception()
            << "Failed to create /Mesh/Attributes group."
            << Alert::Raise;
        }
      }

      const auto& localMesh = static_cast<const Geometry::Mesh<Context::Local>&>(mesh);
      const auto& connectivity = localMesh.getConnectivity();
      const size_t D = connectivity.getDimension();

      for (size_t d = 0; d <= D; ++d)
      {
        std::vector<U64> attrs(connectivity.getCount(d), NullAttributeMarker);
        for (Index i = 0; i < static_cast<Index>(connectivity.getCount(d)); ++i)
        {
          if (const auto attr = localMesh.getAttribute(d, i))
            attrs[static_cast<size_t>(i)] = static_cast<U64>(*attr);
        }
        writeVectorDataset(file, attributePath(d), attrs);
      }
    }

    /**
     * @brief Writes a minimal XDMF visualization mesh file.
     *
     * Creates a new HDF5 file containing only the datasets required for
     * XDMF/ParaView visualization:
     * - `/Mesh/Geometry/Vertices` — vertex coordinates
     * - `/Mesh/XDMF/Topology` — XDMF mixed-topology stream
     * - `/Mesh/XDMF/TopologySize` — length of the topology stream
     * - `/Mesh/Attributes/{d}` — polytope region labels
     *
     * This function does **not** write the full canonical Rodin persistence
     * data (connectivity CSR, incidence, state, transformations). Use
     * `MeshPrinter<FileFormat::HDF5>` for full persistence.
     *
     * @param[in] filename  Output HDF5 file path.
     * @param[in] mesh      Local mesh to export for visualization.
     *
     * @see writeXDMFTopology, writeXDMFVertices, writeXDMFRegionAttribute,
     *      MeshPrinter<FileFormat::HDF5, Context::Local>
     */
    inline
    void writeXDMFMesh(
        const boost::filesystem::path& filename,
        const Geometry::MeshBase& mesh)
    {
      const auto file = File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
      if (!file)
      {
        Alert::Exception()
          << "Failed to create HDF5 file for XDMF mesh: " << filename
          << Alert::Raise;
      }

      {
        const auto g = Group(H5Gcreate2(file.get(), Path::Mesh, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!g)
        {
          Alert::Exception()
            << "Failed to create /Mesh group."
            << Alert::Raise;
        }
      }

      writeXDMFVertices(file.get(), mesh);
      writeXDMFTopology(file.get(), mesh);
      writeXDMFRegionAttribute(file.get(), mesh);
    }

    /**
     * @brief Writes a grid function as vertex-centered (nodal) data to an
     *        HDF5 file for XDMF visualization.
     *
     * Evaluates the grid function at each vertex of the mesh and stores
     * the result in the `/GridFunction/Values/Data` dataset. For scalar
     * functions, a flat `[nv]` vector is written; for vector functions, a
     * `[nv × vdim]` matrix is written.
     *
     * @tparam GridFunctionType  Concrete grid function type.
     * @param[in] gf        Grid function to export.
     * @param[in] filename  Output HDF5 file path.
     */
    template <class GridFunctionType>
    void writeXDMFNodeAttribute(const GridFunctionType& gf, const boost::filesystem::path& filename)
    {
      using FESType = typename FormLanguage::Traits<GridFunctionType>::FESType;
      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      const auto& fes = gf.getFiniteElementSpace();
      const auto& mesh = fes.getMesh();

      const size_t nv = mesh.getVertexCount();
      const size_t vdim = gf.getDimension();

      const auto file = HDF5::File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
      if (!file)
      {
        Alert::Exception()
          << "Failed to create HDF5 XDMF node attribute file: " << filename
          << Alert::Raise;
      }

      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunction, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction group."
            << Alert::Raise;
        }
      }
      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunctionMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction/Meta group."
            << Alert::Raise;
        }
      }
      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunctionValues, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction/Values group."
            << Alert::Raise;
        }
      }

      HDF5::writeScalarDataset(file.get(), Path::GridFunctionMetaSize, static_cast<HDF5::U64>(nv));
      HDF5::writeScalarDataset(file.get(), Path::GridFunctionMetaDimension, static_cast<HDF5::U64>(vdim));

      const Geometry::Polytope::Traits ts(Geometry::Polytope::Type::Point);

      if constexpr (std::is_same_v<RangeType, ScalarType>)
      {
        std::vector<HDF5::F64> values(nv);
        for (auto it = mesh.getVertex(); !it.end(); ++it)
        {
          const Index i = it->getIndex();
          const Geometry::Point p(*it, ts.getVertex(0), it->getCoordinates());
          values[static_cast<size_t>(i)] = static_cast<HDF5::F64>(gf(p));
        }

        HDF5::writeVectorDataset(file.get(), Path::GridFunctionValuesData, values);
      }
      else
      {
        std::vector<HDF5::F64> values(nv * vdim);
        for (auto it = mesh.getVertex(); !it.end(); ++it)
        {
          const Index i = it->getIndex();
          const Geometry::Point p(*it, ts.getVertex(0), it->getCoordinates());
          const auto value = gf(p);

          for (size_t c = 0; c < vdim; ++c)
            values[static_cast<size_t>(i) * vdim + c] = static_cast<HDF5::F64>(value[c]);
        }

        HDF5::writeMatrixDataset(
            file.get(),
            Path::GridFunctionValuesData,
            values,
            static_cast<hsize_t>(nv),
            static_cast<hsize_t>(vdim));
      }
    }

    /**
     * @brief Writes a grid function as cell-centered data to an HDF5 file
     *        for XDMF visualization.
     *
     * Evaluates the grid function at the vertices of each cell and averages
     * the values to produce a single value per cell. The result is stored in
     * the `/GridFunction/Values/Data` dataset.
     *
     * @tparam GridFunctionType  Concrete grid function type.
     * @param[in] gf        Grid function to export.
     * @param[in] filename  Output HDF5 file path.
     */
    template <class GridFunctionType>
    void writeXDMFCellAttribute(
        const GridFunctionType& gf, const boost::filesystem::path& filename)
    {
      using FESType = typename FormLanguage::Traits<GridFunctionType>::FESType;
      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      const auto& fes = gf.getFiniteElementSpace();
      const auto& mesh = fes.getMesh();

      const size_t nc = mesh.getCellCount();
      const size_t vdim = gf.getDimension();
      const size_t D = mesh.getDimension();
      const auto& c2v = mesh.getConnectivity().getIncidence(D, 0);

      const auto file = HDF5::File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
      if (!file)
      {
        Alert::Exception()
          << "Failed to create HDF5 XDMF cell attribute file: " << filename
          << Alert::Raise;
      }

      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunction, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction group."
            << Alert::Raise;
        }
      }
      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunctionMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction/Meta group."
            << Alert::Raise;
        }
      }
      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunctionValues, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction/Values group."
            << Alert::Raise;
        }
      }

      HDF5::writeScalarDataset(file.get(), Path::GridFunctionMetaSize, static_cast<HDF5::U64>(nc));
      HDF5::writeScalarDataset(file.get(), Path::GridFunctionMetaDimension, static_cast<HDF5::U64>(vdim));

      const Geometry::Polytope::Traits ts(Geometry::Polytope::Type::Point);

      if constexpr (std::is_same_v<RangeType, ScalarType>)
      {
        std::vector<HDF5::F64> values(nc, 0.0);

        for (Index cell = 0; cell < static_cast<Index>(nc); ++cell)
        {
          const auto& vertices = c2v[cell];
          if (vertices.size() == 0)
          {
            Alert::Exception()
              << "Cell with no vertices encountered during XDMF cell export."
              << Alert::Raise;
          }

          ScalarType accum = ScalarType(0);
          for (size_t k = 0; k < vertices.size(); ++k)
          {
            const auto vit = mesh.getVertex(vertices[k]);
            const Geometry::Point p(*vit, ts.getVertex(0), vit->getCoordinates());
            accum += gf(p);
          }

          values[static_cast<size_t>(cell)] =
            static_cast<HDF5::F64>(accum / static_cast<Real>(vertices.size()));
        }

        HDF5::writeVectorDataset(file.get(), Path::GridFunctionValuesData, values);
      }
      else
      {
        std::vector<HDF5::F64> values(nc * vdim, 0.0);

        for (Index cell = 0; cell < static_cast<Index>(nc); ++cell)
        {
          const auto& vertices = c2v[cell];
          if (vertices.size() == 0)
          {
            Alert::Exception()
              << "Cell with no vertices encountered during XDMF cell export."
              << Alert::Raise;
          }

          std::vector<ScalarType> accum(vdim, ScalarType(0));
          for (size_t k = 0; k < vertices.size(); ++k)
          {
            const auto vit = mesh.getVertex(vertices[k]);
            const Geometry::Point p(*vit, ts.getVertex(0), vit->getCoordinates());
            const auto value = gf(p);

            for (size_t c = 0; c < vdim; ++c)
              accum[c] += value[c];
          }

          for (size_t c = 0; c < vdim; ++c)
          {
            values[static_cast<size_t>(cell) * vdim + c] =
              static_cast<HDF5::F64>(accum[c] / static_cast<Real>(vertices.size()));
          }
        }

        HDF5::writeMatrixDataset(
            file.get(),
            Path::GridFunctionValuesData,
            values,
            static_cast<hsize_t>(nc),
            static_cast<hsize_t>(vdim));
      }
    }
    /**
     * @brief Writes a grid function as vertex-centered (nodal) data to an
     *        HDF5 file for XDMF visualization, evaluated on an explicit mesh.
     *
     * This overload allows the grid function to be evaluated on a mesh that
     * differs from `gf.getFiniteElementSpace().getMesh()`. This is essential
     * for distributed (MPI) visualization where the grid function is defined
     * on a distributed mesh but visualization data must be written per-shard.
     *
     * @tparam GridFunctionType  Concrete grid function type.
     * @param[in] gf        Grid function to export.
     * @param[in] visMesh   Mesh to iterate over for vertex coordinates.
     * @param[in] filename  Output HDF5 file path.
     */
    template <class GridFunctionType>
    void writeXDMFNodeAttribute(
        const GridFunctionType& gf,
        const Geometry::MeshBase& visMesh,
        const boost::filesystem::path& filename)
    {
      using FESType = typename FormLanguage::Traits<GridFunctionType>::FESType;
      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      const size_t nv = visMesh.getVertexCount();
      const size_t vdim = gf.getDimension();

      const auto file = HDF5::File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
      if (!file)
      {
        Alert::Exception()
          << "Failed to create HDF5 XDMF node attribute file: " << filename
          << Alert::Raise;
      }

      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunction, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction group."
            << Alert::Raise;
        }
      }
      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunctionMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction/Meta group."
            << Alert::Raise;
        }
      }
      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunctionValues, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction/Values group."
            << Alert::Raise;
        }
      }

      HDF5::writeScalarDataset(file.get(), Path::GridFunctionMetaSize, static_cast<HDF5::U64>(nv));
      HDF5::writeScalarDataset(file.get(), Path::GridFunctionMetaDimension, static_cast<HDF5::U64>(vdim));

      const Geometry::Polytope::Traits ts(Geometry::Polytope::Type::Point);

      if constexpr (std::is_same_v<RangeType, ScalarType>)
      {
        std::vector<HDF5::F64> values(nv);
        for (auto it = visMesh.getVertex(); !it.end(); ++it)
        {
          const Index i = it->getIndex();
          const Geometry::Point p(*it, ts.getVertex(0), it->getCoordinates());
          values[static_cast<size_t>(i)] = static_cast<HDF5::F64>(gf(p));
        }

        HDF5::writeVectorDataset(file.get(), Path::GridFunctionValuesData, values);
      }
      else
      {
        std::vector<HDF5::F64> values(nv * vdim);
        for (auto it = visMesh.getVertex(); !it.end(); ++it)
        {
          const Index i = it->getIndex();
          const Geometry::Point p(*it, ts.getVertex(0), it->getCoordinates());
          const auto value = gf(p);

          for (size_t c = 0; c < vdim; ++c)
            values[static_cast<size_t>(i) * vdim + c] = static_cast<HDF5::F64>(value[c]);
        }

        HDF5::writeMatrixDataset(
            file.get(),
            Path::GridFunctionValuesData,
            values,
            static_cast<hsize_t>(nv),
            static_cast<hsize_t>(vdim));
      }
    }

    /**
     * @brief Writes a grid function as cell-centered data to an HDF5 file
     *        for XDMF visualization, evaluated on an explicit mesh.
     *
     * This overload allows the grid function to be evaluated on a mesh that
     * differs from `gf.getFiniteElementSpace().getMesh()`. This is essential
     * for distributed (MPI) visualization where the grid function is defined
     * on a distributed mesh but visualization data must be written per-shard.
     *
     * @tparam GridFunctionType  Concrete grid function type.
     * @param[in] gf        Grid function to export.
     * @param[in] visMesh   Mesh to iterate over for cell topology.
     * @param[in] filename  Output HDF5 file path.
     */
    template <class GridFunctionType>
    void writeXDMFCellAttribute(
        const GridFunctionType& gf,
        const Geometry::MeshBase& visMesh,
        const boost::filesystem::path& filename)
    {
      using FESType = typename FormLanguage::Traits<GridFunctionType>::FESType;
      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      const size_t nc = visMesh.getCellCount();
      const size_t vdim = gf.getDimension();
      const size_t D = visMesh.getDimension();
      const auto& c2v = visMesh.getConnectivity().getIncidence(D, 0);

      const auto file = HDF5::File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
      if (!file)
      {
        Alert::Exception()
          << "Failed to create HDF5 XDMF cell attribute file: " << filename
          << Alert::Raise;
      }

      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunction, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction group."
            << Alert::Raise;
        }
      }
      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunctionMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction/Meta group."
            << Alert::Raise;
        }
      }
      {
        const auto group = HDF5::Group(H5Gcreate2(file.get(), Path::GridFunctionValues, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
        {
          Alert::Exception()
            << "Failed to create /GridFunction/Values group."
            << Alert::Raise;
        }
      }

      HDF5::writeScalarDataset(file.get(), Path::GridFunctionMetaSize, static_cast<HDF5::U64>(nc));
      HDF5::writeScalarDataset(file.get(), Path::GridFunctionMetaDimension, static_cast<HDF5::U64>(vdim));

      const Geometry::Polytope::Traits ts(Geometry::Polytope::Type::Point);

      if constexpr (std::is_same_v<RangeType, ScalarType>)
      {
        std::vector<HDF5::F64> values(nc, 0.0);

        for (Index cell = 0; cell < static_cast<Index>(nc); ++cell)
        {
          const auto& vertices = c2v[cell];
          if (vertices.size() == 0)
          {
            Alert::Exception()
              << "Cell with no vertices encountered during XDMF cell export."
              << Alert::Raise;
          }

          ScalarType accum = ScalarType(0);
          for (size_t k = 0; k < vertices.size(); ++k)
          {
            const auto vit = visMesh.getVertex(vertices[k]);
            const Geometry::Point p(*vit, ts.getVertex(0), vit->getCoordinates());
            accum += gf(p);
          }

          values[static_cast<size_t>(cell)] =
            static_cast<HDF5::F64>(accum / static_cast<Real>(vertices.size()));
        }

        HDF5::writeVectorDataset(file.get(), Path::GridFunctionValuesData, values);
      }
      else
      {
        std::vector<HDF5::F64> values(nc * vdim, 0.0);

        for (Index cell = 0; cell < static_cast<Index>(nc); ++cell)
        {
          const auto& vertices = c2v[cell];
          if (vertices.size() == 0)
          {
            Alert::Exception()
              << "Cell with no vertices encountered during XDMF cell export."
              << Alert::Raise;
          }

          std::vector<ScalarType> accum(vdim, ScalarType(0));
          for (size_t k = 0; k < vertices.size(); ++k)
          {
            const auto vit = visMesh.getVertex(vertices[k]);
            const Geometry::Point p(*vit, ts.getVertex(0), vit->getCoordinates());
            const auto value = gf(p);

            for (size_t c = 0; c < vdim; ++c)
              accum[c] += value[c];
          }

          for (size_t c = 0; c < vdim; ++c)
          {
            values[static_cast<size_t>(cell) * vdim + c] =
              static_cast<HDF5::F64>(accum[c] / static_cast<Real>(vertices.size()));
          }
        }

        HDF5::writeMatrixDataset(
            file.get(),
            Path::GridFunctionValuesData,
            values,
            static_cast<hsize_t>(nc),
            static_cast<hsize_t>(vdim));
      }
    }
  }

  /**
   * @brief HDF5 mesh loader specialization for local (sequential) meshes.
   *
   * Reconstructs a complete Geometry::Mesh from an HDF5 file previously
   * created by MeshPrinter<FileFormat::HDF5, Context::Local>. The loader
   * reads vertex geometry, full connectivity (entity CSR and incidence CSR),
   * attributes, and transformation indices.
   *
   * Stream-based loading is not supported; only file-path loading is available.
   *
   * ## Usage Example
   * ```cpp
   * Geometry::Mesh<Context::Local> mesh;
   * IO::MeshLoader<IO::FileFormat::HDF5, Context::Local> loader(mesh);
   * loader.load("output.h5");
   * ```
   *
   * @see MeshPrinter<FileFormat::HDF5, Context::Local>
   */
  template <>
  class MeshLoader<FileFormat::HDF5, Context::Local>
    : public MeshLoaderBase<Context::Local>
  {
    public:
      /// @brief Context type for this loader.
      using ContextType = Context::Local;

      /// @brief Type of mesh object being loaded.
      using ObjectType = Geometry::Mesh<ContextType>;

      /// @brief Parent loader class type.
      using Parent = MeshLoaderBase<ContextType>;

      /**
       * @brief Constructs a mesh loader for the given mesh object.
       * @param[in,out] mesh  Mesh to be populated with data read from HDF5.
       */
      explicit
      MeshLoader(ObjectType& mesh)
        : Parent(mesh)
      {}

      /**
       * @brief Stream-based loading is not supported for HDF5.
       *
       * Always raises an exception. Use the file-path overload instead.
       */
      void load(std::istream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 mesh loading requires file-path based loading."
          << Alert::Raise;
      }

      /**
       * @brief Loads a mesh from an HDF5 file.
       * @param[in] filename  Path to the HDF5 mesh file.
       *
       * Reads all mesh components (vertices, connectivity, incidences,
       * attributes, transformations) and rebuilds the mesh via the
       * Mesh::Build() pipeline.
       */
      void load(const boost::filesystem::path& filename) override
      {
        const auto file = HDF5::File(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!file)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to open HDF5 mesh file: " << filename
            << Alert::Raise;
        }

        auto& mesh = this->getObject();

        const size_t sdim = static_cast<size_t>(
            HDF5::readScalarDataset<HDF5::U64>(file.get(), HDF5::Path::MeshMetaSpaceDimension));

        const auto vertices = this->readVertices(file.get());
        auto connectivity = this->readConnectivity(file.get());
        auto attributes = this->readAttributes(file.get(), connectivity);

        Geometry::PolytopeTransformationIndex transformations;
        transformations.initialize(connectivity.getMaximalDimension());
        for (size_t d = 0; d <= connectivity.getMaximalDimension(); ++d)
          transformations.resize(d, connectivity.getCount(d));

        mesh = ObjectType::Build()
          .initialize(sdim)
          .setVertices(vertices)
          .setConnectivity(std::move(connectivity))
          .setAttributeIndex(std::move(attributes))
          .setTransformationIndex(std::move(transformations))
          .finalize();
      }

    private:
      Geometry::PointCloud readVertices(hid_t file) const
      {
        const auto [nv, sdim] = HDF5::readMatrixShape(file, HDF5::Path::MeshGeometryVertices);
        const auto packed = HDF5::readVectorDataset<HDF5::F64>(file, HDF5::Path::MeshGeometryVertices);
        if (packed.size() != static_cast<size_t>(nv * sdim))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid vertex payload size in HDF5 mesh file."
            << Alert::Raise;
        }

        Geometry::PointCloud vertices(static_cast<std::uint8_t>(sdim), static_cast<size_t>(nv));
        for (size_t i = 0; i < static_cast<size_t>(nv); ++i)
        {
          for (size_t d = 0; d < static_cast<size_t>(sdim); ++d)
            vertices(static_cast<std::uint8_t>(d), i) = packed[i * static_cast<size_t>(sdim) + d];
        }
        return vertices;
      }

      Geometry::Connectivity<ContextType> readConnectivity(hid_t file) const
      {
        Geometry::Connectivity<ContextType> connectivity;

        const size_t Dmax = static_cast<size_t>(
            HDF5::readScalarDataset<HDF5::U64>(file, HDF5::Path::MeshConnectivityMetaMaximalDimension));
        const size_t D = static_cast<size_t>(
            HDF5::readScalarDataset<HDF5::U64>(file, HDF5::Path::MeshConnectivityMetaDimension));

        const auto byDimension = HDF5::readVectorDataset<HDF5::U64>(
            file,
            HDF5::Path::MeshConnectivityCountsByDimension);
        if (byDimension.size() != Dmax + 1)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid /Mesh/Connectivity/Counts/ByDimension dataset size."
            << Alert::Raise;
        }

        connectivity.initialize(Dmax);
        connectivity.nodes(static_cast<size_t>(byDimension[0]));

        for (size_t d = 1; d <= Dmax; ++d)
        {
          const auto types = HDF5::readVectorDataset<HDF5::I32>(file, HDF5::entityTypesPath(d));
          const auto offsets = HDF5::readVectorDataset<HDF5::U64>(file, HDF5::entityOffsetsPath(d));
          const auto indices = HDF5::readVectorDataset<HDF5::U64>(file, HDF5::entityIndicesPath(d));

          if (types.size() != static_cast<size_t>(byDimension[d]))
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Invalid entity type count for dimension " << d << "."
              << Alert::Raise;
          }

          if (offsets.size() != types.size() + 1)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Invalid CSR offsets for entity dimension " << d << "."
              << Alert::Raise;
          }

          connectivity.reserve(d, static_cast<size_t>(byDimension[d]));
          for (size_t i = 0; i < types.size(); ++i)
          {
            const size_t begin = static_cast<size_t>(offsets[i]);
            const size_t end = static_cast<size_t>(offsets[i + 1]);
            if (begin > end || end > indices.size())
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Invalid CSR entity offsets for dimension " << d << "."
                << Alert::Raise;
            }

            Geometry::Polytope::Key key(end - begin);
            for (size_t k = begin; k < end; ++k)
              key[k - begin] = static_cast<Index>(indices[k]);

            connectivity.polytope(
                static_cast<Geometry::Polytope::Type>(types[i]),
                std::move(key));
          }
        }

        const auto present = HDF5::readVectorDataset<HDF5::U8>(file, HDF5::Path::MeshConnectivityStatePresent);
        const auto dirty = HDF5::readVectorDataset<HDF5::U8>(file, HDF5::Path::MeshConnectivityStateDirty);

        if (present.size() != (Dmax + 1) * (Dmax + 1))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid /Mesh/Connectivity/State/Present dataset size."
            << Alert::Raise;
        }

        if (dirty.size() != (Dmax + 1) * (Dmax + 1))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid /Mesh/Connectivity/State/Dirty dataset size."
            << Alert::Raise;
        }

        for (size_t d = 0; d <= Dmax; ++d)
        {
          for (size_t dp = 0; dp <= Dmax; ++dp)
          {
            const size_t flat = d * (Dmax + 1) + dp;
            connectivity.setDirty(d, dp, dirty[flat] != 0);
            if (!present[flat])
              continue;

            const auto offsets = HDF5::readVectorDataset<HDF5::U64>(file, HDF5::incidenceOffsetsPath(d, dp));
            const auto indices = HDF5::readVectorDataset<HDF5::U64>(file, HDF5::incidenceIndicesPath(d, dp));

            if (offsets.size() != connectivity.getCount(d) + 1)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Invalid CSR offsets for incidence "
                << d << " -> " << dp << "."
                << Alert::Raise;
            }

            Geometry::Incidence inc;
            inc.resize(connectivity.getCount(d));
            for (size_t i = 0; i < connectivity.getCount(d); ++i)
            {
              const size_t begin = static_cast<size_t>(offsets[i]);
              const size_t end = static_cast<size_t>(offsets[i + 1]);
              if (begin > end || end > indices.size())
              {
                Alert::MemberFunctionException(*this, __func__)
                  << "Invalid CSR row bounds for incidence "
                  << d << " -> " << dp << "."
                  << Alert::Raise;
              }

              auto& row = inc[i];
              row.reserve(end - begin);
              for (size_t k = begin; k < end; ++k)
                row.push_back(static_cast<Index>(indices[k]));
            }

            connectivity.setIncidence({ d, dp }, std::move(inc));
          }
        }

        return connectivity;
      }

      Geometry::AttributeIndex readAttributes(
          hid_t file,
          const Geometry::Connectivity<ContextType>& connectivity) const
      {
        Geometry::AttributeIndex attrs;
        const size_t D = connectivity.getDimension();
        attrs.initialize(D);

        for (size_t d = 0; d <= D; ++d)
        {
          const size_t count = connectivity.getCount(d);
          attrs.resize(d, count);

          const auto path = HDF5::attributePath(d);
          if (!HDF5::exists(file, path))
            continue;

          const auto values = HDF5::readVectorDataset<HDF5::U64>(file, path);
          if (values.size() != count)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Invalid attribute dataset size for dimension " << d << "."
              << Alert::Raise;
          }

          for (Index i = 0; i < static_cast<Index>(count); ++i)
          {
            if (values[static_cast<size_t>(i)] != HDF5::NullAttributeMarker)
            {
              attrs.set(
                  { d, i },
                  Optional<Geometry::Attribute>(
                    static_cast<Geometry::Attribute>(values[static_cast<size_t>(i)])));
            }
          }
        }

        return attrs;
      }
  };

  /**
   * @brief HDF5 mesh printer specialization for local (sequential) meshes.
   *
   * Serializes a complete Geometry::Mesh to an HDF5 file including:
   * - Vertex geometry (`/Mesh/Geometry/Vertices`)
   * - Full connectivity (entity CSR per dimension, incidence CSR, state flags)
   * - Polytope attributes and transformations
   *
   * This printer is strictly for canonical Rodin mesh persistence. It does
   * not write any XDMF-specific datasets. The XDMF visualization path
   * (IO::XDMF) handles its own derived topology generation separately.
   *
   * Stream-based printing is not supported; only file-path printing is available.
   *
   * ## Usage Example
   * ```cpp
   * const Geometry::Mesh<Context::Local>& mesh = ...;
   * IO::MeshPrinter<IO::FileFormat::HDF5, Context::Local>(mesh)
   *   .print("output.h5");
   * ```
   *
   * @see MeshLoader<FileFormat::HDF5, Context::Local>, IO::XDMF
   */
  template <>
  class MeshPrinter<FileFormat::HDF5, Context::Local>
    : public MeshPrinterBase<Context::Local>
  {
    public:
      /// @brief Context type for this printer.
      using ContextType = Context::Local;

      /// @brief Type of mesh object being printed.
      using ObjectType = Geometry::Mesh<ContextType>;

      /// @brief Parent printer class type.
      using Parent = MeshPrinterBase<ContextType>;

      /**
       * @brief Constructs a mesh printer for the given mesh.
       * @param[in] mesh  Mesh to serialize to HDF5.
       */
      explicit
      MeshPrinter(const ObjectType& mesh)
        : Parent(mesh)
      {}

      /**
       * @brief Stream-based printing is not supported for HDF5.
       *
       * Always raises an exception. Use the file-path overload instead.
       */
      void print(std::ostream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 mesh printing requires file-path based printing."
          << Alert::Raise;
      }

      /**
       * @brief Writes the mesh to an HDF5 file.
       * @param[in] filename  Output HDF5 file path.
       *
       * Creates the file and writes all canonical mesh groups and datasets
       * (geometry, connectivity, attributes, transformations).
       */
      void print(const boost::filesystem::path& filename) override
      {
        const auto& mesh = this->getObject();

        const auto file = HDF5::File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
        if (!file)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to create HDF5 mesh file: " << filename
            << Alert::Raise;
        }

        this->createBaseGroups(file.get());

        HDF5::writeScalarDataset(
            file.get(),
            HDF5::Path::MeshMetaSpaceDimension,
            static_cast<HDF5::U64>(mesh.getSpaceDimension()));

        this->writeVertices(file.get());
        this->writeConnectivity(file.get());
        this->writeAttributes(file.get());
        this->writeTransformations(file.get());
      }

    private:
      void createBaseGroups(hid_t file) const
      {
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::Mesh, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh group." << Alert::Raise; }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::MeshMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh/Meta group." << Alert::Raise; }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::MeshGeometry, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh/Geometry group." << Alert::Raise; }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::MeshConnectivity, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh/Connectivity group." << Alert::Raise; }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::MeshConnectivityMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh/Connectivity/Meta group." << Alert::Raise; }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::MeshConnectivityCounts, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh/Connectivity/Counts group." << Alert::Raise; }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::MeshConnectivityEntities, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh/Connectivity/Entities group." << Alert::Raise; }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::MeshConnectivityState, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh/Connectivity/State group." << Alert::Raise; }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::MeshConnectivityIncidence, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh/Connectivity/Incidence group." << Alert::Raise; }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::MeshAttributes, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh/Attributes group." << Alert::Raise; }
        }
        {
          const auto g = HDF5::Group(H5Gcreate2(file, HDF5::Path::MeshTransformations, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!g) { Alert::Exception() << "Failed to create /Mesh/Transformations group." << Alert::Raise; }
        }
      }

      void writeVertices(hid_t file) const
      {
        const auto& mesh = this->getObject();
        const size_t nv = mesh.getVertexCount();
        const size_t sdim = mesh.getSpaceDimension();

        std::vector<HDF5::F64> packed(nv * sdim);
        for (Index i = 0; i < static_cast<Index>(nv); ++i)
        {
          const auto x = mesh.getVertexCoordinates(i);
          for (size_t d = 0; d < sdim; ++d)
            packed[static_cast<size_t>(i) * sdim + d] = static_cast<HDF5::F64>(x(d));
        }

        HDF5::writeMatrixDataset(
            file,
            HDF5::Path::MeshGeometryVertices,
            packed,
            static_cast<hsize_t>(nv),
            static_cast<hsize_t>(sdim));
      }

      void writeConnectivity(hid_t file) const
      {
        const auto& mesh = this->getObject();
        const auto& connectivity = mesh.getConnectivity();
        const size_t Dmax = connectivity.getMaximalDimension();
        const size_t D = connectivity.getDimension();

        HDF5::writeScalarDataset(file, HDF5::Path::MeshConnectivityMetaMaximalDimension, static_cast<HDF5::U64>(Dmax));
        HDF5::writeScalarDataset(file, HDF5::Path::MeshConnectivityMetaDimension, static_cast<HDF5::U64>(D));

        std::vector<HDF5::U64> byDimension(Dmax + 1, 0);
        for (size_t d = 0; d <= Dmax; ++d)
          byDimension[d] = static_cast<HDF5::U64>(connectivity.getCount(d));
        HDF5::writeVectorDataset(file, HDF5::Path::MeshConnectivityCountsByDimension, byDimension);

        std::vector<HDF5::U64> byGeometry(HDF5::getGeometryCountArraySize(), 0);
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Point)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Point));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Segment)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Segment));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Triangle)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Triangle));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Quadrilateral)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Quadrilateral));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Tetrahedron)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Tetrahedron));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Wedge)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Wedge));
        byGeometry[static_cast<size_t>(Geometry::Polytope::Type::Hexahedron)] =
          static_cast<HDF5::U64>(connectivity.getCount(Geometry::Polytope::Type::Hexahedron));
        HDF5::writeVectorDataset(file, HDF5::Path::MeshConnectivityCountsByGeometry, byGeometry);

        for (size_t d = 1; d <= Dmax; ++d)
        {
          const auto group = HDF5::Group(H5Gcreate2(file, HDF5::entityGroupPath(d).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!group)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create entity group for dimension " << d << "."
              << Alert::Raise;
          }

          std::vector<HDF5::I32> types;
          std::vector<HDF5::U64> offsets;
          std::vector<HDF5::U64> indices;

          types.reserve(connectivity.getCount(d));
          offsets.reserve(connectivity.getCount(d) + 1);
          offsets.push_back(0);

          for (Index i = 0; i < static_cast<Index>(connectivity.getCount(d)); ++i)
          {
            types.push_back(static_cast<HDF5::I32>(connectivity.getGeometry(d, i)));
            const auto& key = connectivity.getPolytope(d, i);
            for (size_t k = 0; k < key.size(); ++k)
              indices.push_back(static_cast<HDF5::U64>(key[k]));
            offsets.push_back(static_cast<HDF5::U64>(indices.size()));
          }

          HDF5::writeVectorDataset(file, HDF5::entityTypesPath(d), types);
          HDF5::writeVectorDataset(file, HDF5::entityOffsetsPath(d), offsets);
          HDF5::writeVectorDataset(file, HDF5::entityIndicesPath(d), indices);
        }

        std::vector<HDF5::U8> present((Dmax + 1) * (Dmax + 1), 0);
        std::vector<HDF5::U8> dirty((Dmax + 1) * (Dmax + 1), 0);
        for (size_t d = 0; d <= Dmax; ++d)
        {
          for (size_t dp = 0; dp <= Dmax; ++dp)
          {
            const size_t flat = d * (Dmax + 1) + dp;
            const auto& inc = connectivity.getIncidence(d, dp);
            present[flat] = static_cast<HDF5::U8>(inc.size() == connectivity.getCount(d) ? 1 : 0);
            dirty[flat] = static_cast<HDF5::U8>(connectivity.isDirty(d, dp) ? 1 : 0);
          }
        }

        HDF5::writeMatrixDataset(file, HDF5::Path::MeshConnectivityStatePresent, present, static_cast<hsize_t>(Dmax + 1), static_cast<hsize_t>(Dmax + 1));
        HDF5::writeMatrixDataset(file, HDF5::Path::MeshConnectivityStateDirty, dirty, static_cast<hsize_t>(Dmax + 1), static_cast<hsize_t>(Dmax + 1));

        for (size_t d = 0; d <= Dmax; ++d)
        {
          for (size_t dp = 0; dp <= Dmax; ++dp)
          {
            const size_t flat = d * (Dmax + 1) + dp;
            if (!present[flat])
              continue;

            const auto group = HDF5::Group(H5Gcreate2(file, HDF5::incidenceGroupPath(d, dp).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
            if (!group)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Failed to create incidence group for "
                << d << " -> " << dp << "."
                << Alert::Raise;
            }

            const auto& inc = connectivity.getIncidence(d, dp);
            std::vector<HDF5::U64> offsets;
            std::vector<HDF5::U64> indices;
            offsets.reserve(connectivity.getCount(d) + 1);
            offsets.push_back(0);

            for (size_t i = 0; i < connectivity.getCount(d); ++i)
            {
              const auto& row = inc[i];
              for (const auto j : row)
                indices.push_back(static_cast<HDF5::U64>(j));
              offsets.push_back(static_cast<HDF5::U64>(indices.size()));
            }

            HDF5::writeVectorDataset(file, HDF5::incidenceOffsetsPath(d, dp), offsets);
            HDF5::writeVectorDataset(file, HDF5::incidenceIndicesPath(d, dp), indices);
          }
        }
      }

      void writeAttributes(hid_t file) const
      {
        const auto& mesh = this->getObject();
        const auto& connectivity = mesh.getConnectivity();
        const size_t D = connectivity.getDimension();

        for (size_t d = 0; d <= D; ++d)
        {
          std::vector<HDF5::U64> attrs(connectivity.getCount(d), HDF5::NullAttributeMarker);
          for (Index i = 0; i < static_cast<Index>(connectivity.getCount(d)); ++i)
          {
            if (const auto attr = mesh.getAttribute(d, i))
              attrs[static_cast<size_t>(i)] = static_cast<HDF5::U64>(*attr);
          }
          HDF5::writeVectorDataset(file, HDF5::attributePath(d), attrs);
        }
      }

      void writeTransformations(hid_t file) const
      {
        const auto& mesh = this->getObject();
        const auto& connectivity = mesh.getConnectivity();
        const size_t Dmax = connectivity.getMaximalDimension();

        for (size_t d = 0; d <= Dmax; ++d)
        {
          const auto group = HDF5::Group(H5Gcreate2(file, HDF5::transformationGroupPath(d).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!group)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create transformation group for dimension " << d << "."
              << Alert::Raise;
          }

          std::vector<HDF5::I32> kind(connectivity.getCount(d), 0);
          HDF5::writeVectorDataset(file, HDF5::transformationKindPath(d), kind);
        }
      }
  };

  /**
   * @brief HDF5 grid function loader for local grid functions backed by
   *        `Math::Vector<Scalar>`.
   *
   * Reads a standalone grid function field file (one `.h5` per grid
   * function) and populates the DOF vector. The file is expected to
   * contain `/GridFunction/Meta/{Size,Dimension}` and
   * `/GridFunction/Values/Data`.
   *
   * The grid function must already be attached to a finite element space
   * of matching size and dimension before loading.
   *
   * Stream-based loading is not supported; only file-path loading is available.
   *
   * @tparam FES     Finite element space type.
   * @tparam Scalar  Scalar element type of the data vector.
   *
   * ## Usage Example
   * ```cpp
   * P1 Vh(mesh);
   * GridFunction u(Vh);
   * IO::GridFunctionLoader<IO::FileFormat::HDF5, P1, Math::Vector<Real>> loader(u);
   * loader.load("field.h5");
   * ```
   *
   * @see GridFunctionPrinter<FileFormat::HDF5, FES, Math::Vector<Scalar>>
   */
  template <class FES, class Scalar>
  class GridFunctionLoader<FileFormat::HDF5, FES, Math::Vector<Scalar>>
    : public GridFunctionLoaderBase<FES, Math::Vector<Scalar>>
  {
    public:
      /// @brief Data storage type.
      using DataType = Math::Vector<Scalar>;

      /// @brief Grid function type being loaded.
      using ObjectType = Variational::GridFunction<FES, DataType>;

      /// @brief Parent loader class type.
      using Parent = GridFunctionLoaderBase<FES, DataType>;

      /**
       * @brief Constructs a loader for the given grid function.
       * @param[in,out] gf  Grid function to be populated with loaded data.
       */
      explicit
      GridFunctionLoader(ObjectType& gf)
        : Parent(gf)
      {}

      /**
       * @brief Stream-based loading is not supported for HDF5.
       *
       * Always raises an exception. Use the file-path overload instead.
       */
      void load(std::istream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 GridFunction loading is file-path based."
          << Alert::Raise;
      }

      /**
       * @brief Loads a grid function from a standalone HDF5 field file.
       * @param[in] filename  Path to the HDF5 field file.
       *
       * Reads the raw DOF vector from `/GridFunction/Values/Data` and
       * validates size/dimension against the attached finite element space.
       */
      void load(const boost::filesystem::path& filename) override
      {
        auto& gf = this->getObject();
        const auto file = HDF5::File(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!file)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to open HDF5 GridFunction file: " << filename
            << Alert::Raise;
        }

        const auto values = HDF5::readVectorDataset<HDF5::F64>(file.get(), HDF5::Path::GridFunctionValuesData);
        const auto dofCount = static_cast<size_t>(
            HDF5::readScalarDataset<HDF5::U64>(file.get(), HDF5::Path::GridFunctionMetaSize));
        const auto vectorDim = static_cast<size_t>(
            HDF5::readScalarDataset<HDF5::U64>(file.get(), HDF5::Path::GridFunctionMetaDimension));

        if (values.size() != static_cast<size_t>(gf.getData().size()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid GridFunction data size."
            << Alert::Raise;
        }
        if (dofCount != gf.getSize())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "GridFunction size mismatch."
            << Alert::Raise;
        }
        if (vectorDim != gf.getDimension())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "GridFunction dimension mismatch."
            << Alert::Raise;
        }

        auto& data = gf.getData();
        for (size_t i = 0; i < values.size(); ++i)
          data[i] = static_cast<typename std::remove_reference_t<decltype(data[0])>>(values[i]);
      }
  };

  /**
   * @brief HDF5 grid function printer for local grid functions backed by
   *        `Math::Vector<Scalar>`.
   *
   * Serializes a grid function to a standalone HDF5 field file (one file per
   * grid function, no mesh or FES embedding). The output layout is:
   *
   * ```
   * /GridFunction/Meta/Size        — number of DOFs
   * /GridFunction/Meta/Dimension   — vector dimension
   * /GridFunction/Values/Data      — raw DOF vector
   * ```
   *
   * This printer is strictly for canonical Rodin grid function persistence.
   * It writes the raw DOF vector only. The XDMF visualization path
   * (IO::XDMF) handles evaluation at vertices/cells separately using
   * `HDF5::writeXDMFNodeAttribute` / `HDF5::writeXDMFCellAttribute`.
   *
   * Stream-based printing is not supported; only file-path printing is available.
   *
   * @tparam FES     Finite element space type.
   * @tparam Scalar  Scalar element type of the data vector.
   *
   * ## Usage Example
   * ```cpp
   * const GridFunction& u = ...;
   * IO::GridFunctionPrinter<IO::FileFormat::HDF5, P1, Math::Vector<Real>>(u)
   *   .print("field.h5");
   * ```
   *
   * @see GridFunctionLoader<FileFormat::HDF5, FES, Math::Vector<Scalar>>
   */
  template <class FES, class Scalar>
  class GridFunctionPrinter<FileFormat::HDF5, FES, Math::Vector<Scalar>> final
    : public GridFunctionPrinterBase<FileFormat::HDF5, FES, Math::Vector<Scalar>>
  {
    public:
      /// @brief Data storage type.
      using DataType = Math::Vector<Scalar>;

      /// @brief Grid function type being printed.
      using ObjectType = Variational::GridFunction<FES, DataType>;

      /// @brief Parent printer class type.
      using Parent = GridFunctionPrinterBase<FileFormat::HDF5, FES, DataType>;

      /**
       * @brief Constructs a printer for the given grid function.
       * @param[in] gf  Grid function to serialize.
       */
      explicit
      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      /**
       * @brief Stream-based printing is not supported for HDF5.
       *
       * Always raises an exception. Use the file-path overload instead.
       */
      void print(std::ostream&) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "HDF5 GridFunction printing is file-path based."
          << Alert::Raise;
      }

      /**
       * @brief Writes the grid function to a standalone HDF5 field file.
       * @param[in] filename  Output HDF5 file path.
       *
       * Writes the raw DOF vector and metadata.
       */
      void print(const boost::filesystem::path& filename) override
      {
        const auto& gf = this->getObject();

        const auto file = HDF5::File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
        if (!file)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to create HDF5 GridFunction file: " << filename
            << Alert::Raise;
        }

        {
          const auto group = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::GridFunction, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!group)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /GridFunction group."
              << Alert::Raise;
          }
        }
        {
          const auto group = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::GridFunctionMeta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!group)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /GridFunction/Meta group."
              << Alert::Raise;
          }
        }
        {
          const auto group = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::GridFunctionValues, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
          if (!group)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to create /GridFunction/Values group."
              << Alert::Raise;
          }
        }

        const auto& data = gf.getData();
        std::vector<HDF5::F64> values(static_cast<size_t>(data.size()));
        for (size_t i = 0; i < values.size(); ++i)
          values[i] = static_cast<HDF5::F64>(data[i]);

        HDF5::writeVectorDataset(file.get(), HDF5::Path::GridFunctionValuesData, values);
        HDF5::writeScalarDataset(file.get(), HDF5::Path::GridFunctionMetaSize, static_cast<HDF5::U64>(gf.getSize()));
        HDF5::writeScalarDataset(file.get(), HDF5::Path::GridFunctionMetaDimension, static_cast<HDF5::U64>(gf.getDimension()));
      }
  };
}

#endif
