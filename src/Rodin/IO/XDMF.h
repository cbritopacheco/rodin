/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file XDMF.h
 * @brief XDMF 3 domain writer for temporal visualization output.
 *
 * This file defines the IO::XDMF class, which produces XDMF 3 XML documents
 * that reference HDF5 heavy-data files for mesh topology, geometry, and grid
 * function attribute datasets.
 *
 * The XDMF writer follows a split-file design:
 * - One HDF5 file per mesh snapshot
 * - One HDF5 file per grid function attribute per snapshot
 * - One XDMF XML file referencing all HDF5 files
 *
 * In distributed (MPI) mode, each rank writes its own rank-specific HDF5
 * visualization files using `{rank}` in the file patterns. Only the root
 * rank (configurable, default 0) writes the master `.xdmf` XML file that
 * references all rank-local piece files as a Spatial collection. For
 * transient output, the master XDMF uses a Temporal collection of Spatial
 * collections.
 *
 * ## Typical Workflow
 *
 * ```cpp
 * IO::XDMF xdmf("output/simulation");
 *
 * xdmf.setMesh(mesh)
 *      .add("velocity", u, IO::XDMF::Center::Node)
 *      .add("pressure", p, IO::XDMF::Center::Node);
 *
 * for (size_t step = 0; step < nSteps; ++step)
 * {
 *   // ... solve ...
 *   xdmf.write(time);
 * }
 *
 * xdmf.close();  // writes the XDMF XML; also called by destructor
 * ```
 *
 * @see IO::HDF5, IO::MeshPrinter<FileFormat::HDF5, Context::Local>,
 *      IO::GridFunctionPrinter<FileFormat::HDF5, FES, Data>
 * @see <a href="https://www.xdmf.org/index.php/XDMF_Model_and_Format">
 *      XDMF Model and Format Specification</a>
 */
#ifndef RODIN_IO_XDMF_H
#define RODIN_IO_XDMF_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

#include <boost/filesystem/path.hpp>

#include "Rodin/IO/HDF5.h"
#include "Rodin/Types.h"
#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"

#ifdef RODIN_USE_MPI
#include <mpi.h>
#endif

#include "ForwardDecls.h"

namespace Rodin::IO
{
  /**
   * @brief XDMF domain writer for static or transient visualization output.
   *
   * This class builds an XDMF 3 document referencing HDF5 heavy data.
   *
   * Conceptually:
   * - one XDMF object corresponds to one XDMF Domain,
   * - the domain contains one or more named grids,
   * - each grid owns one observed mesh and zero or more observed attributes,
   * - each call to write(...) appends one temporal snapshot for every configured grid.
   *
   * The writer is RAII-based:
   * - write(...) exports one new snapshot,
   * - close() finalizes the XML document,
   * - destruction closes automatically if necessary.
   *
   * Heavy data is stored in HDF5 files.
   * The XML document only references those files.
   */
  class XDMF
  {
    public:
      /**
       * @brief XDMF topology type identifiers.
       *
       * These values correspond to the XDMF specification's topology type
       * codes used in mixed-topology streams. Each cell in the flat topology
       * array is prefixed by its type id.
       *
       * @see <a href="https://www.xdmf.org/index.php/XDMF_Model_and_Format">
       *      XDMF Specification</a>
       */
      enum class Topology
      {
        POLYVERTEX      = 1,   ///< Single vertex (0D).
        POLYLINE        = 2,   ///< Line segment (1D).
        POLYGON         = 3,   ///< General polygon (2D).
        TRIANGLE        = 4,   ///< Triangle (2D, 3 vertices).
        QUADRILATERAL   = 5,   ///< Quadrilateral (2D, 4 vertices).
        TETRAHEDRON     = 6,   ///< Tetrahedron (3D, 4 vertices).
        PYRAMID         = 7,   ///< Pyramid (3D, 5 vertices).
        WEDGE           = 8,   ///< Wedge / triangular prism (3D, 6 vertices).
        HEXAHEDRON      = 9,   ///< Hexahedron (3D, 8 vertices).

        POLYHEDRON      = 16,  ///< General polyhedron (3D).

        EDGE_3          = 34,  ///< Quadratic edge (3 nodes).
        QUADRILATERAL_9 = 35,  ///< Biquadratic quadrilateral (9 nodes).
        TRIANGLE_6      = 36,  ///< Quadratic triangle (6 nodes).
        QUADRILATERAL_8 = 37,  ///< Serendipity quadrilateral (8 nodes).
        TETRAHEDRON_10  = 38,  ///< Quadratic tetrahedron (10 nodes).
        PYRAMID_13      = 39,  ///< Quadratic pyramid (13 nodes).
        WEDGE_15        = 40,  ///< Quadratic wedge (15 nodes).
        WEDGE_18        = 41,  ///< Biquadratic wedge (18 nodes).

        HEXAHEDRON_20   = 48,  ///< Serendipity hexahedron (20 nodes).
        HEXAHEDRON_24   = 49,  ///< Biquadratic hexahedron (24 nodes).
        HEXAHEDRON_27   = 50   ///< Triquadratic hexahedron (27 nodes).
      };

      /**
       * @brief XML element names used in the XDMF 3 document format.
       */
      struct Keyword
      {
        static constexpr const char* Xdmf      = "Xdmf";       ///< Root element.
        static constexpr const char* Domain    = "Domain";      ///< Domain container.
        static constexpr const char* Grid      = "Grid";        ///< Grid element (uniform or collection).
        static constexpr const char* Topology  = "Topology";    ///< Topology element.
        static constexpr const char* Geometry  = "Geometry";    ///< Geometry element.
        static constexpr const char* Attribute = "Attribute";   ///< Attribute element (field data).
        static constexpr const char* DataItem  = "DataItem";    ///< DataItem element (HDF5 reference).
      };

      /**
       * @brief Maps a Rodin polytope type to the corresponding XDMF topology type.
       * @param[in] geometry  Rodin polytope geometry type.
       * @returns The XDMF Topology value, or empty if the type is unsupported.
       */
      static inline
      Optional<Topology> getTopology(Geometry::Polytope::Type geometry)
      {
        using PT = Geometry::Polytope::Type;

        switch (geometry)
        {
          case PT::Point:         return Topology::POLYVERTEX;
          case PT::Segment:       return Topology::POLYLINE;
          case PT::Triangle:      return Topology::TRIANGLE;
          case PT::Quadrilateral: return Topology::QUADRILATERAL;
          case PT::Tetrahedron:   return Topology::TETRAHEDRON;
          case PT::Wedge:         return Topology::WEDGE;
          case PT::Hexahedron:    return Topology::HEXAHEDRON;
        }

        return {};
      }

      /**
       * @brief Maps an XDMF topology type back to a Rodin polytope type.
       * @param[in] gt  XDMF topology type code.
       * @returns The Rodin Polytope::Type, or empty if the type is unsupported.
       */
      static inline
      Optional<Geometry::Polytope::Type> getGeometry(Topology gt)
      {
        using PT = Geometry::Polytope::Type;

        switch (gt)
        {
          case Topology::POLYVERTEX:      return PT::Point;
          case Topology::POLYLINE:        return PT::Segment;
          case Topology::TRIANGLE:        return PT::Triangle;
          case Topology::QUADRILATERAL:   return PT::Quadrilateral;
          case Topology::TETRAHEDRON:     return PT::Tetrahedron;
          case Topology::WEDGE:           return PT::Wedge;
          case Topology::HEXAHEDRON:      return PT::Hexahedron;

          case Topology::PYRAMID:
          case Topology::EDGE_3:
          case Topology::QUADRILATERAL_9:
          case Topology::TRIANGLE_6:
          case Topology::QUADRILATERAL_8:
          case Topology::TETRAHEDRON_10:
          case Topology::PYRAMID_13:
          case Topology::WEDGE_15:
          case Topology::WEDGE_18:
          case Topology::HEXAHEDRON_20:
          case Topology::HEXAHEDRON_24:
          case Topology::HEXAHEDRON_27:
          case Topology::POLYHEDRON:
          case Topology::POLYGON:
            return {};
        }

        return {};
      }

      /**
       * @brief Mesh export policy.
       *
       * Static:
       * - the mesh is exported once and reused by all later snapshots of the same grid.
       *
       * Transient:
       * - the mesh is exported at every snapshot,
       * - intended for moving meshes or remeshing workflows.
       */
      enum class MeshPolicy
      {
        Static,    ///< Mesh exported once, reused by all snapshots of the grid.
        Transient  ///< Mesh exported at every snapshot (for moving meshes).
      };

      /**
       * @brief Attribute location in the XDMF sense.
       *
       * Node:
       * - export as nodal visualization data.
       *
       * Cell:
       * - export as cell-centered data.
       */
      enum class Center
      {
        Node,   ///< Vertex-centered (nodal) visualization data.
        Cell    ///< Cell-centered visualization data.
      };

      /**
       * @brief File naming patterns used by the writer.
       *
       * Supported placeholders:
       * - {stem}  : writer stem
       * - {grid}  : grid name
       * - {name}  : attribute name
       * - {index} : snapshot index with zero padding
       * - {rank}  : MPI rank (only meaningful in distributed mode)
       */
      struct FilePatterns
      {
        std::string xdmf          = "{stem}.xdmf";              ///< XDMF XML output filename pattern.
        std::string staticMesh    = "{stem}.{grid}.mesh.h5";     ///< Static mesh HDF5 filename pattern.
        std::string transientMesh = "{stem}.{grid}.mesh.{index}.h5";  ///< Transient mesh HDF5 filename pattern.
        std::string attribute     = "{stem}.{grid}.{name}.{index}.h5"; ///< Attribute field HDF5 filename pattern.
      };

      /**
       * @brief Per-grid export options.
       *
       * `patterns` overrides the global writer patterns only when set.
       */
      struct GridOptions
      {
        MeshPolicy meshPolicy = MeshPolicy::Static;  ///< Mesh export policy for this grid.
        Optional<FilePatterns> patterns;               ///< Per-grid file pattern overrides (empty = use writer defaults).
      };

      /**
       * @brief Handle to one named grid inside the XDMF domain.
       *
       * A Grid handle provides access to a specific grid within the XDMF
       * writer. Through this handle, the user can set the mesh, configure
       * export options, and register attributes. Handles are lightweight
       * (pointer + index) and are copyable.
       */
      class Grid
      {
        public:
          Grid(const Grid&) = default;
          Grid& operator=(const Grid&) = default;

          /**
           * @brief Returns the name of this grid.
           * @returns A string view of the grid name (empty for the default grid).
           */
          std::string_view getName() const noexcept;

          /**
           * @brief Sets the mesh observed by this grid.
           * @tparam MeshType  Concrete mesh type (must derive from MeshBase).
           * @param[in] mesh    Mesh to observe (must outlive the XDMF writer).
           * @param[in] policy  Export policy: Static or Transient.
           * @returns Reference to `*this` for method chaining.
           */
          template <class MeshType>
          Grid& setMesh(const MeshType& mesh, MeshPolicy policy = MeshPolicy::Static);

          /**
           * @brief Sets per-grid export options.
           * @param[in] options  Grid options including mesh policy and file patterns.
           * @returns Reference to `*this` for method chaining.
           */
          Grid& setOptions(const GridOptions& options);

          /**
           * @brief Returns the current grid options.
           * @returns Const reference to the GridOptions struct.
           */
          const GridOptions& getOptions() const noexcept;

          /**
           * @brief Registers a grid function attribute using its internal name.
           * @tparam GridFunctionType  Concrete grid function type.
           * @param[in] gf      Grid function to export (must have a name set).
           * @param[in] center  Data centering (Node or Cell).
           * @returns Reference to `*this` for method chaining.
           *
           * The grid function must have a name assigned via `setName()`.
           * Use the overload taking an explicit name string to override.
           */
          template <class GridFunctionType>
          Grid& add(const GridFunctionType& gf, Center center = Center::Node);

          /**
           * @brief Registers a grid function attribute with an explicit name.
           * @tparam GridFunctionType  Concrete grid function type.
           * @param[in] name    Attribute name for the XDMF output.
           * @param[in] gf      Grid function to export.
           * @param[in] center  Data centering (Node or Cell).
           * @returns Reference to `*this` for method chaining.
           */
          template <class GridFunctionType>
          Grid& add(
              const std::string& name,
              const GridFunctionType& gf,
              Center center = Center::Node);

          /**
           * @brief Tests whether a mesh has been set on this grid.
           * @returns `true` if setMesh() has been called.
           */
          bool hasMesh() const noexcept;

          /**
           * @brief Returns the number of registered attributes.
           * @returns Attribute count.
           */
          size_t getAttributeCount() const noexcept;

          /**
           * @brief Resets the grid to its initial state.
           *
           * Clears the mesh, options, attributes, and snapshot history.
           * @returns Reference to `*this` for method chaining.
           */
          Grid& reset();

          /**
           * @brief Removes all registered attributes from this grid.
           * @returns Reference to `*this` for method chaining.
           */
          Grid& clear();

        private:
          friend class XDMF;

          Grid(XDMF& owner, size_t index) noexcept;

          XDMF* m_owner = nullptr;
          size_t m_index = 0;
      };

      /**
       * @brief Constructs an XDMF writer with the given stem path (serial mode).
       * @param[in] stem  Base path for output files. The directory component
       *                  determines where files are written; the filename
       *                  component is used as `{stem}` in file patterns.
       *
       * ## Example
       * ```cpp
       * IO::XDMF xdmf("output/simulation");
       * // produces: output/simulation.xdmf, output/simulation.*.h5, ...
       * ```
       */
      explicit
      XDMF(const boost::filesystem::path& stem);

      /**
       * @brief Constructs an XDMF writer in distributed (MPI) mode.
       *
       * In distributed mode:
       * - Each rank writes its own rank-specific HDF5 visualization files
       *   (using `{rank}` in the file patterns).
       * - Only the root rank writes the master `.xdmf` XML file that
       *   references all per-rank piece files as a Spatial collection.
       * - For transient output, the master XDMF contains a Temporal
       *   collection of Spatial collections.
       * - write(), flush(), and close() are collective operations that
       *   must be called by all ranks.
       *
       * @param[in] comm      MPI communicator.
       * @param[in] stem      Base path for output files.
       * @param[in] rootRank  Rank that writes the master XDMF XML (default 0).
       *
       * ## Example
       * ```cpp
       * IO::XDMF xdmf(MPI_COMM_WORLD, "output/Poisson");
       * xdmf.setMesh(mpiMesh);
       * xdmf.add("u", u, IO::XDMF::Center::Node);
       * xdmf.write(0.0);
       * xdmf.close();
       * ```
       */
#ifdef RODIN_USE_MPI
      XDMF(MPI_Comm comm,
           const boost::filesystem::path& stem,
           size_t rootRank = 0);
#endif

      XDMF(const XDMF&) = delete;        ///< Non-copyable.
      XDMF& operator=(const XDMF&) = delete;  ///< Non-copyable.
      XDMF(XDMF&&) = default;           ///< Move constructible.
      XDMF& operator=(XDMF&&) = default; ///< Move assignable.

      ~XDMF() = default;

      /**
       * @brief Returns the stem path set at construction.
       * @returns Const reference to the stem path.
       */
      const boost::filesystem::path& getStem() const noexcept;

      /**
       * @brief Sets the global file naming patterns.
       * @param[in] patterns  File patterns with `{stem}`, `{grid}`,
       *                      `{name}`, and `{index}` placeholders.
       * @returns Reference to `*this` for method chaining.
       */
      XDMF& setFilePatterns(const FilePatterns& patterns);

      /**
       * @brief Returns the current global file patterns.
       * @returns Const reference to the FilePatterns struct.
       */
      const FilePatterns& getFilePatterns() const noexcept;

      /**
       * @brief Sets the zero-padding width for `{index}` expansion.
       * @param[in] digits  Number of digits for zero-padded indices (default: 6).
       * @returns Reference to `*this` for method chaining.
       */
      XDMF& setPadding(size_t digits);

      /**
       * @brief Returns the current index padding width.
       * @returns Number of zero-padding digits.
       */
      size_t getPadding() const noexcept;

      /**
       * @brief Returns a handle to the default (unnamed) grid.
       *
       * Creates the default grid if it does not yet exist.
       * @returns Grid handle for the default grid.
       */
      Grid grid();

      /**
       * @brief Returns a handle to a named grid.
       *
       * Creates the grid if it does not yet exist. Subsequent calls with the
       * same name return a handle to the same grid.
       *
       * @param[in] name  Grid name (used in file patterns as `{grid}`).
       * @returns Grid handle for the named grid.
       */
      Grid grid(const std::string& name);

      /**
       * @brief Sets the mesh on the default grid.
       * @tparam MeshType  Concrete mesh type.
       * @param[in] mesh    Mesh to observe.
       * @param[in] policy  Export policy (Static or Transient).
       * @returns Reference to `*this` for method chaining.
       */
      template <class MeshType>
      XDMF& setMesh(const MeshType& mesh, MeshPolicy policy = MeshPolicy::Static);

      /**
       * @brief Sets export options on the default grid.
       * @param[in] options  Grid options.
       * @returns Reference to `*this` for method chaining.
       */
      XDMF& setOptions(const GridOptions& options);

      /**
       * @brief Adds a named grid function attribute to the default grid.
       * @tparam GridFunctionType  Concrete grid function type.
       * @param[in] gf      Grid function to export (must have a name).
       * @param[in] center  Data centering (Node or Cell).
       * @returns Reference to `*this` for method chaining.
       */
      template <class GridFunctionType>
      XDMF& add(const GridFunctionType& gf, Center center = Center::Node);

      /**
       * @brief Adds a grid function attribute with an explicit name to the
       *        default grid.
       * @tparam GridFunctionType  Concrete grid function type.
       * @param[in] name    Attribute name.
       * @param[in] gf      Grid function to export.
       * @param[in] center  Data centering (Node or Cell).
       * @returns Reference to `*this` for method chaining.
       */
      template <class GridFunctionType>
      XDMF& add(
          const std::string& name,
          const GridFunctionType& gf,
          Center center = Center::Node);

      /**
       * @brief Writes one temporal snapshot for all configured grids.
       *
       * Uses the current snapshot count as the time value.
       * @returns Reference to `*this` for method chaining.
       */
      XDMF& write();

      /**
       * @brief Writes one temporal snapshot at the given time value.
       * @param[in] time  Physical time associated with this snapshot.
       * @returns Reference to `*this` for method chaining.
       *
       * For each grid, exports the mesh (if needed) and all registered
       * attributes to HDF5 files, and records the snapshot for the final
       * XDMF XML output.
       */
      XDMF& write(Real time);

      /**
       * @brief Finalizes the XDMF XML document and writes it to disk.
       *
       * Generates the complete XDMF 3 XML referencing all recorded
       * snapshots and their HDF5 data files. This method is idempotent;
       * subsequent calls are no-ops. The destructor calls close()
       * automatically.
       */
      void close();

      /**
       * @brief Tests whether the writer has been closed.
       * @returns `true` if close() has been called.
       */
      bool isClosed() const noexcept;

      /**
       * @brief Returns the total number of snapshots written so far.
       * @returns Snapshot count.
       */
      size_t getSnapshotCount() const noexcept;

      /**
       * @brief Returns the number of grids in this domain.
       * @returns Grid count.
       */
      size_t getGridCount() const noexcept;

      /**
       * @brief Tests whether the writer is in distributed (MPI) mode.
       * @returns `true` if constructed with rank/numRanks.
       */
      bool isDistributed() const noexcept;

      /**
       * @brief Returns this process's rank.
       * @returns Rank (0 in serial mode).
       */
      size_t getRank() const noexcept;

      /**
       * @brief Returns the total number of ranks.
       * @returns Number of ranks (1 in serial mode).
       */
      size_t getNumRanks() const noexcept;

      /**
       * @brief Returns the root rank that writes the master XDMF XML.
       * @returns Root rank (0 in serial mode).
       */
      size_t getRootRank() const noexcept;

      void flush() const;

    private:
      /// @brief Internal record for one registered attribute.
      struct AttributeRecord
      {
        std::string name;                           ///< Attribute display name.
        Center center = Center::Node;               ///< Data centering.
        size_t dimension = 1;                       ///< Vector dimension of the attribute.
        std::function<void(const boost::filesystem::path&, Center)> write;  ///< Save callback.
      };

      /// @brief Internal record for one temporal snapshot.
      struct SnapshotRecord
      {
        struct AttributeRecord
        {
          std::string name;
          Center center = Center::Node;
          size_t dimension = 1;
          boost::filesystem::path file;
        };

        struct MeshAttributeRecord
        {
          std::string name;
          Center center = Center::Cell;
          size_t topologicalDimension = 0;
        };

        /// @brief Per-rank mesh metadata cached during write().
        struct PieceMeta
        {
          std::uint64_t vertexCount = 0;
          std::uint64_t cellCount = 0;
          std::uint64_t meshDimension = 0;
          std::uint64_t spaceDimension = 0;
          std::uint64_t topologySize = 0;
        };

        Real time = 0;
        boost::filesystem::path meshFile;
        std::vector<AttributeRecord> attributes;
        std::vector<MeshAttributeRecord> meshAttributes;

        size_t vertexCount = 0;
        size_t cellCount = 0;
        size_t meshDimension = 0;
        size_t spaceDimension = 0;
        size_t topologySize = 0;

        /// @brief Per-rank metadata (size 1 in serial, numRanks on root in distributed).
        std::vector<PieceMeta> pieces;
      };

      /// @brief Internal record for one named grid.
      struct GridRecord
      {
        std::string name;                             ///< Grid name.
        const Geometry::MeshBase* mesh = nullptr;     ///< Effective mesh for visualization (shard for MPI).
        const Geometry::MeshBase* sourceMesh = nullptr; ///< Original mesh for identity checks.
        GridOptions options;                          ///< Per-grid export options.
        bool staticMeshWritten = false;               ///< Whether the static mesh has been exported.
        boost::filesystem::path staticMeshFile;       ///< Path to the static mesh file.
        std::vector<AttributeRecord> attributes;      ///< Registered attributes.
        std::vector<SnapshotRecord> snapshots;        ///< Recorded snapshots.
      };

      boost::filesystem::path m_stem;     ///< Output stem path.
      FilePatterns m_patterns;            ///< Global file naming patterns.
      size_t m_padding = 6;               ///< Zero-padding width for index expansion.
      bool m_closed = false;              ///< Whether close() has been called.
      size_t m_snapshotCount = 0;         ///< Total snapshots written.
      std::vector<GridRecord> m_grids;    ///< All grids in this domain.

      // --- MPI distributed mode -----------------------------------------------
      bool m_distributed = false;         ///< Whether the writer is in distributed mode.
      size_t m_rank = 0;                  ///< This process's rank.
      size_t m_numRanks = 1;              ///< Total number of ranks.
      size_t m_rootRank = 0;              ///< Root rank that writes master XDMF XML.
#ifdef RODIN_USE_MPI
      MPI_Comm m_comm = MPI_COMM_NULL;    ///< Non-owning MPI communicator.
#endif

      /// @brief Writes a single Uniform grid XML element.
      void writeUniformGrid(
          std::ostream& os,
          const std::string& gridName,
          const SnapshotRecord& snap,
          size_t baseIndent) const;

      /// @brief Gathers per-rank mesh metadata into snap.pieces.
      void gatherPieceMeta(SnapshotRecord& snap) const;
  };

  // ---- template method implementations ------------------------------------

  template <class MeshType>
  XDMF::Grid& XDMF::Grid::setMesh(const MeshType& mesh, MeshPolicy policy)
  {
    auto& gr = m_owner->m_grids[m_index];
    // For distributed (MPI) meshes, store the shard for visualization
    // since the writeXDMF helpers require a local mesh. The original mesh
    // pointer is kept in sourceMesh for identity checks in add().
    if constexpr (requires { mesh.getShard(); })
      gr.mesh = &mesh.getShard();
    else
      gr.mesh = &mesh;
    gr.sourceMesh = &mesh;
    gr.options.meshPolicy = policy;
    gr.staticMeshWritten = false;
    gr.staticMeshFile.clear();
    return *this;
  }

  template <class GridFunctionType>
  XDMF::Grid& XDMF::Grid::add(const GridFunctionType& gf, Center center)
  {
    const auto name = gf.getName();
    if (!name)
    {
      Alert::Exception()
        << "Grid function has no name. Use the overload that takes an explicit name."
        << Alert::Raise;
    }
    return add(std::string(name->data(), name->size()), gf, center);
  }

  template <class GridFunctionType>
  XDMF::Grid& XDMF::Grid::add(
      const std::string& name,
      const GridFunctionType& gf,
      Center center)
  {
    auto& gr = m_owner->m_grids[m_index];

    if (gr.sourceMesh && gr.sourceMesh != &gf.getFiniteElementSpace().getMesh())
    {
      Alert::Exception()
        << "Attribute mesh does not match the grid mesh."
        << Alert::Raise;
    }

    for (const auto& attr : gr.attributes)
    {
      if (attr.name == name)
      {
        Alert::Exception()
          << "Duplicate XDMF attribute name \"" << name
          << "\" in grid \"" << gr.name << "\"."
          << Alert::Raise;
      }
    }

    AttributeRecord rec;
    rec.name = name;
    rec.center = center;
    rec.dimension = gf.getDimension();
    // Capture the effective visualization mesh (shard for MPI, local mesh
    // for serial). The attribute write helpers iterate over this mesh's
    // vertices/cells, not over gf.getFiniteElementSpace().getMesh(), which
    // may be a distributed MPI mesh with global counts.
    const auto* visMesh = gr.mesh;
    rec.write = [&gf, visMesh](const boost::filesystem::path& path, Center center)
    {
      switch (center)
      {
        case Center::Node:
          if (visMesh)
            HDF5::writeXDMFNodeAttribute(gf, *visMesh, path);
          else
            HDF5::writeXDMFNodeAttribute(gf, path);
          return;
        case Center::Cell:
          if (visMesh)
            HDF5::writeXDMFCellAttribute(gf, *visMesh, path);
          else
            HDF5::writeXDMFCellAttribute(gf, path);
          return;
      }
    };
    gr.attributes.push_back(std::move(rec));
    return *this;
  }

  template <class MeshType>
  XDMF& XDMF::setMesh(const MeshType& mesh, MeshPolicy policy)
  {
    grid().setMesh(mesh, policy);
    return *this;
  }

  template <class GridFunctionType>
  XDMF& XDMF::add(const GridFunctionType& gf, Center center)
  {
    grid().add(gf, center);
    return *this;
  }

  template <class GridFunctionType>
  XDMF& XDMF::add(
      const std::string& name,
      const GridFunctionType& gf,
      Center center)
  {
    grid().add(name, gf, center);
    return *this;
  }
}

#endif
