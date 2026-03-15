/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_XDMF_H
#define RODIN_IO_XDMF_H

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>
#include <functional>
#include <boost/filesystem/path.hpp>

#include "Rodin/Types.h"
#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"

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

      enum class Topology
      {
        POLYVERTEX      = 1,
        POLYLINE        = 2,
        POLYGON         = 3,
        TRIANGLE        = 4,
        QUADRILATERAL   = 5,
        TETRAHEDRON     = 6,
        PYRAMID         = 7,
        WEDGE           = 8,
        HEXAHEDRON      = 9,

        POLYHEDRON      = 16,

        EDGE_3          = 34,
        QUADRILATERAL_9 = 35,
        TRIANGLE_6      = 36,
        QUADRILATERAL_8 = 37,
        TETRAHEDRON_10  = 38,
        PYRAMID_13      = 39,
        WEDGE_15        = 40,
        WEDGE_18        = 41,

        HEXAHEDRON_20   = 48,
        HEXAHEDRON_24   = 49,
        HEXAHEDRON_27   = 50
      };

      struct Keyword
      {
        static constexpr const char* Xdmf = "Xdmf";
        static constexpr const char* Domain = "Domain";
        static constexpr const char* Grid = "Grid";
        static constexpr const char* Topology = "Topology";
        static constexpr const char* Geometry = "Geometry";
        static constexpr const char* Attribute = "Attribute";
        static constexpr const char* DataItem = "DataItem";
      };

      inline
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

      inline
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
        Static,
        Transient
      };

      /**
       * @brief Attribute location in the XDMF sense.
       *
       * Node:
       * - export as nodal visualization data,
       * - by default this means evaluating/interpolating the field at mesh vertices.
       *
       * Cell:
       * - export as cell-centered data.
       */
      enum class Center
      {
        Node,
        Cell
      };

      /**
       * @brief File naming patterns used by the writer.
       *
       * Supported placeholders:
       * - {stem}  : writer stem
       * - {grid}  : grid name
       * - {name}  : attribute name
       * - {index} : snapshot index with zero padding
       *
       * Intended use:
       * - static mesh file pattern,
       * - transient mesh file pattern,
       * - attribute file pattern,
       * - final xdmf file pattern.
       */
      struct FilePatterns
      {
        std::string xdmf          = "{stem}.xdmf";
        std::string staticMesh    = "{stem}.{grid}.mesh.h5";
        std::string transientMesh = "{stem}.{grid}.mesh.{index}.h5";
        std::string attribute     = "{stem}.{grid}.{name}.{index}.h5";
      };

      /**
       * @brief Per-grid export options.
       *
       * These override the global writer defaults for one grid only.
       */
      struct GridOptions
      {
        MeshPolicy meshPolicy = MeshPolicy::Static;
        FilePatterns patterns = {};
      };

      /**
       * @brief Handle to one named grid inside the XDMF domain.
       *
       * A grid corresponds to one XDMF Grid family:
       * - it stores one observed mesh,
       * - it stores a set of observed attributes,
       * - it records temporal snapshots across successive write(...) calls.
       *
       * The handle is lightweight and non-owning.
       */
      class Grid
      {
        public:
          Grid(const Grid&) = default;
          Grid& operator=(const Grid&) = default;

          /**
           * @brief Gets the grid name.
           */
          std::string_view getName() const noexcept;

          /**
           * @brief Sets the mesh observed by this grid.
           * @param[in] mesh Mesh to export
           * @param[in] policy Mesh export policy
           * @returns Reference to this grid handle
           *
           * Behaviour:
           * - stores a non-owning reference to @p mesh,
           * - resets mesh export state for this grid,
           * - later write(...) calls will export this mesh according to @p policy.
           *
           * The mesh name used in XML defaults to:
           * - mesh.getName() if available,
           * - otherwise the grid name.
           */
          template <class MeshType>
          Grid& setMesh(const MeshType& mesh, MeshPolicy policy = MeshPolicy::Static);

          /**
           * @brief Sets per-grid export options.
           * @param[in] options Options to apply
           * @returns Reference to this grid handle
           *
           * Behaviour:
           * - replaces the current per-grid options,
           * - affects subsequent write(...) calls only.
           */
          Grid& setOptions(const GridOptions& options);

          /**
           * @brief Gets the current per-grid options.
           */
          const GridOptions& getOptions() const noexcept;

          /**
           * @brief Adds an attribute using the field's own name.
           * @param[in] gf Grid function to observe
           * @param[in] center Export center
           * @returns Reference to this grid handle
           *
           * Behaviour:
           * - the attribute name is taken from gf.getName(),
           * - throws if the field has no name,
           * - stores a non-owning reference to the field,
           * - validates that the field mesh is the same mesh as this grid.
           *
           * The field is not exported immediately.
           * Its current value is exported on each write(...).
           */
          template <class GridFunctionType>
          Grid& add(const GridFunctionType& gf, Center center = Center::Node);

          /**
           * @brief Adds an attribute with an explicit name.
           * @param[in] name Attribute name
           * @param[in] gf Grid function to observe
           * @param[in] center Export center
           * @returns Reference to this grid handle
           *
           * Behaviour:
           * - same as add(gf, center), but uses @p name instead of gf.getName().
           */
          template <class GridFunctionType>
          Grid& add(
              const std::string& name,
              const GridFunctionType& gf,
              Center center = Center::Node);

          /**
           * @brief Returns true if this grid already has an observed mesh.
           */
          bool hasMesh() const noexcept;

          /**
           * @brief Returns the number of registered attributes.
           */
          size_t getAttributeCount() const noexcept;

        private:
          friend class XDMF;

          Grid(XDMF& owner, size_t index) noexcept;

          XDMF* m_owner = nullptr;
          size_t m_index = 0;
      };

      /**
       * @brief Constructs an XDMF writer.
       * @param[in] stem Output stem
       *
       * The stem is the base used to build the XML and HDF5 filenames.
       */
      explicit
      XDMF(const boost::filesystem::path& stem);

      XDMF(const XDMF&) = delete;
      XDMF& operator=(const XDMF&) = delete;
      XDMF(XDMF&&) = default;
      XDMF& operator=(XDMF&&) = default;

      /**
       * @brief Destructor.
       *
       * Behaviour:
       * - if the writer is not closed yet, close() is called automatically.
       */
      ~XDMF();

      /**
       * @brief Gets the output stem.
       */
      const boost::filesystem::path& getStem() const noexcept;

      /**
       * @brief Sets global file naming patterns.
       * @param[in] patterns New patterns
       * @returns Reference to this writer
       *
       * These defaults are used by all grids unless a grid overrides them
       * through Grid::setOptions(...).
       */
      XDMF& setFilePatterns(const FilePatterns& patterns);

      /**
       * @brief Gets the global file naming patterns.
       */
      const FilePatterns& getFilePatterns() const noexcept;

      /**
       * @brief Sets the zero-padding width used for {index}.
       * @param[in] digits Padding width
       * @returns Reference to this writer
       */
      XDMF& setPadding(size_t digits);

      /**
       * @brief Gets the current zero-padding width.
       */
      size_t getPadding() const noexcept;

      /**
       * @brief Returns the default grid.
       *
       * This is the shorthand for the common single-mesh workflow.
       * The default grid name is implementation-defined but stable.
       */
      Grid grid();

      /**
       * @brief Returns the named grid, creating it if necessary.
       * @param[in] name Grid name
       * @returns Grid handle
       */
      Grid grid(const std::string& name);

      /**
       * @brief Convenience wrapper on the default grid.
       * @param[in] mesh Mesh to export
       * @param[in] policy Mesh policy
       * @returns Reference to this writer
       */
      template <class MeshType>
      XDMF& setMesh(const MeshType& mesh, MeshPolicy policy = MeshPolicy::Static);

      /**
       * @brief Convenience wrapper on the default grid.
       * @param[in] options Per-grid options for the default grid
       * @returns Reference to this writer
       */
      XDMF& setOptions(const GridOptions& options);

      /**
       * @brief Convenience wrapper on the default grid.
       * @param[in] gf Grid function to export
       * @param[in] center Export center
       * @returns Reference to this writer
       */
      template <class GridFunctionType>
      XDMF& add(const GridFunctionType& gf, Center center = Center::Node);

      /**
       * @brief Convenience wrapper on the default grid.
       * @param[in] name Attribute name
       * @param[in] gf Grid function to export
       * @param[in] center Export center
       * @returns Reference to this writer
       */
      template <class GridFunctionType>
      XDMF& add(
          const std::string& name,
          const GridFunctionType& gf,
          Center center = Center::Node);

      /**
       * @brief Writes one snapshot using the current snapshot index as time.
       * @returns Reference to this writer
       *
       * Behaviour:
       * - increments the snapshot index,
       * - for each configured grid:
       *   - verifies that a mesh is set,
       *   - exports the mesh according to its mesh policy,
       *   - exports all registered attributes,
       *   - records one temporal snapshot with time equal to the snapshot index.
       */
      XDMF& write();

      /**
       * @brief Writes one snapshot at the given physical time.
       * @param[in] time Snapshot time
       * @returns Reference to this writer
       *
       * Behaviour:
       * - increments the snapshot index,
       * - for each configured grid:
       *   - verifies that a mesh is set,
       *   - exports the mesh according to its mesh policy,
       *   - exports all registered attributes,
       *   - records one temporal snapshot with the given @p time.
       */
      XDMF& write(Real time);

      /**
       * @brief Finalizes the XML document and writes it to disk.
       *
       * The output filename is generated from the current xdmf file pattern.
       * After close(), further write(...) calls are invalid.
       */
      void close();

      /**
       * @brief Returns true if the writer has been closed.
       */
      bool isClosed() const noexcept;

      /**
       * @brief Returns the number of snapshots already written.
       */
      size_t getSnapshotCount() const noexcept;

      /**
       * @brief Returns the number of configured grids.
       */
      size_t getGridCount() const noexcept;

    private:
      struct AttributeRecord
      {
        std::string name;
        Center center = Center::Node;
        size_t dimension = 1;
        std::function<void(const boost::filesystem::path&)> save;
      };

      struct SnapshotRecord
      {
        Real time = 0;
        boost::filesystem::path meshFile;
        std::vector<boost::filesystem::path> attributeFiles;
      };

      struct GridRecord
      {
        std::string name;
        const Geometry::MeshBase* mesh = nullptr;
        GridOptions options;
        bool staticMeshWritten = false;
        boost::filesystem::path staticMeshFile;
        std::vector<AttributeRecord> attributes;
        std::vector<SnapshotRecord> snapshots;
      };

      boost::filesystem::path m_stem;
      FilePatterns m_patterns;
      size_t m_padding = 6;
      bool m_closed = false;
      size_t m_snapshotCount = 0;
      std::vector<GridRecord> m_grids;
  };

  // ---- template method implementations ------------------------------------

  template <class MeshType>
  XDMF::Grid& XDMF::Grid::setMesh(const MeshType& mesh, MeshPolicy policy)
  {
    auto& gr = m_owner->m_grids[m_index];
    gr.mesh = &mesh;
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
    if (gr.mesh && gr.mesh != &gf.getFiniteElementSpace().getMesh())
    {
      Alert::Exception()
        << "Attribute mesh does not match the grid mesh."
        << Alert::Raise;
    }

    AttributeRecord rec;
    rec.name = name;
    rec.center = center;
    rec.dimension = gf.getDimension();
    rec.save = [&gf](const boost::filesystem::path& path)
    {
      gf.save(path, FileFormat::HDF5);
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
