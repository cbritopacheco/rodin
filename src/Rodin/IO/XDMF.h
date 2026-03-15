/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_XDMF_H
#define RODIN_IO_XDMF_H

#include <cstddef>
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
        static constexpr const char* Xdmf      = "Xdmf";
        static constexpr const char* Domain    = "Domain";
        static constexpr const char* Grid      = "Grid";
        static constexpr const char* Topology  = "Topology";
        static constexpr const char* Geometry  = "Geometry";
        static constexpr const char* Attribute = "Attribute";
        static constexpr const char* DataItem  = "DataItem";
      };

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
        Static,
        Transient
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
       * `patterns` overrides the global writer patterns only when set.
       */
      struct GridOptions
      {
        MeshPolicy meshPolicy = MeshPolicy::Static;
        Optional<FilePatterns> patterns;
      };

      /**
       * @brief Handle to one named grid inside the XDMF domain.
       */
      class Grid
      {
        public:
          Grid(const Grid&) = default;
          Grid& operator=(const Grid&) = default;

          std::string_view getName() const noexcept;

          template <class MeshType>
          Grid& setMesh(const MeshType& mesh, MeshPolicy policy = MeshPolicy::Static);

          Grid& setOptions(const GridOptions& options);

          const GridOptions& getOptions() const noexcept;

          template <class GridFunctionType>
          Grid& add(const GridFunctionType& gf, Center center = Center::Node);

          template <class GridFunctionType>
          Grid& add(
              const std::string& name,
              const GridFunctionType& gf,
              Center center = Center::Node);

          bool hasMesh() const noexcept;

          size_t getAttributeCount() const noexcept;

        private:
          friend class XDMF;

          Grid(XDMF& owner, size_t index) noexcept;

          XDMF* m_owner = nullptr;
          size_t m_index = 0;
      };

      explicit
      XDMF(const boost::filesystem::path& stem);

      XDMF(const XDMF&) = delete;
      XDMF& operator=(const XDMF&) = delete;
      XDMF(XDMF&&) = default;
      XDMF& operator=(XDMF&&) = default;

      ~XDMF();

      const boost::filesystem::path& getStem() const noexcept;

      XDMF& setFilePatterns(const FilePatterns& patterns);

      const FilePatterns& getFilePatterns() const noexcept;

      XDMF& setPadding(size_t digits);

      size_t getPadding() const noexcept;

      Grid grid();

      Grid grid(const std::string& name);

      template <class MeshType>
      XDMF& setMesh(const MeshType& mesh, MeshPolicy policy = MeshPolicy::Static);

      XDMF& setOptions(const GridOptions& options);

      template <class GridFunctionType>
      XDMF& add(const GridFunctionType& gf, Center center = Center::Node);

      template <class GridFunctionType>
      XDMF& add(
          const std::string& name,
          const GridFunctionType& gf,
          Center center = Center::Node);

      XDMF& write();

      XDMF& write(Real time);

      void close();

      bool isClosed() const noexcept;

      size_t getSnapshotCount() const noexcept;

      size_t getGridCount() const noexcept;

    private:
      struct AttributeRecord
      {
        std::string name;
        Center center = Center::Node;
        size_t dimension = 1;
        std::function<void(const boost::filesystem::path&, Center)> write;
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
    rec.write = [&gf](const boost::filesystem::path& path, Center center)
    {
      IO::GridFunctionPrinter<
          IO::FileFormat::HDF5,
          typename std::remove_cvref_t<GridFunctionType>::FESType,
          typename std::remove_cvref_t<GridFunctionType>::DataType>(gf)
        .setXDMF(true)
        .setCenter(center == Center::Node ? IO::HDF5::Center::Node
                                          : IO::HDF5::Center::Cell)
        .print(path);
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
