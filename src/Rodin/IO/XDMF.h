/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_XDMF_H
#define RODIN_IO_XDMF_H

#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <boost/filesystem/path.hpp>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::IO
{
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

      template <class MeshType>
      XDMF(const MeshType& mesh, const boost::filesystem::path& h5MeshFile)
        : m_h5MeshFile(h5MeshFile),
          m_vertexCount(mesh.getVertexCount()),
          m_spaceDimension(mesh.getSpaceDimension()),
          m_cellCount(mesh.getCellCount()),
          m_xdmfTopologySize(0)
      {
        for (auto it = mesh.getCell(); !it.end(); ++it)
          m_xdmfTopologySize += xdmfCellWordCount(it->getGeometry(), it->getVertices().size());
      }

      template <class GridFunctionType>
      XDMF& add(
          const std::string& name,
          const GridFunctionType& gf,
          const boost::filesystem::path& h5File)
      {
        Attribute attr;
        attr.name = name;
        attr.h5File = h5File;
        attr.datasetPath = "/GridFunction/Values/Data";
        attr.dofCount = gf.getSize();
        attr.components = gf.getDimension();
        m_attributes.push_back(std::move(attr));
        return *this;
      }

      void save(const boost::filesystem::path& filename) const
      {
        std::ofstream os(filename.c_str());
        if (!os)
          Alert::Exception() << "Failed to open XDMF output file." << Alert::Raise;

        os << "<?xml version=\"1.0\" ?>\n";
        os << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        os << "<" << Keyword::Xdmf << " Version=\"3.0\">\n";
        os << "  <" << Keyword::Domain << ">\n";
        os << "    <" << Keyword::Grid << " Name=\"RodinMesh\" GridType=\"Uniform\">\n";

        this->writeTopology(os);
        this->writeGeometry(os);
        this->writeAttributes(os);

        os << "    </" << Keyword::Grid << ">\n";
        os << "  </" << Keyword::Domain << ">\n";
        os << "</" << Keyword::Xdmf << ">\n";
      }

    private:
      struct Attribute
      {
        std::string name;
        boost::filesystem::path h5File;
        std::string datasetPath;
        size_t dofCount;
        size_t components;
      };

      const char* geometryType() const
      {
        return m_spaceDimension >= 3 ? "XYZ" : "XY";
      }

      const char* attributeType(size_t components) const
      {
        return components > 1 ? "Vector" : "Scalar";
      }

      void writeTopology(std::ostream& os) const
      {
        os << "      <" << Keyword::Topology
           << " TopologyType=\"Mixed\" NumberOfElements=\"" << m_cellCount << "\">\n";
        os << "        <" << Keyword::DataItem
           << " Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << m_xdmfTopologySize << "\">"
           << m_h5MeshFile.string() << ":/Mesh/XDMF/Topology</" << Keyword::DataItem << ">\n";
        os << "      </" << Keyword::Topology << ">\n";
      }

      void writeGeometry(std::ostream& os) const
      {
        os << "      <" << Keyword::Geometry << " GeometryType=\"" << geometryType() << "\">\n";
        os << "        <" << Keyword::DataItem
           << " Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\""
           << m_vertexCount << " " << m_spaceDimension << "\">"
           << m_h5MeshFile.string() << ":/Mesh/Geometry/Vertices</" << Keyword::DataItem << ">\n";
        os << "      </" << Keyword::Geometry << ">\n";
      }

      static size_t xdmfCellWordCount(Geometry::Polytope::Type g, size_t vertexCount)
      {
        switch (g)
        {
          case Geometry::Polytope::Type::Segment:
            // In XDMF Mixed topology, polyline entries carry the vertex-count
            // field explicitly: [2, count, v0, ..., v(count-1)].
            return 2 + vertexCount;
          default:
            return 1 + vertexCount;
        }
      }

      void writeAttributes(std::ostream& os) const
      {
        for (const auto& attr : m_attributes)
        {
          os << "      <" << Keyword::Attribute << " Name=\"" << attr.name
             << "\" AttributeType=\"" << attributeType(attr.components) << "\" Center=\"Node\">\n";
          os << "        <" << Keyword::DataItem
             << " Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\""
             << attr.dofCount;
          if (attr.components > 1)
            os << " " << attr.components;
          os << "\">" << attr.h5File.string() << ":" << attr.datasetPath
             << "</" << Keyword::DataItem << ">\n";
          os << "      </" << Keyword::Attribute << ">\n";
        }
      }

      boost::filesystem::path m_h5MeshFile;
      size_t m_vertexCount;
      size_t m_spaceDimension;
      size_t m_cellCount;
      size_t m_xdmfTopologySize;
      std::vector<Attribute> m_attributes;
  };
}

#endif
