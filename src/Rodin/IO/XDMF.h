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

      template <class MeshType>
      XDMF(const MeshType& mesh, const boost::filesystem::path& h5MeshFile)
        : m_h5MeshFile(h5MeshFile),
          m_vertexCount(mesh.getVertexCount()),
          m_spaceDimension(mesh.getSpaceDimension()),
          m_cellCount(mesh.getCellCount()),
          m_xdmfTopologySize(0)
      {
        for (auto it = mesh.getCell(); !it.end(); ++it)
          m_xdmfTopologySize += 1 + it->getVertices().size();
      }

      template <class GridFunctionType>
      XDMF& addGridFunction(
          const std::string& name,
          const GridFunctionType& gf,
          const boost::filesystem::path& h5File)
      {
        Attribute attr;
        attr.name = name;
        attr.h5File = h5File;
        attr.datasetPath = "/GridFunction/Values";
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

        writeTopology(os);
        writeGeometry(os);
        writeAttributes(os);

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
           << m_h5MeshFile.string() << ":/Mesh/XDMFTopology</" << Keyword::DataItem << ">\n";
        os << "      </" << Keyword::Topology << ">\n";
      }

      void writeGeometry(std::ostream& os) const
      {
        os << "      <" << Keyword::Geometry << " GeometryType=\"" << geometryType() << "\">\n";
        os << "        <" << Keyword::DataItem
           << " Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\""
           << m_vertexCount << " " << m_spaceDimension << "\">"
           << m_h5MeshFile.string() << ":/Mesh/Vertices</" << Keyword::DataItem << ">\n";
        os << "      </" << Keyword::Geometry << ">\n";
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
