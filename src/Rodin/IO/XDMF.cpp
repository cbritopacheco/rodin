/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cassert>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "XDMF.h"
#include "HDF5.h"
#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Geometry/Mesh.h"

namespace Rodin::IO
{
  // ---- helpers (file-local) ------------------------------------------------

  static
  std::string padIndex(size_t index, size_t width)
  {
    std::ostringstream oss;
    oss << std::setw(static_cast<int>(width)) << std::setfill('0') << index;
    return oss.str();
  }

  static
  std::string expandPattern(
      const std::string& pattern,
      const std::string& stem,
      const std::string& grid,
      const std::string& name,
      const std::string& index)
  {
    std::string result = pattern;

    auto replace = [&](const std::string& placeholder, const std::string& value)
    {
      for (std::string::size_type pos = 0; ;)
      {
        pos = result.find(placeholder, pos);
        if (pos == std::string::npos)
          break;
        result.replace(pos, placeholder.size(), value);
        pos += value.size();
      }
    };

    replace("{stem}", stem);
    replace("{grid}", grid);
    replace("{name}", name);
    replace("{index}", index);

    // Normalize separators when some placeholders are empty.
    while (result.find("..") != std::string::npos)
      result.replace(result.find(".."), 2, ".");

    // Remove "/./" if patterns ever introduce it.
    while (result.find("/./") != std::string::npos)
      result.replace(result.find("/./"), 3, "/");

    // Remove a trailing dot before an extension, e.g. "foo..h5" -> "foo.h5".
    const auto dotExt = result.find(".h5");
    if (dotExt != std::string::npos && dotExt > 0 && result[dotExt - 1] == '.')
      result.erase(dotExt - 1, 1);

    // More generally, avoid trailing separators.
    if (!result.empty() && result.back() == '.')
      result.pop_back();

    return result;
  }

  static
  void writeXMLHeader(std::ostream& os)
  {
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    os << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    os << "<Xdmf Version=\"3.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n";
  }

  static
  void writeXMLFooter(std::ostream& os)
  {
    os << "</Xdmf>\n";
  }

  static
  std::string indent(size_t level)
  {
    return std::string(level * 2, ' ');
  }

  static
  const char* getGeometryType(size_t sdim)
  {
    switch (sdim)
    {
      case 1:
        return "X";
      case 2:
        return "XY";
      case 3:
        return "XYZ";
      default:
        Alert::Exception()
          << "Unsupported space dimension for XDMF geometry: " << sdim
          << Alert::Raise;
    }
    assert(false);
    return nullptr;
  }

  // ---- XDMF::Grid ----------------------------------------------------------

  XDMF::Grid::Grid(XDMF& owner, size_t index) noexcept
    : m_owner(&owner),
      m_index(index)
  {}

  std::string_view XDMF::Grid::getName() const noexcept
  {
    return m_owner->m_grids[m_index].name;
  }

  XDMF::Grid& XDMF::Grid::setOptions(const GridOptions& options)
  {
    m_owner->m_grids[m_index].options = options;
    return *this;
  }

  const XDMF::GridOptions& XDMF::Grid::getOptions() const noexcept
  {
    return m_owner->m_grids[m_index].options;
  }

  bool XDMF::Grid::hasMesh() const noexcept
  {
    return m_owner->m_grids[m_index].mesh != nullptr;
  }

  size_t XDMF::Grid::getAttributeCount() const noexcept
  {
    return m_owner->m_grids[m_index].attributes.size();
  }

  // ---- XDMF ----------------------------------------------------------------

  XDMF::XDMF(const boost::filesystem::path& stem)
    : m_stem(stem)
  {}

  XDMF::~XDMF()
  {
    if (!m_closed)
    {
      try
      {
        close();
      }
      catch (...)
      {
        // Suppress exceptions in destructor.
      }
    }
  }

  const boost::filesystem::path& XDMF::getStem() const noexcept
  {
    return m_stem;
  }

  XDMF& XDMF::setFilePatterns(const FilePatterns& patterns)
  {
    m_patterns = patterns;
    return *this;
  }

  const XDMF::FilePatterns& XDMF::getFilePatterns() const noexcept
  {
    return m_patterns;
  }

  XDMF& XDMF::setPadding(size_t digits)
  {
    m_padding = digits;
    return *this;
  }

  size_t XDMF::getPadding() const noexcept
  {
    return m_padding;
  }

  XDMF::Grid XDMF::grid()
  {
    return grid("");
  }

  XDMF::Grid XDMF::grid(const std::string& name)
  {
    for (size_t i = 0; i < m_grids.size(); ++i)
    {
      if (m_grids[i].name == name)
        return Grid(*this, i);
    }

    m_grids.push_back(GridRecord{});
    m_grids.back().name = name;
    return Grid(*this, m_grids.size() - 1);
  }

  XDMF& XDMF::setOptions(const GridOptions& options)
  {
    grid().setOptions(options);
    return *this;
  }

  XDMF& XDMF::write()
  {
    return write(static_cast<Real>(m_snapshotCount));
  }

  XDMF& XDMF::write(Real time)
  {
    if (m_closed)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Cannot write to a closed XDMF writer."
        << Alert::Raise;
    }

    const std::string stemStr  = m_stem.filename().string();
    const std::string indexStr = padIndex(m_snapshotCount, m_padding);

    for (auto& gr : m_grids)
    {
      if (!gr.mesh)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Grid \"" << gr.name << "\" has no mesh set."
          << Alert::Raise;
      }

      const auto& patterns = gr.options.patterns ? *gr.options.patterns : m_patterns;

      SnapshotRecord snapshot;
      snapshot.time = time;

      // --- export mesh ------------------------------------------------------
      if (gr.options.meshPolicy == MeshPolicy::Static)
      {
        if (!gr.staticMeshWritten)
        {
          const auto meshFile = expandPattern(
              patterns.staticMesh,
              stemStr,
              gr.name,
              "",
              "");
          const auto meshPath = m_stem.parent_path() / meshFile;

          IO::MeshPrinter<IO::FileFormat::HDF5, Context::Local>(
              static_cast<const Geometry::Mesh<Context::Local>&>(*gr.mesh))
            .setXDMF(true)
            .print(meshPath);
          gr.staticMeshFile = meshFile;
          gr.staticMeshWritten = true;
        }

        snapshot.meshFile = gr.staticMeshFile;
      }
      else
      {
        const auto meshFile = expandPattern(
            patterns.transientMesh,
            stemStr,
            gr.name,
            "",
            indexStr);
        const auto meshPath = m_stem.parent_path() / meshFile;

      IO::MeshPrinter<IO::FileFormat::HDF5, Context::Local>(
          static_cast<const Geometry::Mesh<Context::Local>&>(*gr.mesh))
        .setXDMF(true)
        .print(meshPath);
        snapshot.meshFile = meshFile;
      }

      // --- export attributes ------------------------------------------------
      for (const auto& attr : gr.attributes)
      {
        const auto attrFile = expandPattern(
            patterns.attribute,
            stemStr,
            gr.name,
            attr.name,
            indexStr);
        const auto attrPath = m_stem.parent_path() / attrFile;

        attr.write(attrPath, attr.center);
        snapshot.attributeFiles.push_back(attrFile);
      }

      gr.snapshots.push_back(std::move(snapshot));
    }

    ++m_snapshotCount;
    return *this;
  }

  void XDMF::close()
  {
    if (m_closed)
      return;

    m_closed = true;

    const std::string stemStr = m_stem.filename().string();
    const auto xdmfFile = expandPattern(m_patterns.xdmf, stemStr, "", "", "");
    const auto xdmfPath = m_stem.parent_path() / xdmfFile;

    std::ofstream os(xdmfPath.string());
    if (!os)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to open XDMF output file: " << xdmfPath
        << Alert::Raise;
    }

    writeXMLHeader(os);
    os << indent(1) << "<Domain>\n";

    for (const auto& gr : m_grids)
    {
      os << indent(2)
         << "<Grid Name=\"" << gr.name
         << "\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

      for (const auto& snap : gr.snapshots)
      {
        const auto meshH5 = snap.meshFile.string();
        const auto topologySize = HDF5::getXDMFMixedTopologySize(*gr.mesh);

        os << indent(3) << "<Grid Name=\"" << gr.name << "\" GridType=\"Uniform\">\n";

        os << indent(4) << "<Time Value=\"" << snap.time << "\" />\n";

        os << indent(4) << "<Topology TopologyType=\"Mixed\" NumberOfElements=\""
           << gr.mesh->getCellCount() << "\">\n";
        os << indent(5) << "<DataItem Format=\"HDF\" NumberType=\"UInt\" Dimensions=\""
           << topologySize << "\">"
           << meshH5 << ":" << HDF5::Path::MeshXDMFTopology
           << "</DataItem>\n";
        os << indent(4) << "</Topology>\n";

        os << indent(4) << "<Geometry GeometryType=\""
           << getGeometryType(gr.mesh->getSpaceDimension()) << "\">\n";
        os << indent(5) << "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\""
           << gr.mesh->getVertexCount() << " " << gr.mesh->getSpaceDimension()
           << "\">"
           << meshH5 << ":" << HDF5::Path::MeshGeometryVertices
           << "</DataItem>\n";
        os << indent(4) << "</Geometry>\n";

        for (size_t a = 0; a < gr.attributes.size(); ++a)
        {
          const auto& attr = gr.attributes[a];
          const auto attrH5 = snap.attributeFiles[a].string();
          const char* centerStr = (attr.center == Center::Node) ? "Node" : "Cell";
          const char* attrType  = (attr.dimension == 1) ? "Scalar" : "Vector";

          std::ostringstream dimStr;
          if (attr.center == Center::Node)
          {
            if (attr.dimension == 1)
              dimStr << gr.mesh->getVertexCount();
            else
              dimStr << gr.mesh->getVertexCount() << " " << attr.dimension;
          }
          else
          {
            if (attr.dimension == 1)
              dimStr << gr.mesh->getCellCount();
            else
              dimStr << gr.mesh->getCellCount() << " " << attr.dimension;
          }

          os << indent(4) << "<Attribute Name=\"" << attr.name
             << "\" AttributeType=\"" << attrType
             << "\" Center=\"" << centerStr << "\">\n";
          os << indent(5) << "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\""
             << dimStr.str() << "\">"
             << attrH5 << ":" << HDF5::Path::GridFunctionValuesData
             << "</DataItem>\n";
          os << indent(4) << "</Attribute>\n";
        }

        os << indent(3) << "</Grid>\n";
      }

      os << indent(2) << "</Grid>\n";
    }

    os << indent(1) << "</Domain>\n";
    writeXMLFooter(os);
  }

  bool XDMF::isClosed() const noexcept
  {
    return m_closed;
  }

  size_t XDMF::getSnapshotCount() const noexcept
  {
    return m_snapshotCount;
  }

  size_t XDMF::getGridCount() const noexcept
  {
    return m_grids.size();
  }
}
