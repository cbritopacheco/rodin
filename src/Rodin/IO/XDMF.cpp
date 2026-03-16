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

#include "XDMF.h"
#include "HDF5.h"
#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/MemberFunctionException.h"

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
  std::string makeRankSuffix(const std::string& rank)
  {
    return rank.empty() ? std::string() : ".r" + rank;
  }

  static
  std::string expandPattern(
      const std::string& pattern,
      const std::string& stem,
      const std::string& grid,
      const std::string& name,
      const std::string& index,
      const std::string& rank = "")
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
    replace("{rank}", rank);
    replace("{rank_suffix}", makeRankSuffix(rank));

    // Normalize repeated separators introduced by empty placeholders.
    while (result.find("..") != std::string::npos)
      result.replace(result.find(".."), 2, ".");

    while (result.find("/./") != std::string::npos)
      result.replace(result.find("/./"), 3, "/");

    while (result.find("//") != std::string::npos)
      result.replace(result.find("//"), 2, "/");

    // Remove a remaining dot just before ".h5"
    const auto dotExt = result.find(".h5");
    if (dotExt != std::string::npos && dotExt > 0 && result[dotExt - 1] == '.')
      result.erase(dotExt - 1, 1);

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

  struct XDMFMeshMeta
  {
    size_t vertexCount = 0;
    size_t cellCount = 0;
    size_t meshDimension = 0;
    size_t spaceDimension = 0;
    size_t topologySize = 0;
  };

  static
  XDMFMeshMeta readXDMFMeshMeta(const boost::filesystem::path& filename)
  {
    XDMFMeshMeta meta;

    const auto file = HDF5::File(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    if (!file)
    {
      Alert::Exception()
        << "Failed to open HDF5 XDMF mesh file: " << filename
        << Alert::Raise;
    }

    {
      const auto shape = HDF5::readMatrixShape(file.get(), HDF5::Path::MeshGeometryVertices);
      meta.vertexCount = static_cast<size_t>(shape.first);
      meta.spaceDimension = static_cast<size_t>(shape.second);
    }

    meta.topologySize = static_cast<size_t>(
        HDF5::readScalarDataset<HDF5::U64>(file.get(), HDF5::Path::MeshXDMFTopologySize));

    size_t maxDim = 0;
    for (;;)
    {
      const auto path = HDF5::attributePath(maxDim);
      if (!HDF5::exists(file.get(), path))
        break;
      ++maxDim;
    }

    if (maxDim == 0)
    {
      Alert::Exception()
        << "Failed to infer mesh dimension from HDF5 attributes in file: " << filename
        << Alert::Raise;
    }

    meta.meshDimension = maxDim - 1;

    {
      const auto attrs = HDF5::readVectorDataset<HDF5::U64>(
          file.get(),
          HDF5::attributePath(meta.meshDimension));
      meta.cellCount = attrs.size();
    }

    return meta;
  }

  static
  std::string makeGridPieceName(const std::string& gridName, size_t rank)
  {
    if (gridName.empty())
      return "r" + std::to_string(rank);
    return gridName + "_r" + std::to_string(rank);
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

  XDMF::Grid& XDMF::Grid::reset()
  {
    auto& gr = m_owner->m_grids[m_index];
    gr.mesh = nullptr;
    gr.sourceMesh = nullptr;
    gr.options = GridOptions{};
    gr.staticMeshWritten = false;
    gr.staticMeshFile.clear();
    gr.attributes.clear();
    gr.snapshots.clear();
    return *this;
  }

  XDMF::Grid& XDMF::Grid::clear()
  {
    m_owner->m_grids[m_index].attributes.clear();
    return *this;
  }

  // ---- XDMF ----------------------------------------------------------------

  XDMF::XDMF(const boost::filesystem::path& stem)
    : m_stem(stem)
  {}

  XDMF::XDMF(const boost::filesystem::path& stem,
             size_t rank, size_t numRanks, size_t rootRank)
    : m_stem(stem),
      m_distributed(true),
      m_rank(rank),
      m_numRanks(numRanks),
      m_rootRank(rootRank)
  {
    // Use {rank_suffix} so that the same expansion logic can omit the rank
    // cleanly when rank == "".
    m_patterns.staticMesh    = "{stem}.{grid}{rank_suffix}.mesh.h5";
    m_patterns.transientMesh = "{stem}.{grid}{rank_suffix}.mesh.{index}.h5";
    m_patterns.attribute     = "{stem}.{grid}.{name}{rank_suffix}.{index}.h5";
  }

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
    const std::string rankStr  = m_distributed ? std::to_string(m_rank) : std::string();

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
      snapshot.vertexCount = gr.mesh->getVertexCount();
      snapshot.cellCount = gr.mesh->getCellCount();
      snapshot.meshDimension = gr.mesh->getDimension();
      snapshot.spaceDimension = gr.mesh->getSpaceDimension();
      snapshot.topologySize = HDF5::getXDMFMixedTopologySize(*gr.mesh);

      if (gr.options.meshPolicy == MeshPolicy::Static)
      {
        if (!gr.staticMeshWritten)
        {
          const auto meshFile = expandPattern(
              patterns.staticMesh,
              stemStr,
              gr.name,
              "",
              "",
              rankStr);
          const auto meshPath = m_stem.parent_path() / meshFile;

          HDF5::writeXDMFMesh(meshPath, *gr.mesh);

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
            indexStr,
            rankStr);
        const auto meshPath = m_stem.parent_path() / meshFile;

        HDF5::writeXDMFMesh(meshPath, *gr.mesh);

        snapshot.meshFile = meshFile;
      }

      {
        SnapshotRecord::MeshAttributeRecord meshAttr;
        meshAttr.name = "Region";
        meshAttr.center = Center::Cell;
        meshAttr.topologicalDimension = snapshot.meshDimension;
        snapshot.meshAttributes.push_back(std::move(meshAttr));
      }

      snapshot.attributes.clear();
      snapshot.attributes.reserve(gr.attributes.size());

      for (const auto& attr : gr.attributes)
      {
        const auto attrFile = expandPattern(
            patterns.attribute,
            stemStr,
            gr.name,
            attr.name,
            indexStr,
            rankStr);
        const auto attrPath = m_stem.parent_path() / attrFile;

        attr.write(attrPath, attr.center);

        SnapshotRecord::AttributeRecord snapAttr;
        snapAttr.name = attr.name;
        snapAttr.center = attr.center;
        snapAttr.dimension = attr.dimension;
        snapAttr.file = attrFile;
        snapshot.attributes.push_back(std::move(snapAttr));
      }

      gr.snapshots.push_back(std::move(snapshot));
    }

    ++m_snapshotCount;
    return *this;
  }

  void XDMF::writeUniformGrid(
      std::ostream& os,
      const std::string& gridName,
      const SnapshotRecord& snap,
      size_t bi) const
  {
    const auto meshH5 = snap.meshFile.string();

    os << indent(bi) << "<Grid Name=\"" << gridName << "\" GridType=\"Uniform\">\n";
    os << indent(bi + 1) << "<Time Value=\"" << snap.time << "\" />\n";

    os << indent(bi + 1) << "<Topology TopologyType=\"Mixed\" NumberOfElements=\""
       << snap.cellCount << "\">\n";
    // Precision="8" is required here because the HDF5 topology dataset is U64.
    // Without Precision, some XDMF readers may assume 32-bit UInt and misread
    // the dataset.
    os << indent(bi + 2) << "<DataItem Format=\"HDF\" NumberType=\"UInt\" Precision=\"8\" Dimensions=\""
       << snap.topologySize << "\">"
       << meshH5 << ":" << HDF5::Path::MeshXDMFTopology
       << "</DataItem>\n";
    os << indent(bi + 1) << "</Topology>\n";

    os << indent(bi + 1) << "<Geometry GeometryType=\""
       << getGeometryType(snap.spaceDimension) << "\">\n";
    os << indent(bi + 2) << "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\""
       << snap.vertexCount << " " << snap.spaceDimension
       << "\">"
       << meshH5 << ":" << HDF5::Path::MeshGeometryVertices
       << "</DataItem>\n";
    os << indent(bi + 1) << "</Geometry>\n";

    for (const auto& attr : snap.meshAttributes)
    {
      const char* centerStr = (attr.center == Center::Node) ? "Node" : "Cell";
      std::ostringstream dimStr;
      dimStr << ((attr.center == Center::Node) ? snap.vertexCount : snap.cellCount);

      os << indent(bi + 1) << "<Attribute Name=\"" << attr.name
         << "\" AttributeType=\"Scalar\" Center=\"" << centerStr << "\">\n";
      // Same reason as topology: Region is stored as U64 in HDF5.
      os << indent(bi + 2) << "<DataItem Format=\"HDF\" NumberType=\"UInt\" Precision=\"8\" Dimensions=\""
         << dimStr.str() << "\">"
         << meshH5 << ":" << HDF5::attributePath(attr.topologicalDimension)
         << "</DataItem>\n";
      os << indent(bi + 1) << "</Attribute>\n";
    }

    for (const auto& attr : snap.attributes)
    {
      const auto attrH5 = attr.file.string();
      const char* centerStr = (attr.center == Center::Node) ? "Node" : "Cell";
      const char* attrType  = (attr.dimension == 1) ? "Scalar" : "Vector";

      std::ostringstream dimStr;
      const size_t count = (attr.center == Center::Node)
          ? snap.vertexCount : snap.cellCount;
      if (attr.dimension == 1)
        dimStr << count;
      else
        dimStr << count << " " << attr.dimension;

      os << indent(bi + 1) << "<Attribute Name=\"" << attr.name
         << "\" AttributeType=\"" << attrType
         << "\" Center=\"" << centerStr << "\">\n";
      os << indent(bi + 2) << "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\""
         << dimStr.str() << "\">"
         << attrH5 << ":" << HDF5::Path::GridFunctionValuesData
         << "</DataItem>\n";
      os << indent(bi + 1) << "</Attribute>\n";
    }

    os << indent(bi) << "</Grid>\n";
  }

  void XDMF::flush() const
  {
    if (m_distributed && m_rank != m_rootRank)
      return;

    const std::string stemStr = m_stem.filename().string();
    const auto xdmfFile = expandPattern(m_patterns.xdmf, stemStr, "", "", "");
    const auto xdmfPath = m_stem.parent_path() / xdmfFile;

    std::ofstream os(xdmfPath.string(), std::ios::out | std::ios::trunc);
    if (!os)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to open XDMF output file: " << xdmfPath
        << Alert::Raise;
    }

    writeXMLHeader(os);
    os << indent(1) << "<Domain>\n";

    if (!m_distributed)
    {
      for (const auto& gr : m_grids)
      {
        const std::string gridName = gr.name.empty() ? "mesh" : gr.name;

        os << indent(2)
           << "<Grid Name=\"" << gridName
           << "\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

        for (const auto& snap : gr.snapshots)
          writeUniformGrid(os, gridName, snap, 3);

        os << indent(2) << "</Grid>\n";
      }
    }
    else
    {
      for (const auto& gr : m_grids)
      {
        const std::string gridName = gr.name.empty() ? "mesh" : gr.name;
        const auto& patterns = gr.options.patterns ? *gr.options.patterns : m_patterns;

        os << indent(2)
           << "<Grid Name=\"" << gridName
           << "\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

        for (size_t s = 0; s < gr.snapshots.size(); ++s)
        {
          const auto& snap = gr.snapshots[s];
          const std::string indexStr = padIndex(s, m_padding);

          os << indent(3) << "<Grid Name=\"Step_" << s
             << "\" GridType=\"Collection\" CollectionType=\"Spatial\">\n";
          os << indent(4) << "<Time Value=\"" << snap.time << "\" />\n";

          for (size_t r = 0; r < m_numRanks; ++r)
          {
            const std::string rStr = std::to_string(r);

            const std::string meshH5 =
              (gr.options.meshPolicy == MeshPolicy::Static)
              ? expandPattern(patterns.staticMesh, stemStr, gr.name, "", "", rStr)
              : expandPattern(patterns.transientMesh, stemStr, gr.name, "", indexStr, rStr);

            const auto meshPath = m_stem.parent_path() / meshH5;
            const auto meta = readXDMFMeshMeta(meshPath);

            os << indent(5) << "<Grid Name=\"" << makeGridPieceName(gridName, r)
               << "\" GridType=\"Uniform\">\n";

            os << indent(6) << "<Topology TopologyType=\"Mixed\" NumberOfElements=\""
               << meta.cellCount << "\">\n";
            os << indent(7) << "<DataItem Format=\"HDF\" NumberType=\"UInt\" Precision=\"8\" Dimensions=\""
               << meta.topologySize << "\">"
               << meshH5 << ":" << HDF5::Path::MeshXDMFTopology
               << "</DataItem>\n";
            os << indent(6) << "</Topology>\n";

            os << indent(6) << "<Geometry GeometryType=\""
               << getGeometryType(meta.spaceDimension) << "\">\n";
            os << indent(7) << "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\""
               << meta.vertexCount << " " << meta.spaceDimension
               << "\">"
               << meshH5 << ":" << HDF5::Path::MeshGeometryVertices
               << "</DataItem>\n";
            os << indent(6) << "</Geometry>\n";

            for (const auto& attr : snap.meshAttributes)
            {
              const char* centerStr = (attr.center == Center::Node) ? "Node" : "Cell";
              std::ostringstream dimStr;
              dimStr << ((attr.center == Center::Node) ? meta.vertexCount : meta.cellCount);

              os << indent(6) << "<Attribute Name=\"" << attr.name
                 << "\" AttributeType=\"Scalar\" Center=\"" << centerStr << "\">\n";
              os << indent(7) << "<DataItem Format=\"HDF\" NumberType=\"UInt\" Precision=\"8\" Dimensions=\""
                 << dimStr.str() << "\">"
                 << meshH5 << ":" << HDF5::attributePath(attr.topologicalDimension)
                 << "</DataItem>\n";
              os << indent(6) << "</Attribute>\n";
            }

            for (const auto& attr : snap.attributes)
            {
              const auto attrH5 = expandPattern(
                  patterns.attribute, stemStr, gr.name, attr.name, indexStr, rStr);

              const char* centerStr = (attr.center == Center::Node) ? "Node" : "Cell";
              const char* attrType  = (attr.dimension == 1) ? "Scalar" : "Vector";

              std::ostringstream dimStr;
              const size_t count = (attr.center == Center::Node)
                  ? meta.vertexCount : meta.cellCount;
              if (attr.dimension == 1)
                dimStr << count;
              else
                dimStr << count << " " << attr.dimension;

              os << indent(6) << "<Attribute Name=\"" << attr.name
                 << "\" AttributeType=\"" << attrType
                 << "\" Center=\"" << centerStr << "\">\n";
              os << indent(7) << "<DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\""
                 << dimStr.str() << "\">"
                 << attrH5 << ":" << HDF5::Path::GridFunctionValuesData
                 << "</DataItem>\n";
              os << indent(6) << "</Attribute>\n";
            }

            os << indent(5) << "</Grid>\n";
          }

          os << indent(3) << "</Grid>\n";
        }

        os << indent(2) << "</Grid>\n";
      }
    }

    os << indent(1) << "</Domain>\n";
    writeXMLFooter(os);
    os.flush();

    if (!os)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed while writing XDMF output file: " << xdmfPath
        << Alert::Raise;
    }
  }

  void XDMF::close()
  {
    if (m_closed)
      return;

    flush();
    m_closed = true;
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

  bool XDMF::isDistributed() const noexcept
  {
    return m_distributed;
  }

  size_t XDMF::getRank() const noexcept
  {
    return m_rank;
  }

  size_t XDMF::getNumRanks() const noexcept
  {
    return m_numRanks;
  }

  size_t XDMF::getRootRank() const noexcept
  {
    return m_rootRank;
  }
}
