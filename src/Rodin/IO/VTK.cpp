/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "VTK.h"

#include <sstream>
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::IO::VTK
{
  std::istream& getline(std::istream& is, std::string& line, size_t& currentLineNumber)
  {
    currentLineNumber++;
    return std::getline(is, line);
  }

  std::string skipEmptyLinesAndComments(std::istream& is, size_t& currentLineNumber)
  {
    std::string line;
    while (getline(is, line, currentLineNumber))
    {
      // Trim whitespace
      line.erase(0, line.find_first_not_of(" \t\r\n"));
      line.erase(line.find_last_not_of(" \t\r\n") + 1);
      
      // Skip empty lines and comments (starting with #)
      if (!line.empty() && line[0] != '#')
        break;
    }
    return line;
  }
}

namespace Rodin::IO
{
  void MeshLoader<FileFormat::VTKLegacy, Context::Local>::load(std::istream& is)
  {
    readHeader(is);
    readDataset(is);
    readPoints(is);  // This also sets m_spaceDimension and adds vertices
    m_build.initialize(m_spaceDimension);
    readCells(is);
    readCellTypes(is);
    
    // Determine mesh dimension from cell types
    size_t meshDimension = 0;
    for (const auto& cell : m_cells)
    {
      auto geom = VTK::getGeometry(cell.cellType);
      if (geom)
      {
        size_t cellDim = Geometry::Polytope::Traits(*geom).getDimension();
        meshDimension = std::max(meshDimension, cellDim);
      }
    }
    
    // Reserve space for cells of highest dimension
    m_build.reserve(meshDimension, m_numCells);
    
    // Process cells to build connectivity
    for (const auto& cell : m_cells)
    {
      auto geom = VTK::getGeometry(cell.cellType);
      if (geom)
      {
        m_build.polytope(*geom, cell.vertices);
      }
    }
    
    getObject() = m_build.finalize();
  }

  void MeshLoader<FileFormat::VTKLegacy, Context::Local>::readHeader(std::istream& is)
  {
    // Read VTK version header
    std::string line = VTK::skipEmptyLinesAndComments(is, m_currentLineNumber);
    if (line.find("# vtk DataFile Version") != 0)
    {
      Alert::Exception() << "Expected VTK header, got: " << line << Alert::Raise;
    }
    
    // Read title
    line = VTK::skipEmptyLinesAndComments(is, m_currentLineNumber);
    
    // Read data format (ASCII/BINARY)
    line = VTK::skipEmptyLinesAndComments(is, m_currentLineNumber);
    if (line != "ASCII")
    {
      Alert::Exception() << "Only ASCII format is supported, got: " << line << Alert::Raise;
    }
  }

  void MeshLoader<FileFormat::VTKLegacy, Context::Local>::readDataset(std::istream& is)
  {
    std::string line = VTK::skipEmptyLinesAndComments(is, m_currentLineNumber);
    if (line.find("DATASET") != 0)
    {
      Alert::Exception() << "Expected DATASET keyword, got: " << line << Alert::Raise;
    }
    
    if (line.find("UNSTRUCTURED_GRID") == std::string::npos)
    {
      Alert::Exception() << "Only UNSTRUCTURED_GRID datasets are supported" << Alert::Raise;
    }
  }

  void MeshLoader<FileFormat::VTKLegacy, Context::Local>::readPoints(std::istream& is)
  {
    std::string line = VTK::skipEmptyLinesAndComments(is, m_currentLineNumber);
    
    std::istringstream iss(line);
    std::string keyword, dataType;
    iss >> keyword >> m_numPoints >> dataType;
    
    if (keyword != "POINTS")
    {
      Alert::Exception() << "Expected POINTS keyword, got: " << keyword << Alert::Raise;
    }
    
    if (dataType != "float" && dataType != "double")
    {
      Alert::Exception() << "Unsupported point data type: " << dataType << Alert::Raise;
    }

    // Reserve space for nodes
    m_build.nodes(m_numPoints);
    
    // Read point coordinates and add them to the builder
    for (size_t i = 0; i < m_numPoints; i++)
    {
      line = VTK::skipEmptyLinesAndComments(is, m_currentLineNumber);
      std::istringstream coords(line);
      
      Math::SpatialPoint point(3); // VTK always uses 3D coordinates
      coords >> point(0) >> point(1) >> point(2);
      
      // Determine actual space dimension from first point
      if (i == 0)
      {
        m_spaceDimension = 3;
        // For simplicity, assume 3D unless all z-coordinates are zero
        // A more sophisticated approach would check all points
        if (point(2) == 0.0)
        {
          m_spaceDimension = 2;
          if (point(1) == 0.0)
            m_spaceDimension = 1;
        }
      }
      
      // Add vertex to builder
      m_build.vertex(std::move(point));
    }
  }

  void MeshLoader<FileFormat::VTKLegacy, Context::Local>::readCells(std::istream& is)
  {
    std::string line = VTK::skipEmptyLinesAndComments(is, m_currentLineNumber);
    
    std::istringstream iss(line);
    std::string keyword;
    size_t totalSize;
    iss >> keyword >> m_numCells >> totalSize;
    
    if (keyword != "CELLS")
    {
      Alert::Exception() << "Expected CELLS keyword, got: " << keyword << Alert::Raise;
    }
    
    m_cells.resize(m_numCells);
    
    // Read cell connectivity
    for (size_t i = 0; i < m_numCells; i++)
    {
      line = VTK::skipEmptyLinesAndComments(is, m_currentLineNumber);
      auto result = VTK::ParseCell()(line.begin(), line.end());
      if (result)
      {
        m_cells[i] = *result;
      }
      else
      {
        Alert::Exception() << "Failed to parse cell " << i << Alert::Raise;
      }
    }
  }

  void MeshLoader<FileFormat::VTKLegacy, Context::Local>::readCellTypes(std::istream& is)
  {
    std::string line = VTK::skipEmptyLinesAndComments(is, m_currentLineNumber);
    
    std::istringstream iss(line);
    std::string keyword;
    size_t numTypes;
    iss >> keyword >> numTypes;
    
    if (keyword != "CELL_TYPES")
    {
      Alert::Exception() << "Expected CELL_TYPES keyword, got: " << keyword << Alert::Raise;
    }
    
    if (numTypes != m_numCells)
    {
      Alert::Exception() << "Number of cell types (" << numTypes 
                         << ") does not match number of cells (" << m_numCells << ")" << Alert::Raise;
    }
    
    // Read cell types
    for (size_t i = 0; i < m_numCells; i++)
    {
      line = VTK::skipEmptyLinesAndComments(is, m_currentLineNumber);
      auto result = VTK::ParseUnsignedInteger()(line.begin(), line.end());
      if (result)
      {
        m_cells[i].cellType = static_cast<VTK::CellType>(*result);
      }
      else
      {
        Alert::Exception() << "Failed to parse cell type " << i << Alert::Raise;
      }
    }
  }

  void MeshPrinter<FileFormat::VTKLegacy, Context::Local>::print(std::ostream& os)
  {
    printHeader(os);
    printDataset(os);
    printPoints(os);
    printCells(os);
    printCellTypes(os);
  }

  void MeshPrinter<FileFormat::VTKLegacy, Context::Local>::printHeader(std::ostream& os)
  {
    os << "# vtk DataFile Version 3.0\n";
    os << "Rodin Mesh\n";
    os << "ASCII\n";
  }

  void MeshPrinter<FileFormat::VTKLegacy, Context::Local>::printDataset(std::ostream& os)
  {
    os << "DATASET UNSTRUCTURED_GRID\n";
  }

  void MeshPrinter<FileFormat::VTKLegacy, Context::Local>::printPoints(std::ostream& os)
  {
    const auto& mesh = getObject();
    const size_t numVertices = mesh.getVertexCount();
    
    os << "POINTS " << numVertices << " double\n";
    
    for (size_t i = 0; i < numVertices; i++)
    {
      auto coords = mesh.getVertexCoordinates(i);
      
      // Always output 3D coordinates (pad with zeros if needed)
      os << coords(0) << " ";
      if (coords.size() > 1)
        os << coords(1) << " ";
      else
        os << "0.0 ";
      
      if (coords.size() > 2)
        os << coords(2) << "\n";
      else
        os << "0.0\n";
    }
  }

  void MeshPrinter<FileFormat::VTKLegacy, Context::Local>::printCells(std::ostream& os)
  {
    const auto& mesh = getObject();
    const size_t dim = mesh.getDimension();
    const size_t numCells = mesh.getPolytopeCount(dim);
    
    // Calculate total size (sum of vertices per cell + 1 for each cell)
    size_t totalSize = 0;
    for (size_t i = 0; i < numCells; i++)
    {
      auto geom = mesh.getGeometry(dim, i);
      totalSize += Geometry::Polytope::Traits(geom).getVertexCount() + 1;
    }
    
    os << "CELLS " << numCells << " " << totalSize << "\n";
    
    for (size_t i = 0; i < numCells; i++)
    {
      auto polytope = mesh.getPolytope(dim, i);
      auto vertices = polytope.getVertices();
      
      os << vertices.size();
      for (size_t j = 0; j < vertices.size(); j++)
      {
        os << " " << vertices[j];
      }
      os << "\n";
    }
  }

  void MeshPrinter<FileFormat::VTKLegacy, Context::Local>::printCellTypes(std::ostream& os)
  {
    const auto& mesh = getObject();
    const size_t dim = mesh.getDimension();
    const size_t numCells = mesh.getPolytopeCount(dim);
    
    os << "CELL_TYPES " << numCells << "\n";
    
    for (size_t i = 0; i < numCells; i++)
    {
      auto geom = mesh.getGeometry(dim, i);
      auto cellType = VTK::getCellType(geom);
      if (cellType)
      {
        os << static_cast<int>(*cellType) << "\n";
      }
      else
      {
        Alert::Exception() << "Unsupported geometry type for VTK output" << Alert::Raise;
      }
    }
  }
}