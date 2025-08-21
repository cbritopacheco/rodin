/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "GMSH.h"

#include <sstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "Rodin/Alert.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::IO::GMSH
{
  std::istream& getline(std::istream& is, std::string& line, size_t& currentLineNumber)
  {
    currentLineNumber++;
    return std::getline(is, line);
  }

  std::string skipEmptyLines(std::istream& is, size_t& currentLineNumber)
  {
    std::string line;
    while (GMSH::getline(is, line, currentLineNumber))
    {
      boost::algorithm::trim(line);
      if (!line.empty())
        break;
    }
    return line;
  }
}

namespace Rodin::IO
{
  void MeshLoader<FileFormat::GMSH, Context::Local>::load(std::istream& is)
  {
    m_currentLineNumber = 0;
    
    readMeshFormat(is);
    readNodes(is);
    readElements(is);
    
    this->getObject() = m_builder.finalize();
  }

  void MeshLoader<FileFormat::GMSH, Context::Local>::readMeshFormat(std::istream& is)
  {
    std::string line = GMSH::skipEmptyLines(is, m_currentLineNumber);
    
    if (line != "$MeshFormat")
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Expected $MeshFormat section, got: " << line << Alert::Raise;
    }
    
    // Read format line
    line = GMSH::skipEmptyLines(is, m_currentLineNumber);
    std::istringstream iss(line);
    double version;
    int fileType, dataSize;
    
    if (!(iss >> version >> fileType >> dataSize))
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to parse mesh format line: " << line << Alert::Raise;
    }
    
    m_header.version.major = static_cast<int>(version);
    m_header.version.minor = static_cast<int>((version - m_header.version.major) * 10 + 0.5);
    m_header.fileType = fileType;
    m_header.dataSize = dataSize;
    
    if (m_header.version.major != 2)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Unsupported GMSH format version: " << version 
        << ". Only version 2.x is supported." << Alert::Raise;
    }
    
    // Read end section
    line = GMSH::skipEmptyLines(is, m_currentLineNumber);
    if (line != "$EndMeshFormat")
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Expected $EndMeshFormat, got: " << line << Alert::Raise;
    }
  }

  void MeshLoader<FileFormat::GMSH, Context::Local>::readNodes(std::istream& is)
  {
    std::string line = GMSH::skipEmptyLines(is, m_currentLineNumber);
    
    if (line != "$Nodes")
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Expected $Nodes section, got: " << line << Alert::Raise;
    }
    
    // Read number of nodes
    line = GMSH::skipEmptyLines(is, m_currentLineNumber);
    size_t numNodes;
    std::istringstream iss(line);
    if (!(iss >> numNodes))
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to parse number of nodes: " << line << Alert::Raise;
    }
    
    // Determine space dimension from first node
    line = GMSH::skipEmptyLines(is, m_currentLineNumber);
    std::istringstream firstNodeIss(line);
    int nodeId;
    Real x, y, z;
    if (!(firstNodeIss >> nodeId >> x >> y >> z))
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to parse first node: " << line << Alert::Raise;
    }
    
    // Determine space dimension based on coordinates
    size_t spaceDim = (std::abs(z) > 1e-14) ? 3 : (std::abs(y) > 1e-14) ? 2 : 1;
    
    // Initialize builder
    m_builder.initialize(spaceDim).nodes(numNodes);
    
    // Add first vertex
    if (spaceDim == 3)
      m_builder.vertex({x, y, z});
    else if (spaceDim == 2)
      m_builder.vertex({x, y});
    else
      m_builder.vertex({x});
    
    // Read remaining nodes
    for (size_t i = 1; i < numNodes; i++)
    {
      line = GMSH::skipEmptyLines(is, m_currentLineNumber);
      std::istringstream nodeIss(line);
      
      if (!(nodeIss >> nodeId >> x >> y >> z))
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Failed to parse node " << i + 1 << ": " << line << Alert::Raise;
      }
      
      // Add vertex
      if (spaceDim == 3)
        m_builder.vertex({x, y, z});
      else if (spaceDim == 2)
        m_builder.vertex({x, y});
      else
        m_builder.vertex({x});
    }
    
    // Read end section
    line = GMSH::skipEmptyLines(is, m_currentLineNumber);
    if (line != "$EndNodes")
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Expected $EndNodes, got: " << line << Alert::Raise;
    }
  }

  void MeshLoader<FileFormat::GMSH, Context::Local>::readElements(std::istream& is)
  {
    std::string line = GMSH::skipEmptyLines(is, m_currentLineNumber);
    
    if (line != "$Elements")
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Expected $Elements section, got: " << line << Alert::Raise;
    }
    
    // Read number of elements
    line = GMSH::skipEmptyLines(is, m_currentLineNumber);
    size_t numElements;
    std::istringstream iss(line);
    if (!(iss >> numElements))
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to parse number of elements: " << line << Alert::Raise;
    }
    
    // Read elements
    for (size_t i = 0; i < numElements; i++)
    {
      line = GMSH::skipEmptyLines(is, m_currentLineNumber);
      std::istringstream elemIss(line);
      
      int elemId, elemType, numTags;
      if (!(elemIss >> elemId >> elemType >> numTags))
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Failed to parse element header " << i + 1 << ": " << line << Alert::Raise;
      }
      
      // Read tags (physical region, elementary entity, etc.)
      std::vector<int> tags(numTags);
      for (int j = 0; j < numTags; j++)
      {
        if (!(elemIss >> tags[j]))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to parse tag " << j << " for element " << i + 1 << Alert::Raise;
        }
      }
      
      // Convert GMSH element type to Polytope type
      auto polytopeType = GMSH::toPolytopeType(static_cast<GMSH::ElementType>(elemType));
      auto traits = Geometry::Polytope::Traits(polytopeType);
      
      // Read node indices
      std::vector<Index> nodeIndices(traits.getVertexCount());
      for (size_t j = 0; j < nodeIndices.size(); j++)
      {
        int nodeId;
        if (!(elemIss >> nodeId))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to parse node " << j << " for element " << i + 1 << Alert::Raise;
        }
        nodeIndices[j] = nodeId - 1; // Convert to 0-based indexing
      }
      
      // Add polytope to mesh
      m_builder.polytope(polytopeType, std::move(nodeIndices));
      
      // Set attribute if specified
      if (numTags > 0)
      {
        Attribute attribute = tags[0]; // Use first tag as attribute
        auto traits = Geometry::Polytope::Traits(polytopeType);
        size_t dimension = traits.getDimension();
        // The polytope index will be the current count before we added this one
        size_t polytopeIndex = m_builder.getConnectivity().getPolytopeCount(dimension) - 1;
        m_builder.attribute({dimension, polytopeIndex}, attribute);
      }
    }
    
    // Read end section
    line = GMSH::skipEmptyLines(is, m_currentLineNumber);
    if (line != "$EndElements")
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Expected $EndElements, got: " << line << Alert::Raise;
    }
  }

  void MeshPrinter<FileFormat::GMSH, Context::Local>::print(std::ostream& os)
  {
    printMeshFormat(os);
    printNodes(os);
    printElements(os);
  }

  void MeshPrinter<FileFormat::GMSH, Context::Local>::printMeshFormat(std::ostream& os)
  {
    os << "$MeshFormat\n";
    os << "2.2 0 8\n";  // Version 2.2, ASCII format, sizeof(double)
    os << "$EndMeshFormat\n";
  }

  void MeshPrinter<FileFormat::GMSH, Context::Local>::printNodes(std::ostream& os)
  {
    const auto& mesh = this->getObject();
    const size_t numVertices = mesh.getVertexCount();
    
    os << "$Nodes\n";
    os << numVertices << "\n";
    
    for (size_t i = 0; i < numVertices; i++)
    {
      auto coords = mesh.getVertexCoordinates(i);
      os << (i + 1) << " ";  // 1-based indexing
      os << coords(0) << " " << coords(1) << " ";
      if (coords.size() > 2)
        os << coords(2);
      else
        os << "0.0";
      os << "\n";
    }
    
    os << "$EndNodes\n";
  }

  void MeshPrinter<FileFormat::GMSH, Context::Local>::printElements(std::ostream& os)
  {
    const auto& mesh = this->getObject();
    
    // Count total elements across all dimensions
    size_t totalElements = 0;
    for (size_t d = 0; d <= mesh.getDimension(); d++)
    {
      totalElements += mesh.getPolytopeCount(d);
    }
    
    os << "$Elements\n";
    os << totalElements << "\n";
    
    size_t elementId = 1;
    
    // Print elements for each dimension
    for (size_t d = 0; d <= mesh.getDimension(); d++)
    {
      for (auto it = mesh.getPolytope(d); it; ++it)
      {
        auto polytopeType = it->getType();
        auto gmshType = static_cast<int>(GMSH::fromPolytopeType(polytopeType));
        auto attribute = it->getAttribute();
        
        os << elementId << " " << gmshType << " 2 " << attribute << " " << attribute;
        
        // Print node indices (convert to 1-based)
        auto connectivity = it->getConnectivity();
        for (size_t i = 0; i < connectivity.size(); i++)
        {
          os << " " << (connectivity[i] + 1);
        }
        os << "\n";
        
        elementId++;
      }
    }
    
    os << "$EndElements\n";
  }
}