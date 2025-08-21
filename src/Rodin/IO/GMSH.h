/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_GMSH_H
#define RODIN_IO_GMSH_H

#include <boost/bimap.hpp>
#include <boost/spirit/home/x3.hpp>

#include "Rodin/Types.h"
#include "Rodin/Alert.h"
#include "Rodin/Context.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Geometry/Types.h"

#include "ForwardDecls.h"
#include "MeshLoader.h"
#include "MeshPrinter.h"

namespace Rodin::IO::GMSH
{
  /**
   * @brief Section keywords in GMSH mesh file format version 2.
   */
  enum class Section
  {
    MeshFormat,
    Nodes,
    Elements
  };

  /**
   * @brief Element types in GMSH format.
   */
  enum class ElementType
  {
    Line = 1,
    Triangle = 2,
    Quadrangle = 3,
    Tetrahedron = 4,
    Hexahedron = 5,
    Prism = 6,
    Pyramid = 7,
    Line3 = 8,
    Triangle6 = 9,
    Quadrangle9 = 10,
    Tetrahedron10 = 11,
    Hexahedron27 = 12,
    Prism18 = 13,
    Pyramid14 = 14,
    Point = 15
  };

  /**
   * @brief GMSH mesh header structure.
   */
  struct MeshHeader
  {
    struct Version
    {
      int major;
      int minor;
    } version;
    int fileType;
    int dataSize;
  };

  /**
   * @brief Converts a section enum to its string representation.
   */
  inline
  constexpr
  const char* toCharString(Section section)
  {
    switch (section)
    {
      case Section::MeshFormat:
        return "$MeshFormat";
      case Section::Nodes:
        return "$Nodes";
      case Section::Elements:
        return "$Elements";
    }
    return nullptr;
  }

  /**
   * @brief Stream output operator for GMSH sections.
   */
  inline
  std::ostream& operator<<(std::ostream& os, Section section)
  {
    os << toCharString(section);
    return os;
  }

  /**
   * @brief Converts element type to Geometry::Polytope::Type.
   */
  inline
  Geometry::Polytope::Type toPolytopeType(ElementType type)
  {
    switch (type)
    {
      case ElementType::Point:
        return Geometry::Polytope::Type::Point;
      case ElementType::Line:
      case ElementType::Line3:
        return Geometry::Polytope::Type::Segment;
      case ElementType::Triangle:
      case ElementType::Triangle6:
        return Geometry::Polytope::Type::Triangle;
      case ElementType::Quadrangle:
      case ElementType::Quadrangle9:
        return Geometry::Polytope::Type::Quadrilateral;
      case ElementType::Tetrahedron:
      case ElementType::Tetrahedron10:
        return Geometry::Polytope::Type::Tetrahedron;
      case ElementType::Hexahedron:
      case ElementType::Hexahedron27:
        return Geometry::Polytope::Type::Hexahedron;
      case ElementType::Prism:
      case ElementType::Prism18:
        return Geometry::Polytope::Type::Prism;
      case ElementType::Pyramid:
      case ElementType::Pyramid14:
        return Geometry::Polytope::Type::Pyramid;
    }
    return Geometry::Polytope::Type::Point; // fallback
  }

  /**
   * @brief Converts Geometry::Polytope::Type to GMSH element type.
   */
  inline
  ElementType fromPolytopeType(Geometry::Polytope::Type type)
  {
    switch (type)
    {
      case Geometry::Polytope::Type::Point:
        return ElementType::Point;
      case Geometry::Polytope::Type::Segment:
        return ElementType::Line;
      case Geometry::Polytope::Type::Triangle:
        return ElementType::Triangle;
      case Geometry::Polytope::Type::Quadrilateral:
        return ElementType::Quadrangle;
      case Geometry::Polytope::Type::Tetrahedron:
        return ElementType::Tetrahedron;
      case Geometry::Polytope::Type::Hexahedron:
        return ElementType::Hexahedron;
      case Geometry::Polytope::Type::Prism:
        return ElementType::Prism;
      case Geometry::Polytope::Type::Pyramid:
        return ElementType::Pyramid;
    }
    return ElementType::Point; // fallback
  }

  /**
   * @internal
   */
  std::istream& getline(std::istream& is, std::string& line, size_t& currentLineNumber);

  /**
   * @internal
   */
  std::string skipEmptyLines(std::istream& is, size_t& currentLineNumber);
}

namespace Rodin::IO
{
  /**
   * @ingroup MeshLoaderSpecializations
   * @brief Specialization for loading Sequential meshes in the GMSH file format.
   */
  template <>
  class MeshLoader<IO::FileFormat::GMSH, Context::Local>
    : public MeshLoaderBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;
      using ObjectType = Geometry::Mesh<ContextType>;
      using Parent = MeshLoaderBase<ContextType>;

      MeshLoader(ObjectType& mesh)
        : MeshLoaderBase<Context::Local>(mesh)
      {}

      void load(std::istream& is) override;

      void readMeshFormat(std::istream& is);
      void readNodes(std::istream& is);
      void readElements(std::istream& is);

    private:
      size_t m_currentLineNumber;
      GMSH::MeshHeader m_header;
      ObjectType::Builder m_builder;
  };

  /**
   * @brief Specialization for printing Sequential meshes in the GMSH file format.
   */
  template <>
  class MeshPrinter<FileFormat::GMSH, Context::Local>
    : public MeshPrinterBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;
      using ObjectType = Geometry::Mesh<ContextType>;
      using Parent = MeshPrinterBase<ContextType>;

      MeshPrinter(const ObjectType& mesh)
        : MeshPrinterBase(mesh)
      {}

      void print(std::ostream& os) override;

      void printMeshFormat(std::ostream& os);
      void printNodes(std::ostream& os);
      void printElements(std::ostream& os);
  };
}

#endif