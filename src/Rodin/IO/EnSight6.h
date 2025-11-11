/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_ENSIGHT6_H
#define RODIN_IO_ENSIGHT6_H

#include <cassert>
#include <iomanip>
#include <boost/bimap.hpp>
#include <boost/spirit/home/x3.hpp>

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"

#include "Rodin/Alert/MemberFunctionException.h"

#include "ForwardDecls.h"
#include "MeshPrinter.h"
#include "GridFunctionPrinter.h"
#include "Rodin/Utility/IsSpecialization.h"

namespace Rodin::IO::EnSight6
{
  /**
   * @brief Keywords used in EnSight6 file format.
   *
   * These keywords identify different sections and data organization in
   * EnSight6 geometry and result files. EnSight6 is a format used for
   * post-processing and visualization.
   *
   * @see <a href="https://www.ensight.com">EnSight</a>
   */
  enum class Keyword
  {
    node,         ///< Node (vertex) data
    id,           ///< Node/element ID specification
    off,          ///< ID mode off
    given,        ///< IDs are explicitly given
    assign,       ///< Assign IDs automatically
    ignore,       ///< Ignore ID information
    element,      ///< Element data
    coordinates,  ///< Coordinate data section
    part,         ///< Part definition
    block,        ///< Block data organization
    iblanked,     ///< Blanking information
    per           ///< Per-node or per-element data
  };

  /**
   * @brief Converts an EnSight6 keyword to its string representation.
   * @param[in] kw Keyword to convert
   * @returns C-style string representation
   */
  inline
  constexpr
  const char* toCharString(Keyword kw)
  {
    switch (kw)
    {
      case Keyword::node:
        return "node";
      case Keyword::id:
        return "id";
      case Keyword::off:
        return "off";
      case Keyword::given:
        return "given";
      case Keyword::assign:
        return "assign";
      case Keyword::ignore:
        return "ignore";
      case Keyword::element:
        return "element";
      case Keyword::coordinates:
        return "coordinates";
      case Keyword::part:
        return "part";
      case Keyword::block:
        return "block";
      case Keyword::iblanked:
        return "iblanked";
      case Keyword::per:
        return "per";
    }
    return nullptr;
  }

  /**
   * @brief Stream output operator for EnSight6 keywords.
   * @param[in,out] os Output stream
   * @param[in] kw Keyword to output
   * @returns Reference to the output stream
   */
  inline
  std::ostream& operator<<(std::ostream& os, Keyword kw)
  {
    os << toCharString(kw);
    return os;
  }

  /**
   * @brief Element types in EnSight6 format.
   *
   * Identifies the different element geometries supported by EnSight6.
   * The naming follows EnSight conventions where the number indicates
   * the node count (e.g., bar2 = 2-node line, tria3 = 3-node triangle).
   */
  enum class ElementType
  {
    point,      ///< Point (0D)
    bar2,       ///< 2-node line segment (linear)
    bar3,       ///< 3-node line segment (quadratic)
    tria3,      ///< 3-node triangle (linear)
    tria6,      ///< 6-node triangle (quadratic)
    quad4,      ///< 4-node quadrilateral (linear)
    quad8,      ///< 8-node quadrilateral (quadratic)
    tetra4,     ///< 4-node tetrahedron (linear)
    tetra10,    ///< 10-node tetrahedron (quadratic)
    pyramid5,   ///< 5-node pyramid (linear)
    pyramid13,  ///< 13-node pyramid (quadratic)
    hexa8,      ///< 8-node hexahedron (linear)
    hexa20,     ///< 20-node hexahedron (quadratic)
    penta6,     ///< 6-node pentahedron/wedge (linear)
    penta15,    ///< 15-node pentahedron/wedge (quadratic)
  };

  /**
   * @brief Converts an EnSight6 element type to its string representation.
   * @param[in] kw Element type to convert
   * @returns C-style string representation
   */
  inline
  constexpr
  const char* toCharString(ElementType kw)
  {
    switch (kw)
    {
      case ElementType::point:
        return "point";
      case ElementType::bar2:
        return "bar2";
      case ElementType::bar3:
        return "bar3";
      case ElementType::tria3:
        return "tria3";
      case ElementType::tria6:
        return "tria6";
      case ElementType::quad4:
        return "quad4";
      case ElementType::quad8:
        return "quad8";
      case ElementType::tetra4:
        return "tetra4";
      case ElementType::tetra10:
        return "tetra10";
      case ElementType::pyramid5:
        return "pyramid5";
      case ElementType::pyramid13:
        return "pyramid13";
      case ElementType::hexa8:
        return "hexa8";
      case ElementType::hexa20:
        return "hexa20";
      case ElementType::penta6:
        return "penta6";
      case ElementType::penta15:
        return "penta15";
    }
    return nullptr;
  }

  /**
   * @brief Stream output operator for EnSight6 element types.
   * @param[in,out] os Output stream
   * @param[in] kw Element type to output
   * @returns Reference to the output stream
   */
  inline
  std::ostream& operator<<(std::ostream& os, ElementType kw)
  {
    os << toCharString(kw);
    return os;
  }

  /**
   * @brief Converts Rodin polytope type to EnSight6 element type.
   * @param[in] t Rodin polytope type
   * @returns Optional EnSight6 element type, empty if not supported
   */
  inline
  constexpr
  Optional<ElementType> getGeometry(Geometry::Polytope::Type t)
  {
    switch (t)
    {
      case Geometry::Polytope::Type::Point:
        return ElementType::point;
      case Geometry::Polytope::Type::Segment:
        return ElementType::bar2;
      case Geometry::Polytope::Type::Triangle:
        return ElementType::tria3;
      case Geometry::Polytope::Type::Quadrilateral:
        return ElementType::quad4;
      case Geometry::Polytope::Type::Tetrahedron:
        return ElementType::tetra4;
      case Geometry::Polytope::Type::Wedge:
        return ElementType::penta6;
      default:
        return {};
    }
    assert(false);
    return {};
  }

  /**
   * @brief Variable types in EnSight6 result files.
   *
   * Identifies the type of solution variable being output.
   */
  enum class VariableType
  {
    scalar,   ///< Scalar field (e.g., temperature, pressure)
    complex,  ///< Complex-valued field
    vector    ///< Vector field (e.g., velocity, displacement)
  };

  /**
   * @brief Converts a variable type to its string representation.
   * @param[in] kw Variable type to convert
   * @returns C-style string representation
   */
  inline
  constexpr
  const char* toCharString(VariableType kw)
  {
    switch (kw)
    {
      case VariableType::scalar:
        return "scalar";
      case VariableType::complex:
        return "complex";
      case VariableType::vector:
        return "vector";
    }
    return nullptr;
  }

  /**
   * @brief Stream output operator for variable types.
   * @param[in,out] os Output stream
   * @param[in] kw Variable type to output
   * @returns Reference to the output stream
   */
  inline
  std::ostream& operator<<(std::ostream& os, VariableType kw)
  {
    os << toCharString(kw);
    return os;
  }

  /**
   * @brief Data location specification in EnSight6 format.
   *
   * Indicates whether variable data is defined at nodes or elements.
   */
  enum class Location
  {
    node,     ///< Data defined at nodes (vertices)
    element   ///< Data defined at elements (cells)
  };
}

namespace Rodin::IO
{
  /**
   * @ingroup PrinterSpecializations
   * @brief Specialization for printing sequential meshes in EnSight6 format.
   *
   * This printer writes mesh data in the EnSight6 geometry file format, which
   * is widely used for post-processing and visualization in engineering applications.
   *
   * ## EnSight6 Format
   * The format uses ASCII text and includes:
   * - Description line
   * - Node ID specification
   * - Element ID specification  
   * - Coordinate data (X, Y, Z components)
   * - Element connectivity by type
   *
   * ## Usage Example
   * ```cpp
   * const Mesh<Context::Local>& mesh = getMesh();
   * MeshPrinter<FileFormat::ENSIGHT6, Context::Local> printer(mesh);
   * std::ofstream file("mesh.geo");
   * printer.print(file);
   * ```
   *
   * @see MeshLoader
   */
  template <>
  class MeshPrinter<FileFormat::ENSIGHT6, Context::Local>
    : public MeshPrinterBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = MeshPrinterBase<ContextType>;

      MeshPrinter(const ObjectType& mesh);

      void print(std::ostream& os) override;

      void printHeader(std::ostream& os);
      void printCoordinates(std::ostream& os);
      void printParts(std::ostream& os);

    private:
      std::string m_descriptionLine1;
      std::string m_descriptionLine2;
  };

  template <class FES, class Data>
  class GridFunctionPrinterBase<FileFormat::ENSIGHT6, FES, Data>
    : public Printer<Variational::GridFunction<FES, Data>>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = FileFormat::ENSIGHT6;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using DataType = Data;

      using ObjectType = Variational::GridFunction<FESType, DataType>;

      using Parent = Printer<ObjectType>;

      GridFunctionPrinterBase(const ObjectType& gf)
        : m_gf(gf)
      {}

      void print(std::ostream& os) override
      {
        printHeader(os);
        this->printData(os);
      }

      void printHeader(std::ostream& os)
      {
        using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        if constexpr (std::is_same_v<RangeType, Real>)
        {
          os << EnSight6::VariableType::scalar << ' ';
        }
        else if constexpr (std::is_same_v<RangeType, Complex>)
        {
          os << EnSight6::VariableType::complex << ' '
             << EnSight6::VariableType::scalar << ' ';
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<Real>>)
        {
          os << EnSight6::VariableType::vector << ' ';
        }
        else if constexpr (std::is_same_v<RangeType, Math::Vector<Complex>>)
        {
          os << EnSight6::VariableType::complex << ' ' << EnSight6::VariableType::vector << ' ';
        }
        else
        {
          Alert::MemberFunctionException(*this, __func__)
            << "EnSight6 format does not support this RangeType."
            << Alert::Raise;
        }
        os << EnSight6::Keyword::per << ' ' << EnSight6::Keyword::node << '\n';
      }

      const ObjectType& getObject() const override
      {
        return m_gf.get();
      }

      virtual void printData(std::ostream& os) = 0;

    private:
      std::reference_wrapper<const ObjectType> m_gf;
  };

  template <class FES>
  class GridFunctionPrinter<
    FileFormat::ENSIGHT6, FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
    : public GridFunctionPrinterBase<
        FileFormat::ENSIGHT6, FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = FileFormat::ENSIGHT6;

      using RangeType = typename FormLanguage::Traits<FES>::RangeType;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using DataType = Math::Vector<ScalarType>;

      using ObjectType = Variational::GridFunction<FESType, DataType>;

      using Parent = GridFunctionPrinterBase<Format, FESType, DataType>;

      using Parent::Parent;

      void printData(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        os << std::setprecision(5) << std::scientific;
        size_t count = 0;
        using RangeType = typename FormLanguage::Traits<FESType>::RangeType;
        RangeType v;
        for (auto it = mesh.getVertex(); !it.end(); ++it)
        {
          const Geometry::Point p(
            *it,
            Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
            it->getCoordinates()
          );

          v = gf(p);

          if constexpr (std::is_same_v<RangeType, Real>)
          {
            os << std::setw(12) << v;
            if (++count % 6 == 0)
              os << '\n';
          }
          else if constexpr (Utility::IsSpecialization<RangeType, Math::Vector>::Value)
          {
            for (size_t j = 0; j < v.size(); ++j)
            {
              os << std::setw(12) << v[j];
              if (++count % 6 == 0)
                os << '\n';
            }
          }
          else
          {
            assert(false);
          }
        }

        // If we didn’t end exactly on a multiple of 6, finish the last line
        if (count % 6 != 0)
          os << '\n';
      }
  };

}

#endif
