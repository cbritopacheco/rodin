#ifndef RODIN_IO_ENSIGHT6_H
#define RODIN_IO_ENSIGHT6_H

#include <cassert>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <boost/bimap.hpp>
#include <boost/spirit/home/x3.hpp>

#include "Rodin/Types.h"
#include "Rodin/Alert.h"
#include "Rodin/Context.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Geometry/Types.h"

#include "Rodin/Alert/MemberFunctionException.h"

#include "ForwardDecls.h"
#include "MeshLoader.h"
#include "MeshPrinter.h"
#include "GridFunctionLoader.h"
#include "GridFunctionPrinter.h"

namespace Rodin::IO::EnSight6
{
  enum class Keyword
  {
    node,
    id,
    off,
    given,
    assign,
    ignore,
    element,
    coordinates,
    part,
    block,
    iblanked,
    per
  };

  /**
   * @brief Parser for parsing empty lines or comments in EnSight6 files
   */
  struct ParseEmptyLine
  {
    template <class Iterator>
    inline
    bool operator()(Iterator begin, Iterator end) const
    {
      using boost::spirit::x3::space;
      using boost::spirit::x3::char_;
      const auto p = *space >> ('#' >> *(char_ - '\n') | boost::spirit::x3::eoi);
      return boost::spirit::x3::phrase_parse(begin, end, p, boost::spirit::x3::blank);
    }
  };

  /**
   * @brief Parser for parsing keywords in EnSight6 files
   */
  struct ParseKeyword
  {
    template <class Iterator>
    inline
    Optional<std::string> operator()(Iterator begin, Iterator end) const
    {
      using boost::spirit::x3::alpha;
      using boost::spirit::x3::alnum;
      using boost::spirit::x3::char_;
      std::string result;
      const auto p = +(alpha | char_('_')) >> *(alnum | char_('_'));
      const bool r = boost::spirit::x3::phrase_parse(begin, end, p, boost::spirit::x3::space, result);
      if (r)
        return result;
      return {};
    }
  };

  /**
   * @brief Parser for parsing unsigned integers in EnSight6 files
   */
  struct ParseUnsignedInteger
  {
    template <class Iterator>
    inline
    Optional<unsigned int> operator()(Iterator begin, Iterator end) const
    {
      using boost::spirit::x3::uint_;
      unsigned int result;
      const bool r = boost::spirit::x3::phrase_parse(begin, end, uint_, boost::spirit::x3::space, result);
      if (r)
        return result;
      return {};
    }
  };

  /**
   * @brief Parser for parsing floating point numbers in EnSight6 files
   */
  struct ParseReal
  {
    template <class Iterator>
    inline
    Optional<Real> operator()(Iterator begin, Iterator end) const
    {
      using boost::spirit::x3::double_;
      Real result;
      const bool r = boost::spirit::x3::phrase_parse(begin, end, double_, boost::spirit::x3::space, result);
      if (r)
        return result;
      return {};
    }
  };

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

  inline
  std::ostream& operator<<(std::ostream& os, Keyword kw)
  {
    os << toCharString(kw);
    return os;
  }

  enum class ElementType
  {
    point,
    bar2,
    bar3,
    tria3,
    tria6,
    quad4,
    quad8,
    tetra4,
    tetra10,
    pyramid5,
    pyramid13,
    hexa8,
    hexa20,
    penta6,
    penta15,
  };

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

  inline
  std::ostream& operator<<(std::ostream& os, ElementType kw)
  {
    os << toCharString(kw);
    return os;
  }

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

  enum class VariableType
  {
    scalar,
    complex,
    vector
  };

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

  inline
  std::ostream& operator<<(std::ostream& os, VariableType kw)
  {
    os << toCharString(kw);
    return os;
  }

  enum class Location
  {
    node,
    element
  };

  /**
   * @brief Convert string to Keyword enum
   */
  inline
  Optional<Keyword> toKeyword(const char* str)
  {
    if (std::strcmp(str, "node") == 0)
      return Keyword::node;
    else if (std::strcmp(str, "id") == 0)
      return Keyword::id;
    else if (std::strcmp(str, "off") == 0)
      return Keyword::off;
    else if (std::strcmp(str, "given") == 0)
      return Keyword::given;
    else if (std::strcmp(str, "assign") == 0)
      return Keyword::assign;
    else if (std::strcmp(str, "ignore") == 0)
      return Keyword::ignore;
    else if (std::strcmp(str, "element") == 0)
      return Keyword::element;
    else if (std::strcmp(str, "coordinates") == 0)
      return Keyword::coordinates;
    else if (std::strcmp(str, "part") == 0)
      return Keyword::part;
    else if (std::strcmp(str, "block") == 0)
      return Keyword::block;
    else if (std::strcmp(str, "iblanked") == 0)
      return Keyword::iblanked;
    else if (std::strcmp(str, "per") == 0)
      return Keyword::per;
    return {};
  }

  /**
   * @brief Convert string to ElementType enum
   */
  inline
  Optional<ElementType> toElementType(const char* str)
  {
    if (std::strcmp(str, "point") == 0)
      return ElementType::point;
    else if (std::strcmp(str, "bar2") == 0)
      return ElementType::bar2;
    else if (std::strcmp(str, "bar3") == 0)
      return ElementType::bar3;
    else if (std::strcmp(str, "tria3") == 0)
      return ElementType::tria3;
    else if (std::strcmp(str, "tria6") == 0)
      return ElementType::tria6;
    else if (std::strcmp(str, "quad4") == 0)
      return ElementType::quad4;
    else if (std::strcmp(str, "quad8") == 0)
      return ElementType::quad8;
    else if (std::strcmp(str, "tetra4") == 0)
      return ElementType::tetra4;
    else if (std::strcmp(str, "tetra10") == 0)
      return ElementType::tetra10;
    else if (std::strcmp(str, "pyramid5") == 0)
      return ElementType::pyramid5;
    else if (std::strcmp(str, "pyramid13") == 0)
      return ElementType::pyramid13;
    else if (std::strcmp(str, "hexa8") == 0)
      return ElementType::hexa8;
    else if (std::strcmp(str, "hexa20") == 0)
      return ElementType::hexa20;
    else if (std::strcmp(str, "penta6") == 0)
      return ElementType::penta6;
    else if (std::strcmp(str, "penta15") == 0)
      return ElementType::penta15;
    return {};
  }

  /**
   * @brief Convert ElementType to Geometry::Polytope::Type
   */
  inline
  Optional<Geometry::Polytope::Type> toPolytopeType(ElementType et)
  {
    switch (et)
    {
      case ElementType::point:
        return Geometry::Polytope::Type::Point;
      case ElementType::bar2:
      case ElementType::bar3:
        return Geometry::Polytope::Type::Segment;
      case ElementType::tria3:
      case ElementType::tria6:
        return Geometry::Polytope::Type::Triangle;
      case ElementType::quad4:
      case ElementType::quad8:
        return Geometry::Polytope::Type::Quadrilateral;
      case ElementType::tetra4:
      case ElementType::tetra10:
        return Geometry::Polytope::Type::Tetrahedron;
      case ElementType::penta6:
      case ElementType::penta15:
        return Geometry::Polytope::Type::Wedge;
      case ElementType::pyramid5:
      case ElementType::pyramid13:
        // Return Tetrahedron for now as Pyramid might not be supported
        return Geometry::Polytope::Type::Tetrahedron;
      case ElementType::hexa8:
      case ElementType::hexa20:
        // Return Quadrilateral for now as Hexahedron might not be supported  
        return Geometry::Polytope::Type::Quadrilateral;
      default:
        return {};
    }
  }
}

namespace Rodin::IO
{
  template <>
  class MeshLoader<FileFormat::ENSIGHT6, Context::Local>
    : public MeshLoaderBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = MeshLoaderBase<ContextType>;

      MeshLoader(ObjectType& mesh);

      void load(std::istream& is) override;

      std::istream& getline(std::istream& is, std::string& line);
      std::string skipEmptyLines(std::istream& is);
      void readHeader(std::istream& is);
      void readCoordinates(std::istream& is);
      void readParts(std::istream& is);

    private:
      Geometry::Mesh<Context::Local>::Builder m_build;
      size_t m_currentLineNumber;
      size_t m_spaceDimension;
      std::string m_descriptionLine1;
      std::string m_descriptionLine2;
  };

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
  class GridFunctionLoaderBase<FileFormat::ENSIGHT6, FES, Data>
    : public Loader<Variational::GridFunction<FES, Data>>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = FileFormat::ENSIGHT6;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using DataType = Data;

      using ObjectType = Variational::GridFunction<FESType, DataType>;

      using Parent = Loader<ObjectType>;

      GridFunctionLoaderBase(ObjectType& gf)
        : m_gf(gf),
          m_currentLineNumber(0)
      {}

      void load(std::istream& is) override
      {
        readHeader(is);
        this->readData(is);
      }

      void readHeader(std::istream& is)
      {
        auto line = skipEmptyLines(is);
        // Parse variable type (scalar, complex, vector)
        // For now, we'll implement a basic parser that expects the variable type
      }

      const ObjectType& getObject() const override
      {
        return m_gf.get();
      }

      ObjectType& getObject() override
      {
        return m_gf.get();
      }

      std::istream& getline(std::istream& is, std::string& line)
      {
        m_currentLineNumber++;
        return std::getline(is, line);
      }

      std::string skipEmptyLines(std::istream& is)
      {
        std::string line;
        while (getline(is, line))
        {
          if (!EnSight6::ParseEmptyLine()(line.begin(), line.end()))
            break;
        }
        return line;
      }

      virtual void readData(std::istream& is) = 0;

    private:
      std::reference_wrapper<ObjectType> m_gf;
      size_t m_currentLineNumber;
  };

  template <class FES>
  class GridFunctionLoader<
    FileFormat::ENSIGHT6, FES, typename FormLanguage::Traits<FES>::ScalarType>
    : public GridFunctionLoaderBase<
        FileFormat::ENSIGHT6, FES, typename FormLanguage::Traits<FES>::ScalarType>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = FileFormat::ENSIGHT6;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using DataType = ScalarType;

      using ObjectType = Variational::GridFunction<FESType, DataType>;

      using Parent = GridFunctionLoaderBase<Format, FESType, DataType>;

      using Parent::Parent;

      void readData(std::istream& is) override
      {
        auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        if constexpr (Utility::IsSpecialization<FES, Variational::P1>::Value)
        {
          auto& data = gf.getData();
          for (size_t i = 0; i < data.size(); i++)
          {
            Real value;
            is >> value;
            data[i] = value;
          }
        }
        else
        {
          // For other finite element spaces, read nodal values
          size_t nodeIndex = 0;
          for (auto it = mesh.getVertex(); !it.end(); ++it, ++nodeIndex)
          {
            Real value;
            is >> value;
            if (nodeIndex < gf.getData().size())
              gf.getData()[nodeIndex] = value;
          }
        }
      }
  };

  template <class FES>
  class GridFunctionLoader<
    FileFormat::ENSIGHT6, FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
    : public GridFunctionLoaderBase<
        FileFormat::ENSIGHT6, FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = FileFormat::ENSIGHT6;

      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using DataType = Math::Vector<ScalarType>;

      using ObjectType = Variational::GridFunction<FESType, DataType>;

      using Parent = GridFunctionLoaderBase<Format, FESType, DataType>;

      using Parent::Parent;

      void readData(std::istream& is) override
      {
        auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t vdim = fes.getVectorDimension();

        if constexpr (Utility::IsSpecialization<FES, Variational::P1>::Value)
        {
          auto& data = gf.getData();
          for (size_t i = 0; i < data.size(); i += 3)
          {
            Real x0 = 0.0, x1 = 0.0, x2 = 0.0;
            is >> x0 >> x1 >> x2;
            if (i < data.size()) data[i] = x0;
            if (i + 1 < data.size()) data[i + 1] = x1;
            if (i + 2 < data.size()) data[i + 2] = x2;
          }
        }
        else
        {
          // For other finite element spaces, read nodal values
          size_t nodeIndex = 0;
          for (auto it = mesh.getVertex(); !it.end(); ++it, ++nodeIndex)
          {
            const Geometry::Point p(
              *it,
              Geometry::Polytope::Traits(Geometry::Polytope::Type::Point).getVertex(0),
              it->getCoordinates()
            );

            RangeType v;
            if constexpr (std::is_same_v<RangeType, Real>)
            {
              Real value;
              is >> value;
              gf.getData()[nodeIndex] = value;
            }
            else if constexpr (Utility::IsSpecialization<RangeType, Math::Vector>::Value)
            {
              for (size_t j = 0; j < v.size(); ++j)
              {
                Real value;
                is >> value;
                v[j] = value;
              }
              // Set the grid function value at this point
              // Note: This is a simplified approach for demonstration
            }
          }
        }
      }
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
    FileFormat::ENSIGHT6, FES, typename FormLanguage::Traits<FES>::ScalarType>
    : public GridFunctionPrinterBase<
        FileFormat::ENSIGHT6, FES, typename FormLanguage::Traits<FES>::ScalarType>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = FileFormat::ENSIGHT6;

      using RangeType = typename FormLanguage::Traits<FES>::RangeType;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using DataType = ScalarType;

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
        
        if constexpr (Utility::IsSpecialization<FES, Variational::P1>::Value)
        {
          const auto& data = gf.getData();
          for (size_t i = 0; i < data.size(); i++)
          {
            os << std::setw(12) << data[i];
            if (++count % 6 == 0)
              os << '\n';
          }
        }
        else
        {
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
            os << std::setw(12) << v;
            if (++count % 6 == 0)
              os << '\n';
          }
        }

        // If we didn't end exactly on a multiple of 6, finish the last line
        if (count % 6 != 0)
          os << '\n';
      }
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
        const size_t vdim = fes.getVectorDimension();
        os << std::setprecision(5) << std::scientific;
        size_t count = 0;
        if constexpr (Utility::IsSpecialization<FES, Variational::P1>::Value)
        {
          const auto& data = gf.getData();
          while (count < fes.getSize())
          {
            Real x0 = 0.0, x1 = 0.0, x2 = 0.0;
            if (vdim > 0) x0 = data[count];
            if (vdim > 1) x1 = data[count + 1];
            if (vdim > 2) x2 = data[count + 2];

            // Always write three components: X, Y, Z
            os << std::setw(12) << x0
               << std::setw(12) << x1
               << std::setw(12) << x2;

            count += 3;
            // os << std::setw(12) << data[i];
            if (count % 6 == 0)
              os << '\n';
          }
        }
        else
        {
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
        }

        // If we didn’t end exactly on a multiple of 6, finish the last line
        if (count % 6 != 0)
          os << '\n';
      }
  };

}

#endif
