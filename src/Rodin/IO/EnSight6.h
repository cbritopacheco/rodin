#ifndef RODIN_IO_ENSIGHT6_H
#define RODIN_IO_ENSIGHT6_H

#include <cassert>
#include <iomanip>
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
}

namespace Rodin::IO
{
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
