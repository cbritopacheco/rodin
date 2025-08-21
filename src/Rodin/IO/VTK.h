/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_VTK_H
#define RODIN_IO_VTK_H

#include <iomanip>
#include <boost/spirit/home/x3.hpp>

#include "Rodin/Types.h"
#include "Rodin/Alert.h"
#include "Rodin/Context.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Geometry/Types.h"
#include "Rodin/FormLanguage/Traits.h"

#include "ForwardDecls.h"
#include "MeshLoader.h"
#include "MeshPrinter.h"

namespace Rodin::IO::VTK
{
  enum class Keyword
  {
    DATASET,
    POINTS,
    CELLS,
    CELL_TYPES,
    POINT_DATA,
    CELL_DATA,
    SCALARS,
    VECTORS,
    UNSTRUCTURED_GRID
  };

  inline
  constexpr
  const char* toCharString(Keyword kw)
  {
    switch (kw)
    {
      case Keyword::DATASET:
        return "DATASET";
      case Keyword::POINTS:
        return "POINTS";
      case Keyword::CELLS:
        return "CELLS";
      case Keyword::CELL_TYPES:
        return "CELL_TYPES";
      case Keyword::POINT_DATA:
        return "POINT_DATA";
      case Keyword::CELL_DATA:
        return "CELL_DATA";
      case Keyword::SCALARS:
        return "SCALARS";
      case Keyword::VECTORS:
        return "VECTORS";
      case Keyword::UNSTRUCTURED_GRID:
        return "UNSTRUCTURED_GRID";
    }
    return nullptr;
  }

  inline
  std::ostream& operator<<(std::ostream& os, Keyword kw)
  {
    os << toCharString(kw);
    return os;
  }

  /**
   * @brief VTK cell types
   */
  enum CellType
  {
    VTK_VERTEX        = 1,
    VTK_LINE          = 3,
    VTK_TRIANGLE      = 5,
    VTK_QUAD          = 9,
    VTK_TETRA         = 10,
    VTK_HEXAHEDRON    = 12,
    VTK_WEDGE         = 13,
    VTK_PYRAMID       = 14
  };

  inline
  constexpr
  Optional<CellType> getCellType(Geometry::Polytope::Type type)
  {
    switch (type)
    {
      case Geometry::Polytope::Type::Point:
        return VTK_VERTEX;
      case Geometry::Polytope::Type::Segment:
        return VTK_LINE;
      case Geometry::Polytope::Type::Triangle:
        return VTK_TRIANGLE;
      case Geometry::Polytope::Type::Quadrilateral:
        return VTK_QUAD;
      case Geometry::Polytope::Type::Tetrahedron:
        return VTK_TETRA;
      case Geometry::Polytope::Type::Wedge:
        return VTK_WEDGE;
    }
    return {};
  }

  inline
  constexpr
  Optional<Geometry::Polytope::Type> getGeometry(CellType type)
  {
    switch (type)
    {
      case VTK_VERTEX:
        return Geometry::Polytope::Type::Point;
      case VTK_LINE:
        return Geometry::Polytope::Type::Segment;
      case VTK_TRIANGLE:
        return Geometry::Polytope::Type::Triangle;
      case VTK_QUAD:
        return Geometry::Polytope::Type::Quadrilateral;
      case VTK_TETRA:
        return Geometry::Polytope::Type::Tetrahedron;
      case VTK_WEDGE:
        return Geometry::Polytope::Type::Wedge;
      case VTK_HEXAHEDRON:
      case VTK_PYRAMID:
        // These types are not currently supported by Rodin
        return {};
    }
    return {};
  }

  /**
   * @internal
   */
  std::istream& getline(std::istream& is, std::string& line, size_t& currentLineNumber);

  /**
   * @internal
   */
  std::string skipEmptyLinesAndComments(std::istream& is, size_t& currentLineNumber);

  struct ParseUnsignedInteger
  {
    template <class Iterator>
    inline
    Optional<unsigned int> operator()(Iterator begin, Iterator end) const
    {
      using boost::spirit::x3::space;
      using boost::spirit::x3::uint_;
      using boost::spirit::x3::_attr;

      unsigned int v;
      const auto get_unsigned_integer = [&](auto& ctx) { v = _attr(ctx); };
      const auto p = uint_[get_unsigned_integer];
      const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
      if (begin != end)
        return {};
      else if (r)
        return v;
      else
        return {};
    }
  };

  struct ParseDouble
  {
    template <class Iterator>
    inline
    Optional<double> operator()(Iterator begin, Iterator end) const
    {
      using boost::spirit::x3::space;
      using boost::spirit::x3::double_;
      using boost::spirit::x3::_attr;

      double v;
      const auto get_double = [&](auto& ctx) { v = _attr(ctx); };
      const auto p = double_[get_double];
      const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
      if (begin != end)
        return {};
      else if (r)
        return v;
      else
        return {};
    }
  };

  class ParseVertex
  {
    public:
      struct Data
      {
        Math::SpatialPoint vertex;
      };

      constexpr
      ParseVertex(size_t sdim)
        : m_sdim(sdim)
      {}

      template <class Iterator>
      inline
      Optional<Data> operator()(Iterator begin, Iterator end) const
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::double_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::repeat;
        
        size_t i = 0;
        Data res{ Math::SpatialPoint(m_sdim) };
        const auto get_x = [&](auto& ctx) { assert(i < m_sdim); res.vertex(i++) = _attr(ctx); };
        const auto p = repeat(m_sdim)[double_[get_x]];
        const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
        if (begin != end)
          return {};
        else if (r)
          return { res };
        else
          return {};
      }

    private:
      const size_t m_sdim;
  };

  class ParseCell
  {
    public:
      struct Data
      {
        Array<Index> vertices;
        CellType cellType;
      };

      template <class Iterator>
      inline
      Optional<Data> operator()(Iterator begin, Iterator end) const
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::repeat;

        Data res;
        unsigned int numVertices;
        const auto get_num_vertices = [&](auto& ctx) { numVertices = _attr(ctx); };
        const auto p1 = uint_[get_num_vertices];
        bool r1 = boost::spirit::x3::phrase_parse(begin, end, p1, space);
        
        if (!r1)
          return {};

        res.vertices.resize(numVertices);
        size_t i = 0;
        const auto get_vertex = [&](auto& ctx) { res.vertices(i++) = _attr(ctx); };
        const auto p2 = repeat(numVertices)[uint_[get_vertex]];
        bool r2 = boost::spirit::x3::phrase_parse(begin, end, p2, space);

        if (begin != end)
          return {};
        else if (r2)
          return { res };
        else
          return {};
      }
  };
}

namespace Rodin::IO
{
  /**
   * @brief Template specialization to load VTKLegacy format meshes.
   * @ingroup MeshLoaderSpecializations
   */
  template <>
  class MeshLoader<FileFormat::VTKLegacy, Context::Local>
    : public MeshLoaderBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = MeshLoaderBase<ContextType>;

      MeshLoader(ObjectType& mesh)
        : Parent(mesh),
          m_currentLineNumber(0)
      {}

      void load(std::istream& is) override;

    private:
      void readHeader(std::istream& is);
      void readDataset(std::istream& is);
      void readPoints(std::istream& is);
      void readCells(std::istream& is);
      void readCellTypes(std::istream& is);

      Rodin::Geometry::Mesh<Rodin::Context::Local>::Builder m_build;
      
      size_t m_currentLineNumber;
      size_t m_spaceDimension;
      size_t m_numPoints;
      size_t m_numCells;
      std::vector<VTK::ParseCell::Data> m_cells;
  };

  /**
   * @brief Template specialization to print VTKLegacy format meshes.
   * @ingroup MeshPrinterSpecializations
   */
  template <>
  class MeshPrinter<FileFormat::VTKLegacy, Context::Local>
    : public MeshPrinterBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = MeshPrinterBase<ContextType>;

      MeshPrinter(const ObjectType& mesh)
        : Parent(mesh)
      {}

      void print(std::ostream& os) override;

    private:
      void printHeader(std::ostream& os);
      void printDataset(std::ostream& os);
      void printPoints(std::ostream& os);
      void printCells(std::ostream& os);
      void printCellTypes(std::ostream& os);
  };
}

#endif