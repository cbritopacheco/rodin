/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_MFEM_H
#define RODIN_IO_MFEM_H

#include <iomanip>
#include <boost/bimap.hpp>
#include <boost/spirit/home/x3.hpp>

#include "Rodin/Types.h"
#include "Rodin/Alert.h"
#include "Rodin/Context.h"
#include "Rodin/Geometry/Types.h"
#include "Rodin/Geometry/Polytope.h"

#include "ForwardDecls.h"
#include "MeshLoader.h"
#include "MeshPrinter.h"
#include "GridFunctionLoader.h"
#include "GridFunctionPrinter.h"

#include "Rodin/Variational/P0/P0.h"
#include "Rodin/Variational/P1/P1.h"

namespace Rodin::IO::MFEM
{
  /**
   * @internal
   */
  std::istream& getline(std::istream& is, std::string& line, size_t& currentLineNumber);

  /**
   * @internal
   */
  std::string skipEmptyLinesAndComments(std::istream& is, size_t& currentLineNumber);

  /**
   * @internal
   */
  enum MeshType
  {
    LEGACY,
    NONCONFORMING,
    NURBS
  };

  /**
   * @internal
   */
  enum GeometryType
  {
    POINT       = 0,
    SEGMENT     = 1,
    TRIANGLE    = 2,
    SQUARE      = 3,
    TETRAHEDRON = 4,
    CUBE        = 5,
    PRISM       = 6,
    PYRAMID     = 7
  };

  struct MeshVersion
  {
    size_t major;
    size_t minor;
  };

  struct MeshHeader
  {
    MeshType type;
    MeshVersion version;
  };

  enum Ordering
  {
    /// XXX..., YYY..., ZZZ...
    Nodes = 0,

    /// XYZ, XYZ, ...
    VectorDimension = 1
  };

  struct GridFunctionHeader
  {
    std::string fec;
    size_t vdim;
    Ordering ordering;
  };

  enum class Keyword
  {
    dimension,
    elements,
    boundary,
    vertices
  };

  inline
  constexpr
  std::optional<Rodin::Geometry::Polytope::Type> getGeometry(GeometryType t)
  {
    switch (t)
    {
      case GeometryType::POINT:
      {
        return Rodin::Geometry::Polytope::Type::Point;
      }
      case GeometryType::SEGMENT:
      {
        return Rodin::Geometry::Polytope::Type::Segment;
      }
      case GeometryType::TRIANGLE:
      {
        return Rodin::Geometry::Polytope::Type::Triangle;
      }
      case GeometryType::TETRAHEDRON:
      {
        return Rodin::Geometry::Polytope::Type::Tetrahedron;
      }
      case GeometryType::SQUARE:
      {
        return Rodin::Geometry::Polytope::Type::Quadrilateral;
      }
      default:
        return {};
    }
    return {};
  }

  inline
  constexpr
  std::optional<GeometryType> getGeometry(Geometry::Polytope::Type t)
  {
    switch (t)
    {
      case Geometry::Polytope::Type::Point:
        return GeometryType::POINT;
      case Geometry::Polytope::Type::Segment:
        return GeometryType::SEGMENT;
      case Geometry::Polytope::Type::Triangle:
        return GeometryType::TRIANGLE;
      case Geometry::Polytope::Type::Quadrilateral:
        return GeometryType::SQUARE;
      case Geometry::Polytope::Type::Tetrahedron:
        return GeometryType::TETRAHEDRON;
      default:
        return {};
    }
    assert(false);
    return {};
  }

  inline
  constexpr
  const char* toCharString(Keyword kw)
  {
    switch (kw)
    {
      case Keyword::dimension:
        return "dimension";
      case Keyword::elements:
        return "elements";
      case Keyword::boundary:
        return "boundary";
      case Keyword::vertices:
        return "vertices";
    }
    return nullptr;
  }

  inline
  bool operator==(const std::string& str, Keyword kw)
  {
    return str == toCharString(kw);
  }

  inline
  bool operator==(Keyword kw, const std::string& str)
  {
    return str == toCharString(kw);
  }

  inline
  bool operator==(Keyword kw, const char* str)
  {
    return strcmp(toCharString(kw), str) == 0;
  }

  inline
  bool operator==(const char* str, Keyword kw)
  {
    return strcmp(toCharString(kw), str) == 0;
  }

  inline
  bool operator!=(const char* str, Keyword kw)
  {
    return !operator==(str, kw);
  }

  inline
  bool operator!=(const std::string& str, Keyword kw)
  {
    return !operator==(str, kw);
  }

  inline
  bool operator!=(Keyword kw, const std::string& str)
  {
    return !operator==(str, kw);
  }

  inline
  std::ostream& operator<<(std::ostream& os, Keyword kw)
  {
    os << toCharString(kw);
    return os;
  }

  inline
  std::optional<Keyword> toKeyword(const char* str)
  {
    Keyword res;
    if (str == Keyword::boundary)
      res = Keyword::boundary;
    else if (str == Keyword::dimension)
      res = Keyword::dimension;
    else if (str == Keyword::elements)
      res = Keyword::elements;
    else if (str == Keyword::vertices)
      res = Keyword::vertices;
    else
      return {};
    assert(res == str);
    return res;
  }

  struct ParseUnsignedInteger
  {
    template <class Iterator>
    inline
    std::optional<unsigned int> operator()(Iterator begin, Iterator end) const
    {
      using boost::spirit::x3::space;
      using boost::spirit::x3::blank;
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

  class ParseVertex
  {
    public:
      ParseVertex(size_t sdim)
        : m_sdim(sdim)
      {}

      template <class Iterator>
      inline
      std::optional<Math::Vector<Scalar>> operator()(Iterator begin, Iterator end) const
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::blank;
        using boost::spirit::x3::repeat;
        using boost::spirit::x3::double_;
        using boost::spirit::x3::_attr;

        size_t i = 0;
        Math::Vector<Scalar> res(m_sdim);
        const auto get_double = [&](auto& ctx) { assert(i < m_sdim); res(i++) = _attr(ctx); };
        const auto p = repeat(m_sdim)[double_[get_double]];
        const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
        if (begin != end)
          return {};
        else if (r)
          return res;
        else
          return {};
      }

    private:
      size_t m_sdim;
  };

  struct ParseGeometry
  {
    struct Data
    {
      Geometry::Attribute attribute;
      Geometry::Polytope::Type geometry;
      Array<Index> vertices;
    };

    template <class Iterator>
    inline
    std::optional<Data> operator()(Iterator begin, Iterator end) const
    {
      using boost::spirit::x3::space;
      using boost::spirit::x3::blank;
      using boost::spirit::x3::uint_;
      using boost::spirit::x3::_attr;
      using boost::spirit::x3::repeat;

      Data res;
      GeometryType geometry;
      const auto get_attribute = [&](auto& ctx) { res.attribute = _attr(ctx); };
      const auto get_geometry = [&](auto& ctx) { geometry = static_cast<GeometryType>(_attr(ctx)); };
      const auto p = uint_[get_attribute] >> uint_[get_geometry];
      const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);

      if (!r)
        return {};
      auto g = getGeometry(geometry);
      if (!g)
        return {};
      res.geometry = *g;

     res.vertices.resize(Geometry::Polytope::getVertexCount(res.geometry));
     size_t i = 0;
     const auto get_vertex = [&](auto& ctx) { res.vertices(i++) = _attr(ctx); };
     const auto pvs = repeat(res.vertices.size())[uint_[get_vertex]];
     const bool rvs = boost::spirit::x3::phrase_parse(begin, end, pvs, space);

      if (begin != end)
        return {};
      else if (rvs)
        return res;
      else
        return {};
    }
  };

  struct ParseEmptyLine
  {
    template <class Iterator>
    inline
    bool operator()(Iterator begin, Iterator end) const
    {
      if (begin == end)
        return true;
      using boost::spirit::x3::space;
      using boost::spirit::x3::blank;
      using boost::spirit::x3::char_;
      const auto p = *blank;
      const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
      if (begin != end)
        return false;
      return r;
    }
  };

  struct ParseEmptyLineOrComment
  {
    template <class Iterator>
    inline
    bool operator()(Iterator begin, Iterator end) const
    {
      if (begin == end)
        return true;
      using boost::spirit::x3::space;
      using boost::spirit::x3::blank;
      using boost::spirit::x3::char_;
      const auto comment = boost::spirit::x3::char_('#') >> *char_;
      const auto p = comment | *blank;
      const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
      if (begin != end)
        return false;
      return r;
    }
  };

  struct ParseKeyword
  {
    template <class Iterator>
    inline
    std::optional<std::string> operator()(Iterator begin, Iterator end) const
    {
      using boost::spirit::x3::space;
      using boost::spirit::x3::blank;
      using boost::spirit::x3::uint_;
      using boost::spirit::x3::_attr;
      using boost::spirit::x3::alpha;

      std::string kw;
      const auto get_keyword = [&](auto& ctx) { kw = _attr(ctx); };
      const auto p = (+alpha)[get_keyword];
      const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
      if (begin != end)
        return {};
      else if (r)
        return kw;
      else
        return {};
    }
  };

  class ParseMeshHeader
  {
    public:
      template <class Iterator>
      inline
      std::optional<MeshHeader> operator()(Iterator begin, Iterator end) const
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::blank;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::alpha;

        MeshHeader h;
        const auto get_major_version = [&](auto& ctx) { h.version.major = _attr(ctx); };
        const auto get_minor_version = [&](auto& ctx) { h.version.minor = _attr(ctx); };
        const auto p =
          boost::spirit::x3::string("MFEM")
            >> boost::spirit::x3::string("mesh")
            >> boost::spirit::x3::char_('v') >> uint_[get_major_version]
            >> boost::spirit::x3::char_('.') >> uint_[get_minor_version];
        const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
        h.type = MeshType::LEGACY;
        if (begin != end)
          return {};
        else if (r)
          return h;
        else
          return {};
      }
  };
}

namespace Rodin::IO
{
  template <>
  class MeshLoader<IO::FileFormat::MFEM, Context::Sequential>
    : public MeshLoaderBase<Context::Sequential>
  {
    public:
      using ContextType = Context::Sequential;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = MeshPrinterBase<ContextType>;

      MeshLoader(ObjectType& mesh)
        : MeshLoaderBase<Context::Sequential>(mesh)
      {}

      void load(std::istream& is) override;

      void readHeader(std::istream& is);

      void readDimension(std::istream& is);

      void readMesh(std::istream& is);

    private:
      size_t m_dimension;
      size_t m_spaceDimension;
      size_t m_currentLineNumber;
      MFEM::MeshHeader m_header;
      ObjectType::Builder m_build;
  };

  template <>
  class MeshPrinter<FileFormat::MFEM, Context::Sequential>
    : public MeshPrinterBase<Context::Sequential>
  {
    public:
      using ContextType = Context::Sequential;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = MeshPrinterBase<ContextType>;

      MeshPrinter(const ObjectType& mesh)
        : MeshPrinterBase(mesh)
      {}

      void print(std::ostream& os) override;

      void printHeader(std::ostream& os);
      void printDimension(std::ostream& os);
      void printMesh(std::ostream& os);
  };

  template <class Range>
  class GridFunctionPrinter<FileFormat::MFEM, Variational::P0<Range, Geometry::Mesh<Context::Sequential>>>
    : public GridFunctionPrinterBase<Variational::P0<Range, Geometry::Mesh<Context::Sequential>>>
  {
    public:
      using FESType = Variational::P0<Range, Geometry::Mesh<Context::Sequential>>;

      using ObjectType = Variational::GridFunction<FESType>;

      using Parent = GridFunctionPrinterBase<FESType>;

      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      void print(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        os << "FiniteElementSpace\n"
           << "FiniteElementCollection: " << "L2_" << fes.getMesh().getDimension() << "D_P0\n"
           << "VDim: " << fes.getVectorDimension() << '\n'
           << "Ordering: " << MFEM::Ordering::VectorDimension
           << "\n\n";
        const auto& matrix = gf.getData();
        const Scalar* data = matrix.data();
        assert(matrix.size() >= 0);
        for (size_t i = 0; i < static_cast<size_t>(matrix.size()); i++)
          os << data[i] << '\n';
      }
  };

  template <class Range>
  class GridFunctionLoader<FileFormat::MFEM, Variational::P1<Range, Geometry::Mesh<Context::Sequential>>>
    : public GridFunctionLoaderBase<Variational::P1<Range, Geometry::Mesh<Context::Sequential>>>
  {
    public:
      using FESType = Variational::P1<Range, Geometry::Mesh<Context::Sequential>>;

      using ObjectType = Variational::GridFunction<FESType>;

      using Parent = GridFunctionLoaderBase<FESType>;

      GridFunctionLoader(ObjectType& gf)
        : Parent(gf)
      {}

      void load(std::istream& is) override;

    private:
      size_t m_dimension;
      size_t m_spaceDimension;
      size_t m_currentLineNumber;
  };

  template <class Range>
  class GridFunctionPrinter<FileFormat::MFEM, Variational::P1<Range, Geometry::Mesh<Context::Sequential>>>
    : public GridFunctionPrinterBase<Variational::P1<Range, Geometry::Mesh<Context::Sequential>>>
  {
    public:
      using FESType = Variational::P1<Range, Geometry::Mesh<Context::Sequential>>;

      using ObjectType = Variational::GridFunction<FESType>;

      using Parent = GridFunctionPrinterBase<FESType>;

      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      void print(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        os << "FiniteElementSpace\n"
           << "FiniteElementCollection: " << "H1_" << fes.getMesh().getDimension() << "D_P1\n"
           << "VDim: " << fes.getVectorDimension() << '\n'
           << "Ordering: " << MFEM::Ordering::VectorDimension
           << "\n\n";
        const auto& matrix = gf.getData();
        const Scalar* data = matrix.data();
        assert(matrix.size() >= 0);
        for (size_t i = 0; i < static_cast<size_t>(matrix.size()); i++)
          os << data[i] << '\n';
      }
  };
}

#include "MFEM.hpp"

#endif
