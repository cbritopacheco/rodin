/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_MEDIT_H
#define RODIN_IO_MEDIT_H

#include <iomanip>
#include <boost/bimap.hpp>
#include <boost/spirit/home/x3.hpp>

#include "Rodin/Types.h"
#include "Rodin/Alert.h"
#include "Rodin/Context.h"
#include "Rodin/Geometry/Types.h"

#include "ForwardDecls.h"
#include "MeshLoader.h"
#include "MeshPrinter.h"
#include "GridFunctionLoader.h"
#include "GridFunctionPrinter.h"

namespace Rodin::IO::MEDIT
{
  enum SolutionType
  {
    Scalar = 1,
    Vector = 2,
    Tensor = 3
  };

  enum class Keyword
  {
    MeshVersionFormatted,
    Dimension,
    Vertices,
    Triangles,
    Quadrilaterals,
    Tetrahedra,
    Corners,
    Ridges,
    Edges,
    SolAtVertices,
    SolAtEdges,
    SolAtTriangles,
    SolAtQuadrilaterals,
    SolAtTetrahedra,
    SolAtPentahedra,
    SolAtHexahedra,
    RequiredVertices,
    Normals,
    NormalAtVertices,
    Tangents,
    TangentAtVertices,
    End
  };

  inline
  constexpr
  const char* toCharString(Keyword kw)
  {
    switch (kw)
    {
      case Keyword::MeshVersionFormatted:
        return "MeshVersionFormatted";
      case Keyword::Dimension:
        return "Dimension";
      case Keyword::Vertices:
        return "Vertices";
      case Keyword::Triangles:
        return "Triangles";
      case Keyword::Quadrilaterals:
        return "Quadrilaterals";
      case Keyword::Tetrahedra:
        return "Tetrahedra";
      case Keyword::Corners:
        return "Corners";
      case Keyword::Ridges:
        return "Ridges";
      case Keyword::Edges:
        return "Edges";
      case Keyword::SolAtVertices:
        return "SolAtVertices";
      case Keyword::SolAtEdges:
        return "SolAtEdges";
      case Keyword::SolAtTriangles:
        return "SolAtTriangles";
      case Keyword::SolAtQuadrilaterals:
        return "SolAtQuadrilaterals";
      case Keyword::SolAtTetrahedra:
        return "SolAtTetrahedra";
      case Keyword::SolAtPentahedra:
        return "SolAtPentahedra";
      case Keyword::SolAtHexahedra:
        return "SolAtHexahedra";
      case Keyword::RequiredVertices:
        return "RequiredVertices";
      case Keyword::Normals:
        return "Normals";
      case Keyword::NormalAtVertices:
        return "NormalAtVertices";
      case Keyword::Tangents:
        return "Tangents";
      case Keyword::TangentAtVertices:
        return "TangentAtVertices";
      case Keyword::End:
        return "End";
    }
    return nullptr;
  }

  inline
  bool operator==(const std::string& str, Keyword kw)
  {
    return str == toCharString(kw);
  }

  inline
  bool operator!=(const std::string& str, Keyword kw)
  {
    return str != toCharString(kw);
  }

  inline
  bool operator==(Keyword kw, const std::string& str)
  {
    return str == toCharString(kw);
  }

  inline
  bool operator!=(Keyword kw, const std::string& str)
  {
    return str != toCharString(kw);
  }

  inline
  bool operator==(Keyword kw, const char* str)
  {
    return strcmp(toCharString(kw), str) == 0;
  }

  inline
  bool operator!=(Keyword kw, const char* str)
  {
    return strcmp(toCharString(kw), str) != 0;
  }

  inline
  bool operator==(const char* str, Keyword kw)
  {
    return strcmp(toCharString(kw), str) == 0;
  }

  inline
  bool operator!=(const char* str, Keyword kw)
  {
    return strcmp(toCharString(kw), str) != 0;
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
    if (str == Keyword::MeshVersionFormatted)
      res = Keyword::MeshVersionFormatted;
    else if (str == Keyword::Dimension)
      res = Keyword::Dimension;
    else if (str == Keyword::Vertices)
      res = Keyword::Vertices;
    else if (str == Keyword::Triangles)
      res = Keyword::Triangles;
    else if (str == Keyword::Quadrilaterals)
      res = Keyword::Quadrilaterals;
    else if (str == Keyword::Tetrahedra)
      res = Keyword::Tetrahedra;
    else if (str == Keyword::Corners)
      res = Keyword::Corners;
    else if (str == Keyword::Ridges)
      res = Keyword::Ridges;
    else if (str == Keyword::Edges)
      res = Keyword::Edges;
    else if (str == Keyword::SolAtVertices)
      res = Keyword::SolAtVertices;
    else if (str == Keyword::SolAtEdges)
      res = Keyword::SolAtEdges;
    else if (str == Keyword::SolAtTriangles)
      res = Keyword::SolAtTriangles;
    else if (str == Keyword::SolAtQuadrilaterals)
      res = Keyword::SolAtQuadrilaterals;
    else if (str == Keyword::SolAtTetrahedra)
      res = Keyword::SolAtTetrahedra;
    else if (str == Keyword::SolAtPentahedra)
      res = Keyword::SolAtPentahedra;
    else if (str == Keyword::SolAtHexahedra)
      res = Keyword::SolAtHexahedra;
    else if (str == Keyword::RequiredVertices)
      res = Keyword::RequiredVertices;
    else if (str == Keyword::Normals)
      res = Keyword::Normals;
    else if (str == Keyword::NormalAtVertices)
      res = Keyword::NormalAtVertices;
    else if (str == Keyword::Tangents)
      res = Keyword::Tangents;
    else if (str == Keyword::TangentAtVertices)
      res = Keyword::TangentAtVertices;
    else if (str == Keyword::End)
      res = Keyword::End;
    else
      return {};
    assert(res == str);
    return res;
  }

  class ParseEntity
  {
    public:
      struct Data
      {
        Array<Index> vertices;
        Geometry::Attribute attribute;
      };

      constexpr
      ParseEntity(size_t n)
        : m_n(n)
      {}

      template <class Iterator>
      inline
      std::optional<Data> operator()(Iterator begin, Iterator end) const
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::repeat;
        size_t i = 0;
        Data res{ Array<Index>(m_n), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE };
        const auto get_vertex = [&](auto& ctx) { assert(i < m_n); res.vertices(i++) = _attr(ctx); };
        const auto get_attribute = [&](auto& ctx) { res.attribute = _attr(ctx); };
        const auto p = uint_[get_vertex] >> repeat(m_n - 1)[uint_[get_vertex]] >> uint_[get_attribute];
        const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
        if (begin != end)
          return {};
        else if (r)
          return res;
        else
          return {};
      }

    private:
      size_t m_n;
  };

  class ParseVertex
  {
    public:
      struct Data
      {
        Math::Vector vertex;
        Geometry::Attribute attribute;
      };

      constexpr
      ParseVertex(size_t sdim)
        : m_sdim(sdim)
      {}

      template <class Iterator>
      inline
      std::optional<Data> operator()(Iterator begin, Iterator end) const
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::double_;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::repeat;
        size_t i = 0;
        Data res{ Math::Vector(m_sdim), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE };
        const auto get_x = [&](auto& ctx) { assert(i < m_sdim); res.vertex(i++) = _attr(ctx); };
        const auto get_attribute = [&](auto& ctx) { res.attribute = _attr(ctx); };
        const auto p = double_[get_x] >> repeat(m_sdim - 1)[double_[get_x]] >> uint_[get_attribute];
        const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
        if (begin != end)
          return {};
        else if (r)
          return { res, };
        else
          return {};
      }

    private:
      const size_t m_sdim;
  };

  class UnexpectedKeywordException : public Alert::Exception
  {
    public:
      UnexpectedKeywordException(const std::string& actual, const std::string& expected)
      {
        *this << "Unexpected keyword: " << std::quoted(actual) << ". "
              << "Expected: " << std::quoted(expected) << ".";
      }
  };

  std::ostream& operator<<(std::ostream& os, Keyword kw);

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
      const auto p = *blank;
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

  struct ParseInteger
  {
    template <class Iterator>
    inline
    std::optional<int> operator()(Iterator begin, Iterator end) const
    {
      using boost::spirit::x3::space;
      using boost::spirit::x3::blank;
      using boost::spirit::x3::int_;
      using boost::spirit::x3::_attr;

      int v;
      const auto get_integer = [&](auto& ctx) { v = _attr(ctx); };
      const auto p = int_[get_integer];
      const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);
      if (begin != end)
        return {};
      else if (r)
        return v;
      else
        return {};
    }
  };

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

  struct ParseMeshVersionFormatted
  {
    template <class Iterator>
    inline
    std::optional<unsigned int> operator()(Iterator begin, Iterator end) const
    {
      static constexpr const char* expected = toCharString(Keyword::MeshVersionFormatted);
      using boost::spirit::x3::space;
      using boost::spirit::x3::blank;
      using boost::spirit::x3::uint_;
      using boost::spirit::x3::_attr;
      using boost::spirit::x3::alpha;
      std::string kw;
      const auto get_keyword = [&](auto& ctx) { kw = _attr(ctx); };
      const auto pkw = (+alpha)[get_keyword];
      const bool rkw = boost::spirit::x3::phrase_parse(begin, end, pkw, space);
      if (rkw)
      {
        assert(kw.size() > 0);
        if (kw == expected)
        {
          unsigned int version;
          const auto get_version = [&](auto& ctx) { version = _attr(ctx); };
          const auto pversion = uint_[get_version];
          const bool rversion = boost::spirit::x3::phrase_parse(begin, end, pversion, space);
          if (begin != end)
            return {};
          else if (rversion)
            return version;
          else
            return {};
        }
        else
        {
          throw UnexpectedKeywordException(kw, expected);
        }
      }
      else
      {
        Alert::Exception() << "Failed to parse keyword: "
                           << std::quoted(expected) << '.'
                           << Alert::Raise;
      }
      return {};
    }
  };

  struct ParseDimension
  {
    template <class Iterator>
    inline
    std::optional<unsigned int> operator()(Iterator begin, Iterator end) const
    {
      static constexpr const char* expected = toCharString(Keyword::Dimension);
      using boost::spirit::x3::space;
      using boost::spirit::x3::blank;
      using boost::spirit::x3::uint_;
      using boost::spirit::x3::_attr;
      using boost::spirit::x3::alpha;
      std::string kw;
      const auto get_keyword = [&](auto& ctx) { kw = _attr(ctx); };
      const auto pkw = (+alpha)[get_keyword];
      const bool rkw = boost::spirit::x3::phrase_parse(begin, end, pkw, space);
      if (rkw)
      {
        assert(kw.size() > 0);
        if (kw == expected)
        {
          unsigned int dimension;
          const auto get_dimension = [&](auto& ctx) { dimension = _attr(ctx); };
          const auto pdimension = uint_[get_dimension];
          const bool rdimension = boost::spirit::x3::phrase_parse(begin, end, pdimension, space);
          if (begin != end)
            return {};
          else if (rdimension)
            return dimension;
          else
            return {};
        }
        else
        {
          throw UnexpectedKeywordException(kw, expected);
        }
      }
      else
      {
        Alert::Exception() << "Failed to parse keyword: "
                           << std::quoted(expected) << '.'
                           << Alert::Raise;
      }
      return {};
    }
  };
}

namespace Rodin::IO
{
  /**
   * @ingroup MeshLoaderSpecializations
   * @brief Specialization for loading Sequential meshes in the MEDIT file format.
   *
   * The MEDIT file format specification can be found by visiting
   * <a href="https://www.ljll.math.upmc.fr/frey/logiciels/Docmedit.dir/index.html">this
   * link</a>.
   */
  template <>
  class MeshLoader<IO::FileFormat::MEDIT, Context::Sequential>
    : public MeshLoaderBase<Context::Sequential>
  {
    public:
      using Object = Rodin::Geometry::Mesh<Context::Sequential>;

      MeshLoader(Object& mesh)
        : MeshLoaderBase<Context::Sequential>(mesh),
          m_currentLineNumber(0)
      {}

      void load(std::istream& is) override;

      std::istream& getline(std::istream& is, std::string& line);
      std::string skipEmptyLines(std::istream& is);
      void readVersion(std::istream& is);
      void readDimension(std::istream& is);
      void readEntities(std::istream& is);

      std::unordered_map<MEDIT::Keyword, size_t>& getCountMap()
      {
        return m_count;
      }

      const std::unordered_map<MEDIT::Keyword, size_t>& getCountMap() const
      {
        return m_count;
      }

      std::unordered_map<MEDIT::Keyword, std::istream::pos_type>& getPositionMap()
      {
        return m_pos;
      }

      const std::unordered_map<MEDIT::Keyword, std::istream::pos_type>& getPositionMap() const
      {
        return m_pos;
      }

    private:
      Rodin::Geometry::Mesh<Rodin::Context::Sequential>::Builder m_build;

      size_t m_version;
      size_t m_spaceDimension;
      size_t m_currentLineNumber;

      std::unordered_map<MEDIT::Keyword, std::istream::pos_type> m_pos;
      std::unordered_map<MEDIT::Keyword, size_t> m_count;
  };

  template <>
  class MeshPrinter<FileFormat::MEDIT, Context::Sequential>
    : public MeshPrinterBase<Context::Sequential>
  {
    public:
      MeshPrinter(const Rodin::Geometry::Mesh<Context::Sequential>& mesh)
        : MeshPrinterBase(mesh)
      {}

      void print(std::ostream& os) override
      {
        print(os, true);
      }

      void print(std::ostream& os, bool printEnd);

      void printVersion(std::ostream& os);
      void printDimension(std::ostream& os);
      void printEntities(std::ostream& os);
      void printEnd(std::ostream& os);
  };

  template <class Range>
  class GridFunctionLoader<FileFormat::MEDIT,
        Variational::P1<Range, Context::Sequential, Geometry::Mesh<Context::Sequential>>>
    : public GridFunctionLoaderBase<
        Variational::P1<Range, Context::Sequential, Geometry::Mesh<Context::Sequential>>>
  {
    public:
      /// Type of finite element space
      using FES = Variational::P1<Range, Context::Sequential, Geometry::Mesh<Context::Sequential>>;

      GridFunctionLoader(Variational::GridFunction<FES>& gf)
        : GridFunctionLoaderBase<FES>(gf),
          m_currentLineNumber(0)
      {}

      void load(std::istream& is) override;

      std::istream& getline(std::istream& is, std::string& line);
      std::string skipEmptyLines(std::istream& is);

      void readVersion(std::istream& is);
      void readDimension(std::istream& is);
      void readData(std::istream& is);

    private:
      size_t m_version;
      size_t m_spaceDimension;
      size_t m_currentLineNumber;
  };

  template <class FES>
  class GridFunctionPrinter<FileFormat::MEDIT, FES>
    : public GridFunctionPrinterBase<FES>
  {
    public:
      GridFunctionPrinter(const Variational::GridFunction<FES>& gf)
        : GridFunctionPrinterBase<FES>(gf)
      {}

      void print(std::ostream& os) override
      {
        printVersion(os);
        printDimension(os);
        printData(os);
        printEnd(os);
      }

      void printVersion(std::ostream& os)
      {
        os << MEDIT::Keyword::MeshVersionFormatted << "\n2" << "\n\n";
      }

      void printDimension(std::ostream& os)
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        os << MEDIT::Keyword::Dimension << '\n' << mesh.getSpaceDimension() << "\n\n";
      }

      void printData(std::ostream& os)
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t vdim = fes.getVectorDimension();
        os << MEDIT::Keyword::SolAtVertices << '\n'
           << mesh.getVertexCount() << '\n'
           << 1 // Only one solution
           << " " << ((vdim > 1) ? MEDIT::SolutionType::Vector : MEDIT::SolutionType::Scalar)
           << '\n';

        if constexpr (Utility::IsSpecialization<FES, Variational::P1>::Value)
        {
          os << gf.getData().reshaped();
        }
        else
        {
          for (auto it = mesh.getVertex(); !it.end(); ++it)
          {
            const Geometry::Point p(*it, it->getTransformation(),
                Geometry::Polytope::getVertices(Geometry::Polytope::Type::Point).col(0),
                it->getCoordinates());
            os << gf(p) << '\n';
          }
        }
        os << '\n';
      }

      void printEnd(std::ostream& os)
      {
        os << '\n' << IO::MEDIT::Keyword::End;
      }
  };
}

#include "MEDIT.hpp"
#endif
