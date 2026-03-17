/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_MEDIT_H
#define RODIN_IO_MEDIT_H

#include <iomanip>
#include <unordered_map>
#include <boost/bimap.hpp>
#include <boost/spirit/home/x3.hpp>

#include "Rodin/Geometry/AttributeIndex.h"
#include "Rodin/Types.h"
#include "Rodin/Context/Local.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Geometry/Types.h"

#include "ForwardDecls.h"
#include "MeshLoader.h"
#include "MeshPrinter.h"
#include "GridFunctionLoader.h"
#include "GridFunctionPrinter.h"
#include "Rodin/Variational/P1/ForwardDecls.h"
#include "Rodin/Variational/H1/ForwardDecls.h"

#define RODIN_IO_MEDIT_DEFAULT_POLYTOPE_ATTRIBUTE 0

namespace Rodin::IO::MEDIT
{
  /**
   * @brief Keywords used in MEDIT mesh and solution file formats.
   *
   * These keywords identify different sections in MEDIT files (.mesh and .sol).
   * The MEDIT format is a text-based format used by the MMG remeshing software.
   *
   * @see <a href="https://www.ljll.math.upmc.fr/frey/logiciels/Docmedit.dir/index.html">MEDIT Format Specification</a>
   */
  enum class Keyword
  {
    MeshVersionFormatted,  ///< Format version declaration
    Dimension,             ///< Spatial dimension
    Vertices,              ///< Vertex coordinates section
    Triangles,             ///< Triangle elements section
    Quadrilaterals,        ///< Quadrilateral elements section
    Tetrahedra,            ///< Tetrahedral elements section
    Hexahedra,             ///< Hexahedral elements section
    Wedges,                ///< Wedge (prism) elements section
    Corners,               ///< Corner vertices section
    Ridges,                ///< Ridge edges section
    Edges,                 ///< Edge elements section
    SolAtVertices,         ///< Solution at vertices
    SolAtEdges,            ///< Solution at edges
    SolAtTriangles,        ///< Solution at triangles
    SolAtQuadrilaterals,   ///< Solution at quadrilaterals
    SolAtTetrahedra,       ///< Solution at tetrahedra
    SolAtPentahedra,       ///< Solution at pentahedra
    SolAtHexahedra,        ///< Solution at hexahedra
    RequiredVertices,      ///< Required vertices section
    RequiredEdges,         ///< Required edges section
    Normals,               ///< Normal vectors section
    NormalAtVertices,      ///< Normals at vertices
    Tangents,              ///< Tangent vectors section
    TangentAtVertices,     ///< Tangents at vertices
    End                    ///< End of file marker
  };

  /**
   * @brief Converts a MEDIT keyword enum to its string representation.
   * @param[in] kw Keyword to convert
   * @returns C-style string representation of the keyword
   */
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
      case Keyword::Hexahedra:
        return "Hexahedra";
      case Keyword::Wedges:
        return "Wedges";
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
      case Keyword::RequiredEdges:
        return "RequiredEdges";
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
  Optional<Keyword> toKeyword(const char* str)
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
    else if (str == Keyword::Hexahedra)
      res = Keyword::Hexahedra;
    else if (str == Keyword::Wedges)
      res = Keyword::Wedges;
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
    else if (str == Keyword::RequiredEdges)
      res = Keyword::RequiredEdges;
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

  /**
   * @brief Solution data types in MEDIT solution files.
   *
   * Identifies the type of solution data stored in .sol files.
   */
  enum SolutionType
  {
    Real = 1,    ///< Scalar (real-valued) solution
    Vector = 2,  ///< Vector-valued solution
    Tensor = 3   ///< Tensor-valued solution
  };

  /**
   * @brief Parser for mesh entities (elements) in MEDIT format.
   * @internal
   *
   * Parses element connectivity and attribute information.
   */
  class ParseEntity
  {
    public:
      /**
       * @brief Parsed entity data.
       */
      struct Data
      {
        Array<Index> vertices;       ///< Vertex indices defining the entity
        Geometry::Attribute attribute;  ///< Entity attribute (material ID)
      };

      /**
       * @brief Constructs an entity parser for @p n vertices.
       * @param[in] n Number of vertices in the entity
       */
      constexpr
      ParseEntity(size_t n)
        : m_n(n)
      {}

      /**
       * @brief Parses entity data from an iterator range.
       * @tparam Iterator Iterator type
       * @param[in] begin Start of input range
       * @param[in] end End of input range
       * @returns Optional entity data if parsing succeeds, empty otherwise
       */
      template <class Iterator>
      Optional<Data> operator()(Iterator begin, Iterator end) const
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::repeat;
        size_t i = 0;
        Data res{ Array<Index>(m_n), ~Geometry::Attribute(0) };
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
        Math::SpatialPoint vertex;
        Geometry::Attribute attribute;
      };

      constexpr
      ParseVertex(size_t sdim)
        : m_sdim(sdim)
      {}

      template <class Iterator>
      Optional<Data> operator()(Iterator begin, Iterator end) const
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::double_;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::repeat;
        size_t i = 0;
        Data res{ Math::SpatialPoint(m_sdim), ~Geometry::Attribute(0) };
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
    Optional<std::string> operator()(Iterator begin, Iterator end) const
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
    Optional<int> operator()(Iterator begin, Iterator end) const
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
    Optional<unsigned int> operator()(Iterator begin, Iterator end) const
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
    Optional<unsigned int> operator()(Iterator begin, Iterator end) const
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
    Optional<unsigned int> operator()(Iterator begin, Iterator end) const
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
  class MeshLoader<IO::FileFormat::MEDIT, Context::Local>
    : public MeshLoaderBase<Context::Local>
  {
    public:
      using ObjectType = Rodin::Geometry::Mesh<Context::Local>;

      using Parent = MeshLoaderBase<Context::Local>;

      MeshLoader(ObjectType& mesh)
        : Parent(mesh),
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
      Rodin::Geometry::Mesh<Rodin::Context::Local>::Builder m_build;

      size_t m_version;
      size_t m_spaceDimension;
      size_t m_currentLineNumber;

      std::unordered_map<MEDIT::Keyword, std::istream::pos_type> m_pos;
      std::unordered_map<MEDIT::Keyword, size_t> m_count;
  };

  template <>
  class MeshPrinter<FileFormat::MEDIT, Context::Local>
    : public MeshPrinterBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = MeshPrinterBase<ContextType>;

      MeshPrinter(const ObjectType& mesh)
        : MeshPrinterBase(mesh)
      {}

      void print(std::ostream& os) override
      {
        printMesh(os, true);
      }

      void printMesh(std::ostream& os, bool printEnd);
      void printVersion(std::ostream& os);
      void printDimension(std::ostream& os);
      void printEntities(std::ostream& os);
      void printEnd(std::ostream& os);
  };

  template <class Range>
  class GridFunctionLoader<
    FileFormat::MEDIT,
    Variational::P1<Range, Geometry::Mesh<Context::Local>>,
    Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
    : public GridFunctionLoaderBase<
        Variational::P1<Range, Geometry::Mesh<Context::Local>>,
        Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
  {
    public:
      using FESType = Variational::P1<Range, Geometry::Mesh<Context::Local>>;

      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;

      using DataType = Math::Vector<ScalarType>;

      using ObjectType = Variational::GridFunction<FESType, DataType>;

      using Parent = GridFunctionLoaderBase<FESType, DataType>;

      GridFunctionLoader(ObjectType& gf)
        : Parent(gf),
          m_currentLineNumber(0)
      {}

      void load(std::istream& is) override
      {
        readVersion(is);
        readDimension(is);
        readData(is);
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
          if (!MEDIT::ParseEmptyLine()(line.begin(), line.end()))
            break;
        }
        return line;
      }

      void readVersion(std::istream& is)
      {
        auto line = skipEmptyLines(is);
        Optional<unsigned int> version =
          MEDIT::ParseMeshVersionFormatted()(line.begin(), line.end());
        if (version) // Version was on the same line
        {
          m_version = *version;
        }
        else // Version is not on the same line
        {
          auto line = skipEmptyLines(is);
          version = MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
          if (version)
            m_version = *version;
          else
            Alert::Exception() << "Failed to parse version number of mesh." << Alert::Raise;
        }
      }

      void readDimension(std::istream& is)
      {
        auto line = skipEmptyLines(is);
        Optional<unsigned int> dimension = MEDIT::ParseDimension()(line.begin(), line.end());
        if (dimension) // Version was on the same line
          m_spaceDimension = *dimension;
        else // Version is not on the same line
        {
          auto line = skipEmptyLines(is);
          dimension = MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
          if (dimension)
            m_spaceDimension = *dimension;
          else
            Alert::Exception() << "Failed to parse dimension of mesh." << Alert::Raise;
        }
      }

      void readData(std::istream& is)
      {
        auto& gf = this->getObject();

        auto line = skipEmptyLines(is);
        Optional<std::string> kw =
          MEDIT::ParseKeyword()(line.begin(), line.end());
        if (!kw || *kw != MEDIT::Keyword::SolAtVertices)
        {
          Alert::Exception() << "Expected keyword " << MEDIT::Keyword::SolAtVertices
                             << " on line " << m_currentLineNumber 
                             << Alert::Raise;
        }

        line = skipEmptyLines(is);
        Optional<unsigned int> size = MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
        if (!size)
        {
          Alert::Exception() << "Failed to parse solution size at line "
                             << m_currentLineNumber
                             << Alert::Raise;
        }

        line = skipEmptyLines(is);
        size_t solCount, vdim;
        using boost::spirit::x3::space;
        using boost::spirit::x3::blank;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::repeat;
        const auto get_sol_count = [&](auto& ctx) { solCount = _attr(ctx); };
        const auto get_vdim = [&](auto& ctx) { vdim = _attr(ctx); };
        const auto p = uint_[get_sol_count] >> uint_[get_vdim];
        auto it = line.begin();
        const bool r = boost::spirit::x3::phrase_parse(it, line.end(), p, space);

        (void) solCount;
        assert(solCount == 1);
        if (it != line.end() || !r)
        {
          Alert::Exception() << "Failed to parse solution count and vector dimension at line "
                             << m_currentLineNumber
                             << Alert::Raise;
        }

        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t count = mesh.getVertexCount();
        for (size_t i = 0; i < count; ++i)
          for (size_t d = 0; d < vdim; ++d)
            is >> gf[d * count + i];
      }

    private:
      size_t m_version;
      size_t m_spaceDimension;
      size_t m_currentLineNumber;
  };

  template <size_t K, class Range>
  class GridFunctionLoader<
    FileFormat::MEDIT,
    Variational::H1<K, Range, Geometry::Mesh<Context::Local>>,
    Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
    : public GridFunctionLoaderBase<
        Variational::H1<K, Range, Geometry::Mesh<Context::Local>>,
        Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
  {
    public:
      using FESType    = Variational::H1<K, Range, Geometry::Mesh<Context::Local>>;
      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
      using DataType   = Math::Vector<ScalarType>;
      using ObjectType = Variational::GridFunction<FESType, DataType>;
      using Parent     = GridFunctionLoaderBase<FESType, DataType>;

      GridFunctionLoader(ObjectType& gf)
        : Parent(gf),
          m_version(0),
          m_spaceDimension(0),
          m_currentLineNumber(0)
      {}

      void load(std::istream& is) override
      {
        readVersion(is);
        readDimension(is);
        readData(is);
      }

    private:
      // -------------------------------------------------------------
      // Line helpers (same style as P1 loader)
      // -------------------------------------------------------------
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
          if (!MEDIT::ParseEmptyLine()(line.begin(), line.end()))
            break;
        }
        return line;
      }

      void readVersion(std::istream& is)
      {
        auto line = skipEmptyLines(is);
        Optional<unsigned int> version =
          MEDIT::ParseMeshVersionFormatted()(line.begin(), line.end());
        if (version)
        {
          m_version = *version;
        }
        else
        {
          auto line2 = skipEmptyLines(is);
          version = MEDIT::ParseUnsignedInteger()(line2.begin(), line2.end());
          if (version)
            m_version = *version;
          else
            Alert::Exception()
              << "Failed to parse version number of mesh."
              << Alert::Raise;
        }
      }

      void readDimension(std::istream& is)
      {
        auto line = skipEmptyLines(is);
        Optional<unsigned int> dimension =
          MEDIT::ParseDimension()(line.begin(), line.end());
        if (dimension)
        {
          m_spaceDimension = *dimension;
        }
        else
        {
          auto line2 = skipEmptyLines(is);
          dimension = MEDIT::ParseUnsignedInteger()(line2.begin(), line2.end());
          if (dimension)
            m_spaceDimension = *dimension;
          else
            Alert::Exception()
              << "Failed to parse dimension of mesh."
              << Alert::Raise;
        }
      }

      void readData(std::istream& is)
      {
        auto& gf   = this->getObject();
        auto& data = gf.getData();
        const auto& fes  = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        const size_t vdim       = fes.getVectorDimension();
        const size_t D          = mesh.getDimension();
        const size_t scalarSize = fes.getSize() / vdim;

        // 1. SolAtVertices keyword
        auto line = skipEmptyLines(is);
        Optional<std::string> kw =
          MEDIT::ParseKeyword()(line.begin(), line.end());
        if (!kw || *kw != MEDIT::Keyword::SolAtVertices)
        {
          Alert::Exception()
            << "Expected keyword " << MEDIT::Keyword::SolAtVertices
            << " on line " << m_currentLineNumber
            << Alert::Raise;
        }

        // 2. Number of vertices
        line = skipEmptyLines(is);
        Optional<unsigned int> size =
          MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
        if (!size)
        {
          Alert::Exception()
            << "Failed to parse solution size at line "
            << m_currentLineNumber
            << Alert::Raise;
        }

        const size_t nVertices = mesh.getVertexCount();
        if (static_cast<size_t>(*size) != nVertices)
        {
          Alert::Exception()
            << "MEDIT SolAtVertices size (" << *size
            << ") does not match mesh vertex count (" << nVertices << ")."
            << Alert::Raise;
        }

        // 3. solCount, vdim line
        line = skipEmptyLines(is);
        size_t solCount = 0;
        size_t fileVdim = 0;

        using boost::spirit::x3::space;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;

        const auto get_sol_count = [&](auto& ctx) { solCount  = _attr(ctx); };
        const auto get_vdim      = [&](auto& ctx) { fileVdim = _attr(ctx); };

        const auto p = uint_[get_sol_count] >> uint_[get_vdim];
        auto it = line.begin();
        const bool r = boost::spirit::x3::phrase_parse(it, line.end(), p, space);

        if (!r || it != line.end())
        {
          Alert::Exception()
            << "Failed to parse solution count and vector dimension at line "
            << m_currentLineNumber
            << Alert::Raise;
        }

        if (solCount != 1)
        {
          Alert::Exception()
            << "Only a single solution is supported in SolAtVertices, got "
            << solCount
            << Alert::Raise;
        }

        if (fileVdim != vdim)
        {
          Alert::Exception()
            << "Vector dimension mismatch: file vdim = "
            << fileVdim << ", H1 vdim = " << vdim
            << Alert::Raise;
        }

        if (data.size() == 0)
          return;
        assert(static_cast<size_t>(data.size()) == vdim * scalarSize);

        // 4. Read vertex-based P1 data
        std::vector<ScalarType> vertexValues(vdim * nVertices);
        for (size_t v = 0; v < nVertices; ++v)
        {
          for (size_t c = 0; c < vdim; ++c)
          {
            ScalarType val;
            is >> val;
            vertexValues[c * nVertices + v] = val;
          }
        }

        // -------------------------------------------------------------
        // Special case: 0D mesh (only points).
        //
        // No edges, no cells. Only vertex DOFs exist in H1, so we
        // simply copy SolAtVertices into those DOFs.
        // -------------------------------------------------------------
        if (D == 0)
        {
          for (size_t v = 0; v < nVertices; ++v)
          {
            const auto& vdofs = fes.getDOFs(0, static_cast<Index>(v));
            assert(vdofs.size() == 1);

            Index d = vdofs(0);
            for (size_t c = 0; c < vdim; ++c)
            {
              data.coeffRef(d + static_cast<Index>(c * scalarSize)) =
                vertexValues[c * nVertices + v];
            }
          }
          return;
        }

        // -------------------------------------------------------------
        // Generic D >= 1: element-based injection
        // -------------------------------------------------------------

        const auto& connD0  = mesh.getConnectivity().getIncidence(D, 0);
        const size_t nCells = mesh.getConnectivity().getCount(D);

        // P1/Q1 evaluators on reference elements

        auto evalP1_on_segment = [&](Real x,
                                     const ScalarType* u) -> ScalarType
        {
          // Reference segment [0,1]; vertex 0 at x=0, vertex 1 at x=1
          return (ScalarType(1) - x) * u[0] + x * u[1];
        };

        auto evalP1_on_triangle = [&](Real x, Real y,
                                      const ScalarType* u) -> ScalarType
        {
          // Ref triangle (0,0)-(1,0)-(0,1)
          const Real l0 = Real(1) - x - y;
          const Real l1 = x;
          const Real l2 = y;
          return l0 * u[0] + l1 * u[1] + l2 * u[2];
        };

        auto evalQ1_on_quad = [&](Real x, Real y,
                                  const ScalarType* u) -> ScalarType
        {
          // Ref quad [0,1]^2 with vertices:
          // 0:(0,0), 1:(1,0), 2:(1,1), 3:(0,1)
          const Real one_x = Real(1) - x;
          const Real one_y = Real(1) - y;

          const Real psi0 = one_x * one_y;
          const Real psi1 = x      * one_y;
          const Real psi2 = x      * y;
          const Real psi3 = one_x  * y;

          return psi0 * u[0] + psi1 * u[1] + psi2 * u[2] + psi3 * u[3];
        };

        auto evalP1_on_tet = [&](Real x, Real y, Real z,
                                 const ScalarType* u) -> ScalarType
        {
          // Ref tet (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1)
          const Real l0 = Real(1) - x - y - z;
          const Real l1 = x;
          const Real l2 = y;
          const Real l3 = z;
          return l0 * u[0] + l1 * u[1] + l2 * u[2] + l3 * u[3];
        };

        auto evalP1_on_wedge = [&](Real x, Real y, Real z,
                                   const ScalarType* u) -> ScalarType
        {
          // Ref wedge = triangle(x,y) × segment(z)
          //
          // Triangle barycentric:
          //   λ0 = 1 - x - y, λ1 = x, λ2 = y
          // Segment z ∈ [0,1].
          //
          // Vertices:
          //   v0: (λ0, z=0), v1: (λ1, 0), v2: (λ2, 0)
          //   v3: (λ0, 1),   v4: (λ1, 1), v5: (λ2, 1)
          const Real l0 = Real(1) - x - y;
          const Real l1 = x;
          const Real l2 = y;

          const Real one_z = Real(1) - z;

          const Real psi0 = l0 * one_z;
          const Real psi1 = l1 * one_z;
          const Real psi2 = l2 * one_z;
          const Real psi3 = l0 * z;
          const Real psi4 = l1 * z;
          const Real psi5 = l2 * z;

          return psi0 * u[0] + psi1 * u[1] + psi2 * u[2]
               + psi3 * u[3] + psi4 * u[4] + psi5 * u[5];
        };

        // Element loop
        for (Index c = 0; c < static_cast<Index>(nCells); ++c)
        {
          const auto geom   = mesh.getGeometry(D, c);
          const auto& verts = connD0[c];
          const auto& cdofs = fes.getDOFs(D, c);

          const size_t nLocal = static_cast<size_t>(cdofs.size());

          // H1 local nodes in reference element
          const auto& nodes =
            Variational::H1Element<K, typename FormLanguage::Traits<Range>::ScalarType>::getNodes(geom);
          assert(nodes.size() == nLocal);

          // Number of vertices per cell (topology)
          size_t nCellVertices = 0;
          switch (geom)
          {
            case Geometry::Polytope::Type::Point:
              nCellVertices = 1;
              break;
            case Geometry::Polytope::Type::Segment:
              nCellVertices = 2;
              break;
            case Geometry::Polytope::Type::Triangle:
              nCellVertices = 3;
              break;
            case Geometry::Polytope::Type::Quadrilateral:
              nCellVertices = 4;
              break;
            case Geometry::Polytope::Type::Tetrahedron:
              nCellVertices = 4;
              break;
            case Geometry::Polytope::Type::Wedge:
              nCellVertices = 6;
              break;
            default:
              Alert::Exception()
                << "Unsupported cell geometry in MEDIT H1 loader."
                << Alert::Raise;
          }

          assert(verts.size() == nCellVertices);

          // Per-component injection
          for (size_t comp = 0; comp < vdim; ++comp)
          {
            // Gather vertex values u_vert[0..nCellVertices-1]
            ScalarType u_vert[8];
            for (size_t lv = 0; lv < nCellVertices; ++lv)
            {
              const Index vIdx = verts[lv];
              assert(static_cast<size_t>(vIdx) < nVertices);
              u_vert[lv] =
                vertexValues[comp * nVertices + static_cast<size_t>(vIdx)];
            }

            // Interpolate at each H1 nodal point of this cell
            for (size_t j = 0; j < nLocal; ++j)
            {
              const auto& pt = nodes[j];
              ScalarType u_val{};

              switch (geom)
              {
                case Geometry::Polytope::Type::Point:
                {
                  // Single vertex, nothing to interpolate
                  u_val = u_vert[0];
                  break;
                }
                case Geometry::Polytope::Type::Segment:
                {
                  const Real x = pt.x();
                  u_val = evalP1_on_segment(x, u_vert);
                  break;
                }
                case Geometry::Polytope::Type::Triangle:
                {
                  const Real x = pt.x();
                  const Real y = pt.y();
                  u_val = evalP1_on_triangle(x, y, u_vert);
                  break;
                }
                case Geometry::Polytope::Type::Quadrilateral:
                {
                  const Real x = pt.x();
                  const Real y = pt.y();
                  u_val = evalQ1_on_quad(x, y, u_vert);
                  break;
                }
                case Geometry::Polytope::Type::Tetrahedron:
                {
                  const Real x = pt.x();
                  const Real y = pt.y();
                  const Real z = pt.z();
                  u_val = evalP1_on_tet(x, y, z, u_vert);
                  break;
                }
                case Geometry::Polytope::Type::Wedge:
                {
                  const Real x = pt.x();
                  const Real y = pt.y();
                  const Real z = pt.z();
                  u_val = evalP1_on_wedge(x, y, z, u_vert);
                  break;
                }
                default:
                  assert(false);
              }

              const Index gdof = cdofs(static_cast<Index>(j));
              assert(static_cast<size_t>(gdof) < scalarSize);
              data.coeffRef(gdof + static_cast<Index>(comp * scalarSize)) = u_val;
            }
          }
        }
      }

      size_t m_version;
      size_t m_spaceDimension;
      size_t m_currentLineNumber;
  };

  template <class FES, class Data>
  class GridFunctionPrinterBase<FileFormat::MEDIT, FES, Data>
    : public Printer<Variational::GridFunction<FES, Data>>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = FileFormat::MEDIT;

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
        printVersion(os);
        printDimension(os);

        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t vdim = fes.getVectorDimension();

        os << MEDIT::Keyword::SolAtVertices << '\n'
           << mesh.getVertexCount() << '\n'
           << 1 // Only one solution
           << " " << ((vdim > 1) ? MEDIT::SolutionType::Vector : MEDIT::SolutionType::Real)
           << '\n';

        this->printData(os);

        os << '\n';

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

      void printEnd(std::ostream& os)
      {
        os << '\n' << IO::MEDIT::Keyword::End;
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
    FileFormat::MEDIT, FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
    : public GridFunctionPrinterBase<
        FileFormat::MEDIT, FES, Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
  {
    public:
      using FESType = FES;

      static constexpr FileFormat Format = FileFormat::MEDIT;

      using RangeType = typename FormLanguage::Traits<FES>::RangeType;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using DataType = Math::Vector<ScalarType>;

      using Parent = GridFunctionPrinterBase<Format, FES, DataType>;

      using Parent::Parent;

      void printData(std::ostream& os)
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const Geometry::Polytope::Traits ts(Geometry::Polytope::Type::Point);
        for (auto it = mesh.getVertex(); !it.end(); ++it)
        {
          const Geometry::Point p(
              *it,
              ts.getVertex(0),
              it->getCoordinates());
          os << gf(p) << '\n';
        }
        os << '\n';
      }
  };
}

#endif
