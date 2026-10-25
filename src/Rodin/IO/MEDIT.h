/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_MEDIT_H
#define RODIN_IO_MEDIT_H

#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <boost/bimap.hpp>
#include <boost/spirit/home/x3.hpp>

#include "Rodin/Alert/MemberFunctionException.h"
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
    Pyramids,              ///< Pyramid elements section
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
    RequiredTriangles,     ///< Required triangles section
    RequiredTetrahedra,    ///< Required tetrahedra section
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
  inline constexpr const char* toCharString(Keyword kw)
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
      case Keyword::Pyramids:
        return "Pyramids";
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
      case Keyword::RequiredTriangles:
        return "RequiredTriangles";
      case Keyword::RequiredTetrahedra:
        return "RequiredTetrahedra";
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

  /**
   * @brief Compares the operands for equality.
   * @param str Text to compare or convert to a format keyword.
   * @param kw Format keyword to compare or write.
   * @returns True when the operands compare equal.
   */

  inline bool operator==(const std::string& str, Keyword kw)
  {
    return str == toCharString(kw);
  }

  /**
   * @brief Compares the operands for inequality.
   * @param str Text to compare or convert to a format keyword.
   * @param kw Format keyword to compare or write.
   * @returns True when the operands compare unequal.
   */

  inline bool operator!=(const std::string& str, Keyword kw)
  {
    return str != toCharString(kw);
  }

  /**
   * @brief Compares the operands for equality.
   * @param kw Format keyword to compare or write.
   * @param str Text to compare or convert to a format keyword.
   * @returns True when the operands compare equal.
   */

  inline bool operator==(Keyword kw, const std::string& str)
  {
    return str == toCharString(kw);
  }

  /**
   * @brief Compares the operands for inequality.
   * @param kw Format keyword to compare or write.
   * @param str Text to compare or convert to a format keyword.
   * @returns True when the operands compare unequal.
   */

  inline bool operator!=(Keyword kw, const std::string& str)
  {
    return str != toCharString(kw);
  }

  /**
   * @brief Compares the operands for equality.
   * @param kw Format keyword to compare or write.
   * @param str Text to compare or convert to a format keyword.
   * @returns True when the operands compare equal.
   */

  inline bool operator==(Keyword kw, const char* str)
  {
    return strcmp(toCharString(kw), str) == 0;
  }

  /**
   * @brief Compares the operands for inequality.
   * @param kw Format keyword to compare or write.
   * @param str Text to compare or convert to a format keyword.
   * @returns True when the operands compare unequal.
   */

  inline bool operator!=(Keyword kw, const char* str)
  {
    return strcmp(toCharString(kw), str) != 0;
  }

  /**
   * @brief Compares the operands for equality.
   * @param str Text to compare or convert to a format keyword.
   * @param kw Format keyword to compare or write.
   * @returns True when the operands compare equal.
   */

  inline bool operator==(const char* str, Keyword kw)
  {
    return strcmp(toCharString(kw), str) == 0;
  }

  /**
   * @brief Compares the operands for inequality.
   * @param str Text to compare or convert to a format keyword.
   * @param kw Format keyword to compare or write.
   * @returns True when the operands compare unequal.
   */

  inline bool operator!=(const char* str, Keyword kw)
  {
    return strcmp(toCharString(kw), str) != 0;
  }

  /**
   * @brief Writes the format keyword to a stream.
   * @param os Output stream receiving the formatted data.
   * @param kw Format keyword to compare or write.
   * @returns Reference to the output stream after insertion.
   */

  inline std::ostream& operator<<(std::ostream& os, Keyword kw)
  {
    os << toCharString(kw);
    return os;
  }

  /**
   * @brief Converts text to a format keyword.
   * @param str Text to compare or convert to a format keyword.
   * @returns Matching keyword, or an empty optional if the text is not recognized.
   */

  inline Optional<Keyword> toKeyword(const char* str)
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
    else if (str == Keyword::Pyramids)
      res = Keyword::Pyramids;
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
    else if (str == Keyword::RequiredTriangles)
      res = Keyword::RequiredTriangles;
    else if (str == Keyword::RequiredTetrahedra)
      res = Keyword::RequiredTetrahedra;
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
      /// @brief Parsed entity data.
      struct Data
      {
        Array<Index> vertices;       ///< Vertex indices defining the entity
        Geometry::Attribute attribute;  ///< Entity attribute (material ID)
      };

      /**
       * @brief Constructs an entity parser for @p n vertices.
       * @param[in] n Number of vertices in the entity
       */
      constexpr ParseEntity(size_t n)
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
        Data res{Array<Index>(m_n), ~Geometry::Attribute(0)};
        const auto getVertex = [&](auto& ctx) {
          assert(i < m_n);
          res.vertices(i++) = _attr(ctx);
        };
        const auto getAttribute = [&](auto& ctx) { res.attribute = _attr(ctx); };
        const auto p =
          uint_[getVertex] >> repeat(m_n - 1)[uint_[getVertex]] >> uint_[getAttribute];
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

  /// @cond
  class ParseVertex
  {
    public:
      struct Data
      {
        Math::SpatialPoint vertex;
        Geometry::Attribute attribute;
      };

      constexpr ParseVertex(size_t sdim)
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
        Data res{Math::SpatialPoint(m_sdim), ~Geometry::Attribute(0)};
        const auto getX = [&](auto& ctx) {
          assert(i < m_sdim);
          res.vertex(i++) = _attr(ctx);
        };
        const auto getAttribute = [&](auto& ctx) { res.attribute = _attr(ctx); };
        const bool r = [&]() {
          if (m_sdim == 0)
          {
            return boost::spirit::x3::phrase_parse(
              begin, end, uint_[getAttribute], space);
          }
          const auto p =
            double_[getX] >> repeat(m_sdim - 1)[double_[getX]] >> uint_[getAttribute];
          return boost::spirit::x3::phrase_parse(begin, end, p, space);
        }();
        if (begin != end)
          return {};
        else if (r)
        {
          return {
            res,
          };
        }
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
      const auto getKeyword = [&](auto& ctx) { kw = _attr(ctx); };
      const auto p = boost::spirit::x3::lexeme[+alpha][getKeyword];
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
      const auto getInteger = [&](auto& ctx) { v = _attr(ctx); };
      const auto p = int_[getInteger];
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
      const auto getUnsignedInteger = [&](auto& ctx) { v = _attr(ctx); };
      const auto p = uint_[getUnsignedInteger];
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
      static constexpr const char* expected =
        toCharString(Keyword::MeshVersionFormatted);
      using boost::spirit::x3::space;
      using boost::spirit::x3::blank;
      using boost::spirit::x3::uint_;
      using boost::spirit::x3::_attr;
      using boost::spirit::x3::alpha;
      std::string kw;
      const auto getKeyword = [&](auto& ctx) { kw = _attr(ctx); };
      const auto pkw = boost::spirit::x3::lexeme[+alpha][getKeyword];
      const bool rkw = boost::spirit::x3::phrase_parse(begin, end, pkw, space);
      if (rkw)
      {
        assert(kw.size() > 0);
        if (kw == expected)
        {
          unsigned int version;
          const auto getVersion = [&](auto& ctx) { version = _attr(ctx); };
          const auto pversion = uint_[getVersion];
          const bool rversion =
            boost::spirit::x3::phrase_parse(begin, end, pversion, space);
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
        Alert::MemberFunctionException(*this, __func__)
          << "Failed to parse keyword: " << std::quoted(expected) << '.'
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
      const auto getKeyword = [&](auto& ctx) { kw = _attr(ctx); };
      const auto pkw = boost::spirit::x3::lexeme[+alpha][getKeyword];
      const bool rkw = boost::spirit::x3::phrase_parse(begin, end, pkw, space);
      if (rkw)
      {
        assert(kw.size() > 0);
        if (kw == expected)
        {
          unsigned int dimension;
          const auto getDimension = [&](auto& ctx) { dimension = _attr(ctx); };
          const auto pdimension = uint_[getDimension];
          const bool rdimension =
            boost::spirit::x3::phrase_parse(begin, end, pdimension, space);
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
        Alert::MemberFunctionException(*this, __func__)
          << "Failed to parse keyword: " << std::quoted(expected) << '.'
          << Alert::Raise;
      }
      return {};
    }
  };

  /// @endcond
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
      /// @brief Mesh type loaded by this specialization.
      using ObjectType = Rodin::Geometry::Mesh<Context::Local>;

      /// @brief Parent class type.
      using Parent = MeshLoaderBase<Context::Local>;

      /**
       * @brief Constructs a MEDIT loader bound to a local mesh.
       * @param[in,out] mesh Mesh object populated by load().
       */
      MeshLoader(ObjectType& mesh)
        : Parent(mesh),
          m_currentLineNumber(0)
      {}

      /**
       * @brief Loads a local mesh from a MEDIT input stream.
       * @param[in,out] is Input stream containing MEDIT data.
       */
      void load(std::istream& is) override;

      /// @brief Reads one logical line while tracking line numbers.
      /// @param is Input stream containing the MEDIT mesh data.
      /// @param line Output buffer receiving the next input line.
      /// @returns Input stream after reading the next logical line.
      std::istream& getline(std::istream& is, std::string& line);

      /// @brief Skips blank lines and returns the next non-empty line.
      /// @param is Input stream containing the MEDIT mesh data.
      /// @returns Next non-empty input line.
      std::string skipEmptyLines(std::istream& is);

      /// @brief Reads and validates the MEDIT version section.
      /// @param is Input stream containing the MEDIT mesh data.
      void readVersion(std::istream& is);

      /// @brief Reads and validates the MEDIT dimension section.
      /// @param is Input stream containing the MEDIT mesh data.
      void readDimension(std::istream& is);

      /// @brief Reads vertex and element entities from the MEDIT stream.
      /// @param is Input stream containing the MEDIT mesh data.
      void readEntities(std::istream& is);

      /// @brief Returns mutable counts collected for each MEDIT keyword.
      /// @returns Map from MEDIT keywords to parsed entity counts.
      std::unordered_map<MEDIT::Keyword, size_t>& getCountMap()
      {
        return m_count;
      }

      /// @brief Returns counts collected for each MEDIT keyword.
      /// @returns Map from MEDIT keywords to parsed entity counts.
      const std::unordered_map<MEDIT::Keyword, size_t>& getCountMap() const
      {
        return m_count;
      }

      /// @brief Returns mutable stream positions for parsed MEDIT sections.
      /// @returns Map from MEDIT keywords to their input-stream positions.
      std::unordered_map<MEDIT::Keyword, std::istream::pos_type>& getPositionMap()
      {
        return m_pos;
      }

      /// @brief Returns stream positions for parsed MEDIT sections.
      /// @returns Map from MEDIT keywords to their input-stream positions.
      const std::unordered_map<MEDIT::Keyword, std::istream::pos_type>&
      getPositionMap() const
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

  /**
   * @ingroup PrinterSpecializations
   *
   * @brief MEDIT mesh printer for local meshes.
   */
  template <>
  class MeshPrinter<FileFormat::MEDIT, Context::Local>
    : public MeshPrinterBase<Context::Local>
  {
    public:
      /// @brief Execution context type.
      using ContextType = Context::Local;

      /// @brief Mesh type being printed.
      using ObjectType = Geometry::Mesh<ContextType>;

      /// @brief Parent class type.
      using Parent = MeshPrinterBase<ContextType>;

      /**
       * @brief Constructs a MEDIT mesh printer.
       * @param[in] mesh Mesh to print.
       */
      MeshPrinter(const ObjectType& mesh)
        : MeshPrinterBase(mesh)
      {}

      void print(std::ostream& os) override
      {
        printMesh(os, true);
      }

      /**
       * @brief Prints the full MEDIT mesh body.
       * @param os Output stream.
       * @param printEnd Whether to write the mesh terminator.
       */
      void printMesh(std::ostream& os, bool printEnd);

      /**
       * @brief Prints the MEDIT file version section.
       * @param os Output stream.
       */
      void printVersion(std::ostream& os);

      /**
       * @brief Prints the MEDIT dimension section.
       * @param os Output stream.
       */
      void printDimension(std::ostream& os);

      /**
       * @brief Prints all supported mesh entity sections.
       * @param os Output stream.
       */
      void printEntities(std::ostream& os);

      /**
       * @brief Prints the MEDIT end marker.
       * @param os Output stream.
       */
      void printEnd(std::ostream& os);
  };

  /// @brief MEDIT grid-function loader for local P1 finite element spaces.
  template <class Range>
    requires(!FormLanguage::IsMatrixRange<Range>::Value)
  class GridFunctionLoader<FileFormat::MEDIT,
    Variational::P1<Range, Geometry::Mesh<Context::Local>>,
    Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
    : public GridFunctionLoaderBase<
        Variational::P1<Range, Geometry::Mesh<Context::Local>>,
        Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = Variational::P1<Range, Geometry::Mesh<Context::Local>>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;

      /// @brief Coefficient data storage type.
      using DataType = Math::Vector<ScalarType>;

      /// @brief Grid function type being loaded.
      using ObjectType = Variational::GridFunction<FESType, DataType>;

      /// @brief Parent class type.
      using Parent = GridFunctionLoaderBase<FESType, DataType>;

      /**
       * @brief Constructs a MEDIT P1 grid-function loader.
       * @param[in,out] gf Grid function to populate.
       */
      GridFunctionLoader(ObjectType& gf)
        : Parent(gf),
          m_version(0),
          m_spaceDimension(0),
          m_currentLineNumber(0)
      {}

      /**
       * @brief Loads a complete vertex payload before replacing the field coefficients.
       * @param is Input solution stream.
       */
      void load(std::istream& is) override
      {
        m_currentLineNumber = 0;
        readVersion(is);
        readDimension(is);
        readData(is);
      }

      /**
       * @brief Reads one input line and advances the line counter.
       * @param is Input stream.
       * @param line Input line to parse.
       * @returns The line.
       */
      std::istream& getline(std::istream& is, std::string& line)
      {
        m_currentLineNumber++;
        return std::getline(is, line);
      }

      /**
       * @brief Skips blank lines and returns the first nonblank line.
       * @param is Input stream.
       * @returns First nonblank input line, or an empty string when the stream ends.
       */
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

      /**
       * @brief Reads the MEDIT solution-file version section.
       * @param is Input stream.
       */
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
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to parse version number of mesh." << Alert::Raise;
          }
        }
      }

      /**
       * @brief Reads the MEDIT solution-file dimension section.
       * @param is Input stream.
       */
      void readDimension(std::istream& is)
      {
        auto line = skipEmptyLines(is);
        Optional<unsigned int> dimension =
          MEDIT::ParseDimension()(line.begin(), line.end());
        if (dimension) // Version was on the same line
          m_spaceDimension = *dimension;
        else // Version is not on the same line
        {
          auto line = skipEmptyLines(is);
          dimension = MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
          if (dimension)
            m_spaceDimension = *dimension;
          else
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to parse dimension of mesh." << Alert::Raise;
          }
        }
        if (m_spaceDimension !=
          this->getObject().getFiniteElementSpace().getMesh().getSpaceDimension())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MEDIT dimension does not match the finite element mesh." << Alert::Raise;
        }
      }

      /**
       * @brief Reads the MEDIT SolAtVertices data section.
       * @param is Input stream.
       */
      void readData(std::istream& is)
      {
        auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        auto line = skipEmptyLines(is);
        Optional<std::string> kw = MEDIT::ParseKeyword()(line.begin(), line.end());
        if (!kw || *kw != MEDIT::Keyword::SolAtVertices)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Expected keyword " << MEDIT::Keyword::SolAtVertices << " on line "
            << m_currentLineNumber << Alert::Raise;
        }

        line = skipEmptyLines(is);
        Optional<unsigned int> size =
          MEDIT::ParseUnsignedInteger()(line.begin(), line.end());
        if (!size)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Failed to parse solution size at line " << m_currentLineNumber
            << Alert::Raise;
        }

        if (static_cast<size_t>(*size) != mesh.getVertexCount())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MEDIT vertex count does not match the finite element mesh."
            << Alert::Raise;
        }

        line = skipEmptyLines(is);
        const size_t vdim = fes.getVectorDimension();
        size_t solCount = 0;
        std::istringstream descriptor(line);
        bool valid = static_cast<bool>(descriptor >> solCount) &&
          (solCount == 1 || solCount == vdim) && solCount > 0;
        for (size_t field = 0; valid && field < solCount; ++field)
        {
          size_t type = 0;
          valid = static_cast<bool>(descriptor >> type) &&
            ((solCount == vdim && type == MEDIT::SolutionType::Real) ||
              (solCount == 1 && type == MEDIT::SolutionType::Vector &&
                vdim == mesh.getSpaceDimension()));
        }
        if (!valid || !(descriptor >> std::ws).eof())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MEDIT solution fields do not match their declared component layout."
            << Alert::Raise;
        }
        const size_t count = mesh.getVertexCount();
        DataType values(gf.getData().size());
        for (size_t i = 0; i < count; ++i)
        {
          const auto& dofs = fes.getDOFs(0, i);
          for (size_t d = 0; d < vdim; ++d)
          {
            if (!(is >> values[dofs[d]]))
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Invalid MEDIT coefficient at vertex " << i << ", component " << d
                << "." << Alert::Raise;
            }
          }
        }
        gf.setData(values);
      }

    private:
      size_t m_version;
      size_t m_spaceDimension;
      size_t m_currentLineNumber;
  };

  /// @brief MEDIT grid-function loader for local H1 finite element spaces.
  template <size_t K, class Range>
    requires(!FormLanguage::IsMatrixRange<Range>::Value)
  class GridFunctionLoader<FileFormat::MEDIT,
    Variational::H1<K, Range, Geometry::Mesh<Context::Local>>,
    Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
    : public GridFunctionLoaderBase<
        Variational::H1<K, Range, Geometry::Mesh<Context::Local>>,
        Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = Variational::H1<K, Range, Geometry::Mesh<Context::Local>>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;

      /// @brief Coefficient data storage type.
      using DataType = Math::Vector<ScalarType>;

      /// @brief Grid function type being loaded.
      using ObjectType = Variational::GridFunction<FESType, DataType>;

      /// @brief Parent loader base type.
      using Parent = GridFunctionLoaderBase<FESType, DataType>;

      /**
       * @brief Constructs a MEDIT H1 grid-function loader.
       * @param[in,out] gf Grid function to populate.
       */
      GridFunctionLoader(ObjectType& gf)
        : Parent(gf)
      {}

      /**
       * @brief Interpolates the parsed P1 vertex field using the native H1 interpolation.
       * @param is Input solution stream.
       */
      void load(std::istream& is) override
      {
        auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        using VertexSpace = Variational::P1<Range, Geometry::Mesh<Context::Local>>;
        const auto vertexSpace = [&] {
          if constexpr (FormLanguage::IsVectorRange<Range>::Value)
            return VertexSpace(fes.getMesh(), fes.getVectorDimension());
          else
            return VertexSpace(fes.getMesh());
        }();
        Variational::GridFunction vertexField(vertexSpace);
        GridFunctionLoader<FileFormat::MEDIT, VertexSpace, DataType>(vertexField)
          .load(is);
        gf = vertexField;
      }
  };

  /**
   * @brief Base class for MEDIT grid-function printers.
   *
   * Emits the common MEDIT solution-file header and footer and delegates the
   * actual coefficient stream to @ref printData().
   */
  template <class FES, class Data>
  class GridFunctionPrinterBase<FileFormat::MEDIT, FES, Data>
    : public Printer<Variational::GridFunction<FES, Data>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief File format handled by this printer base.
      static constexpr FileFormat Format = FileFormat::MEDIT;

      /// @brief Range (evaluation value) type.
      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      /// @brief Coefficient data storage type.
      using DataType = Data;

      /// @brief Grid function type being printed.
      using ObjectType = Variational::GridFunction<FESType, DataType>;

      /// @brief Parent class type.
      using Parent = Printer<ObjectType>;

      /**
       * @brief Constructs a MEDIT grid-function printer base.
       * @param[in] gf Grid function to print.
       */
      GridFunctionPrinterBase(const ObjectType& gf)
        : m_gf(gf)
      {}

      void print(std::ostream& os) override
      {
        if constexpr (std::is_same_v<ScalarType, Complex>)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "MEDIT solution files encode real fields; complex fields require "
            << "separate real and imaginary exports." << Alert::Raise;
        }
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t vdim = fes.getVectorDimension();
        printVersion(os);
        printDimension(os);

        os << MEDIT::Keyword::SolAtVertices << '\n' << mesh.getVertexCount() << '\n';
        if (FormLanguage::IsMatrixRange<RangeType>::Value ||
          (vdim > 1 && vdim != mesh.getSpaceDimension()))
        {
          os << vdim;
          for (size_t c = 0; c < vdim; ++c)
            os << " " << MEDIT::SolutionType::Real;
          os << '\n';
        }
        else
        {
          os << 1 << " "
             << ((vdim > 1) ? MEDIT::SolutionType::Vector : MEDIT::SolutionType::Real)
             << '\n';
        }

        this->printData(os);

        os << '\n';

        printEnd(os);
      }

      /**
       * @brief Prints the MEDIT solution-file version section.
       * @param os Output stream.
       */
      void printVersion(std::ostream& os)
      {
        os << MEDIT::Keyword::MeshVersionFormatted << "\n2" << "\n\n";
      }

      /**
       * @brief Prints the MEDIT solution-file dimension section.
       * @param os Output stream.
       */
      void printDimension(std::ostream& os)
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        os << MEDIT::Keyword::Dimension << '\n' << mesh.getSpaceDimension() << "\n\n";
      }

      /**
       * @brief Prints the MEDIT end marker.
       * @param os Output stream.
       */
      void printEnd(std::ostream& os)
      {
        os << '\n' << IO::MEDIT::Keyword::End;
      }

      const ObjectType& getObject() const override
      {
        return m_gf.get();
      }

      /**
       * @brief Prints the concrete grid-function coefficient data.
       * @param os Output stream.
       */
      virtual void printData(std::ostream& os) = 0;

    private:
      std::reference_wrapper<const ObjectType> m_gf;
  };

  /// @brief MEDIT grid-function printer for vector-backed grid functions.
  template <class FES>
  class GridFunctionPrinter<FileFormat::MEDIT, FES,
    Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
    : public GridFunctionPrinterBase<FileFormat::MEDIT, FES,
        Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief File format handled by this printer.
      static constexpr FileFormat Format = FileFormat::MEDIT;

      /// @brief Range (evaluation value) type.
      using RangeType = typename FormLanguage::Traits<FES>::RangeType;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      /// @brief Coefficient data storage type.
      using DataType = Math::Vector<ScalarType>;

      /// @brief Parent class type.
      using Parent = GridFunctionPrinterBase<Format, FES, DataType>;

      using Parent::Parent;

      /**
       * @brief Prints one value per mesh vertex in MEDIT solution order.
       * @param os Output stream.
       */
      void printData(std::ostream& os)
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const Geometry::Polytope::Traits ts(Geometry::Polytope::Type::Point);
        for (auto it = mesh.getVertex(); !it.end(); ++it)
        {
          const Geometry::Point p(*it, ts.getVertex(0), it->getCoordinates());
          if constexpr (FormLanguage::IsMatrixRange<RangeType>::Value)
          {
            const auto value = gf(p);
            for (size_t r = 0; r < fes.getRows(); ++r)
            {
              for (size_t c = 0; c < fes.getColumns(); ++c)
                os << value(r, c) << ' ';
            }
            os << '\n';
          }
          else
            os << gf(p) << '\n';
        }
        os << '\n';
      }
  };

  /** @brief Reads matrix entries stored as independent MEDIT scalar vertex fields.
   * Higher-order spaces interpolate the vertex field, as required by SolAtVertices.
   */
  template <class FES>
    requires(
      FormLanguage::IsMatrixRange<typename FormLanguage::Traits<FES>::RangeType>::Value)
  class GridFunctionLoader<FileFormat::MEDIT, FES,
    Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
    : public GridFunctionLoaderBase<FES,
        Math::Vector<typename FormLanguage::Traits<FES>::ScalarType>>
  {
    public:
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;

      /// @brief CRTP or finite element base class.
      using Parent = GridFunctionLoaderBase<FES, Math::Vector<ScalarType>>;
      using Parent::Parent;
      void load(std::istream& is) override
      {
        auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t components = fes.getVectorDimension();
        std::string keyword;
        size_t version, dimension, vertices, fields;
        auto fail = [this, function = __func__] {
          Alert::MemberFunctionException(*this, function)
            << "Invalid MEDIT matrix vertex data." << Alert::Raise;
        };
        if (!(is >> keyword >> version) || keyword != "MeshVersionFormatted" ||
          version != 2)
          fail();
        if (!(is >> keyword >> dimension) || keyword != "Dimension" ||
          dimension != mesh.getSpaceDimension())
          fail();
        if (!(is >> keyword >> vertices >> fields) || keyword != "SolAtVertices" ||
          vertices != mesh.getVertexCount() || fields != components)
          fail();
        for (size_t c = 0; c < fields; ++c)
        {
          size_t type;
          if (!(is >> type) || type != MEDIT::SolutionType::Real)
            fail();
        }
        Variational::P1<typename FormLanguage::Traits<FES>::RangeType,
          Geometry::Mesh<Context::Local>>
          vertexSpace(mesh, fes.getRows(), fes.getColumns());
        Variational::GridFunction vertexField(vertexSpace);
        for (size_t v = 0; v < vertices; ++v)
        {
          const auto& dofs = vertexSpace.getDOFs(0, v);
          for (size_t c = 0; c < components; ++c)
          {
            if (!(is >> vertexField.getData().coeffRef(dofs[c])))
              fail();
          }
        }
        gf = vertexField;
      }
  };
}

#endif
