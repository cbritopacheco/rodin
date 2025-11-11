/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_MFEM_H
#define RODIN_IO_MFEM_H

#include <boost/bimap.hpp>
#include <boost/spirit/home/x3.hpp>
#include <ostream>

#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"
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
   * @brief Keywords used in MFEM mesh file format.
   *
   * These keywords identify different sections of an MFEM mesh file.
   */
  enum class Keyword
  {
    dimension, ///< Dimension section keyword
    elements,  ///< Elements section keyword  
    boundary,  ///< Boundary section keyword
    vertices   ///< Vertices section keyword
  };

  /**
   * @brief Converts a keyword enum to its string representation.
   * @param kw Keyword to convert
   * @return C-style string representation of the keyword
   */
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

  /**
   * @brief Stream output operator for MFEM keywords.
   * @param os Output stream
   * @param kw Keyword to output
   * @return Reference to the output stream
   */
  inline
  std::ostream& operator<<(std::ostream& os, Keyword kw)
  {
    os << toCharString(kw);
    return os;
  }

  /**
   * @brief Converts a C-style string to a keyword enum.
   * @param str String to convert
   * @return Optional keyword if conversion succeeds, empty otherwise
   */
  inline
  Optional<Keyword> toKeyword(const char* str)
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

  /**
   * @brief Version information for MFEM mesh format.
   *
   * Stores major and minor version numbers for the MFEM mesh file format.
   */
  struct MeshVersion
  {
    size_t major;  ///< Major version number
    size_t minor;  ///< Minor version number
  };

  /**
   * @brief Header information for MFEM mesh files.
   *
   * Contains the mesh type and version information parsed from the mesh file header.
   */
  struct MeshHeader
  {
    MeshType type;        ///< Type of mesh (LEGACY, NONCONFORMING, or NURBS)
    MeshVersion version;  ///< Version of the mesh format
  };

  /**
   * @brief Data ordering for grid function storage.
   *
   * Specifies how vector-valued grid function data is organized in memory.
   */
  enum Ordering
  {
    /// Node-major ordering: XXX..., YYY..., ZZZ... (all x-components, then all y-components, etc.)
    Nodes = 0,

    /// Vector dimension-major ordering: XYZ, XYZ, ... (interleaved components)
    VectorDimension = 1
  };

  /**
   * @brief Header information for MFEM grid function files.
   *
   * Contains metadata about the finite element collection, vector dimension,
   * and data ordering for grid functions.
   */
  struct GridFunctionHeader
  {
    std::string fec;     ///< Finite element collection name (e.g., "H1_2D_P1")
    size_t vdim;         ///< Vector dimension of the grid function
    Ordering ordering;   ///< Data ordering (Nodes or VectorDimension)
  };

  /**
   * @brief Converts MFEM geometry type to Rodin polytope type.
   * @param[in] t MFEM geometry type
   * @returns Optional Rodin polytope type, empty if conversion not supported
   *
   * Maps MFEM's internal geometry type enumeration to Rodin's polytope type system.
   */
  inline
  constexpr
  Optional<Rodin::Geometry::Polytope::Type> getGeometry(GeometryType t)
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
      case GeometryType::PRISM:
      {
        return Rodin::Geometry::Polytope::Type::Wedge;
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

  /**
   * @brief Converts Rodin polytope type to MFEM geometry type.
   * @param[in] t Rodin polytope type
   * @returns Optional MFEM geometry type, empty if conversion not supported
   *
   * Maps Rodin's polytope type system to MFEM's internal geometry type enumeration.
   */
  inline
  constexpr
  Optional<GeometryType> getGeometry(Geometry::Polytope::Type t)
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
      case Geometry::Polytope::Type::Wedge:
        return GeometryType::PRISM;
      default:
        return {};
    }
    assert(false);
    return {};
  }

  /**
   * @brief Parser for unsigned integers in MFEM format.
   * @internal
   *
   * Uses Boost.Spirit.X3 to parse unsigned integer values from input text.
   */
  struct ParseUnsignedInteger
  {
    /**
     * @brief Parses an unsigned integer from an iterator range.
     * @tparam Iterator Iterator type
     * @param[in] begin Start of input range
     * @param[in] end End of input range
     * @returns Optional unsigned integer if parsing succeeds, empty otherwise
     */
    template <class Iterator>
    inline
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

  /**
   * @brief Parser for vertex coordinates in MFEM format.
   * @internal
   *
   * Parses spatial coordinates for a single vertex from input text.
   */
  class ParseVertex
  {
    public:
      /**
       * @brief Constructs a vertex parser for the given spatial dimension.
       * @param[in] sdim Spatial dimension (2D or 3D)
       */
      ParseVertex(size_t sdim)
        : m_sdim(sdim)
      {}

      /**
       * @brief Parses vertex coordinates from an iterator range.
       * @tparam Iterator Iterator type
       * @param[in] begin Start of input range
       * @param[in] end End of input range
       * @returns Optional spatial point if parsing succeeds, empty otherwise
       */
      template <class Iterator>
      Optional<Math::SpatialPoint> operator()(Iterator begin, Iterator end) const
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::blank;
        using boost::spirit::x3::repeat;
        using boost::spirit::x3::double_;
        using boost::spirit::x3::_attr;

        size_t i = 0;
        Math::SpatialPoint res(m_sdim);
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

  /**
   * @brief Parser for geometry (element) information in MFEM format.
   * @internal
   *
   * Parses element data including attribute, geometry type, and vertex indices.
   */
  struct ParseGeometry
  {
    /**
     * @brief Parsed geometry data.
     */
    struct Data
    {
      Geometry::Attribute attribute;       ///< Element attribute (material ID)
      Geometry::Polytope::Type geometry;   ///< Element geometry type
      Array<Index> vertices;               ///< Vertex indices defining the element
    };

    /**
     * @brief Parses geometry data from an iterator range.
     * @tparam Iterator Iterator type
     * @param[in] begin Start of input range
     * @param[in] end End of input range
     * @returns Optional geometry data if parsing succeeds, empty otherwise
     */
    template <class Iterator>
    inline
    Optional<Data> operator()(Iterator begin, Iterator end) const
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

     res.vertices.resize(Geometry::Polytope::Traits(res.geometry).getVertexCount());
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

  /**
   * @brief Parser for empty lines.
   * @internal
   *
   * Checks if a line contains only whitespace.
   */
  struct ParseEmptyLine
  {
    /**
     * @brief Checks if the iterator range contains only whitespace.
     * @tparam Iterator Iterator type
     * @param[in] begin Start of input range
     * @param[in] end End of input range
     * @returns True if line is empty or contains only whitespace, false otherwise
     */
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

  /**
   * @brief Parser for empty lines or comment lines.
   * @internal
   *
   * Checks if a line is empty, contains only whitespace, or is a comment (starts with #).
   */
  struct ParseEmptyLineOrComment
  {
    /**
     * @brief Checks if the iterator range is empty, whitespace, or a comment.
     * @tparam Iterator Iterator type
     * @param[in] begin Start of input range
     * @param[in] end End of input range
     * @returns True if line is empty, whitespace, or comment, false otherwise
     */
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

  /**
   * @brief Parser for MFEM keywords.
   * @internal
   *
   * Parses keyword strings (alphabetic text) from MFEM format input.
   */
  struct ParseKeyword
  {
    /**
     * @brief Parses a keyword from an iterator range.
     * @tparam Iterator Iterator type
     * @param[in] begin Start of input range
     * @param[in] end End of input range
     * @returns Optional keyword string if parsing succeeds, empty otherwise
     */
    template <class Iterator>
    inline
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

  /**
   * @brief Parser for MFEM mesh file headers.
   * @internal
   *
   * Parses the mesh header line which contains format identifier and version.
   * Expected format: "MFEM mesh v1.0" or similar.
   */
  class ParseMeshHeader
  {
    public:
      /**
       * @brief Parses a mesh header from an iterator range.
       * @tparam Iterator Iterator type
       * @param[in] begin Start of input range
       * @param[in] end End of input range
       * @returns Optional mesh header if parsing succeeds, empty otherwise
       */
      template <class Iterator>
      inline
      Optional<MeshHeader> operator()(Iterator begin, Iterator end) const
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
  /**
   * @ingroup MeshLoaderSpecializations
   * @brief Specialization for loading sequential meshes in MFEM format.
   *
   * This loader reads mesh data from the MFEM mesh file format, a text-based
   * format used by the MFEM library. It supports various element types and
   * can handle 2D and 3D meshes.
   *
   * ## MFEM Format Structure
   * The MFEM mesh format consists of several sections:
   * - Header: "MFEM mesh v1.0" with version information
   * - dimension: Spatial dimension of the mesh
   * - elements: Element connectivity and attributes
   * - boundary: Boundary element information
   * - vertices: Vertex coordinates
   *
   * ## Usage Example
   * ```cpp
   * Mesh<Context::Local> mesh;
   * MeshLoader<FileFormat::MFEM, Context::Local> loader(mesh);
   * loader.load("mesh.mfem");
   * ```
   *
   * @see MeshPrinter
   */
  template <>
  class MeshLoader<IO::FileFormat::MFEM, Context::Local>
    : public MeshLoaderBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = MeshPrinterBase<ContextType>;

      /**
       * @brief Constructs an MFEM mesh loader for the given mesh.
       * @param[in,out] mesh Mesh object to populate with loaded data
       */
      MeshLoader(ObjectType& mesh)
        : MeshLoaderBase<Context::Local>(mesh)
      {}

      /**
       * @brief Loads mesh from an input stream.
       * @param[in] is Input stream containing MFEM mesh data
       *
       * Reads and parses the complete mesh data including header, dimension,
       * elements, boundaries, and vertices.
       */
      void load(std::istream& is) override;

      /**
       * @brief Reads the mesh file header.
       * @param[in] is Input stream
       *
       * Parses the MFEM header line (e.g., "MFEM mesh v1.0").
       */
      void readHeader(std::istream& is);

      /**
       * @brief Reads the mesh dimension.
       * @param[in] is Input stream
       *
       * Parses the "dimension" section to determine spatial dimension.
       */
      void readDimension(std::istream& is);

      /**
       * @brief Reads the complete mesh data.
       * @param[in] is Input stream
       *
       * Parses elements, boundaries, and vertices sections.
       */
      void readMesh(std::istream& is);

    private:
      size_t m_dimension;
      size_t m_spaceDimension;
      size_t m_currentLineNumber;
      MFEM::MeshHeader m_header;
      ObjectType::Builder m_build;
  };

  /**
   * @ingroup PrinterSpecializations
   * @brief Specialization for printing sequential meshes in MFEM format.
   *
   * This printer writes mesh data in the MFEM mesh file format, a text-based
   * format compatible with the MFEM library.
   *
   * ## Usage Example
   * ```cpp
   * const Mesh<Context::Local>& mesh = getMesh();
   * MeshPrinter<FileFormat::MFEM, Context::Local> printer(mesh);
   * std::ofstream file("output.mfem");
   * printer.print(file);
   * ```
   *
   * @see MeshLoader
   */
  template <>
  class MeshPrinter<FileFormat::MFEM, Context::Local>
    : public MeshPrinterBase<Context::Local>
  {
    public:
      using ContextType = Context::Local;

      using ObjectType = Geometry::Mesh<ContextType>;

      using Parent = MeshPrinterBase<ContextType>;

      /**
       * @brief Constructs an MFEM mesh printer for the given mesh.
       * @param[in] mesh Mesh object to write to output
       */
      MeshPrinter(const ObjectType& mesh)
        : MeshPrinterBase(mesh)
      {}

      /**
       * @brief Prints mesh to an output stream.
       * @param[in,out] os Output stream to write to
       *
       * Writes the complete mesh data in MFEM format including header,
       * dimension, elements, boundaries, and vertices.
       */
      void print(std::ostream& os) override;

      /**
       * @brief Prints the mesh file header.
       * @param[in,out] os Output stream
       */
      void printHeader(std::ostream& os);

      /**
       * @brief Prints the mesh dimension.
       * @param[in,out] os Output stream
       */
      void printDimension(std::ostream& os);

      /**
       * @brief Prints the mesh connectivity and vertex data.
       * @param[in,out] os Output stream
       */
      void printMesh(std::ostream& os);
  };

  /**
   * @brief Specialization for loading P1 grid functions from MFEM format.
   *
   * Loads finite element solution data for continuous Lagrange (P1) elements
   * from MFEM grid function files.
   *
   * @tparam Range Range type for the finite element space
   *
   * ## MFEM Grid Function Format
   * The format includes:
   * - FiniteElementSpace header
   * - FiniteElementCollection name (e.g., "H1_2D_P1")
   * - VDim: vector dimension
   * - Ordering: data layout (Nodes=0 or VectorDimension=1)
   * - Coefficient data values
   *
   * ## Usage Example
   * ```cpp
   * P1<Real> Vh(mesh);
   * GridFunction<P1<Real>> u(Vh);
   * GridFunctionLoader<FileFormat::MFEM, P1<Real>, Vector<Real>> loader(u);
   * loader.load("solution.gf");
   * ```
   *
   * @see GridFunctionPrinter
   */
  template <class Range>
  class GridFunctionLoader<
    FileFormat::MFEM,
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

      /**
       * @brief Constructs a grid function loader.
       * @param[in,out] gf Grid function to populate with loaded data
       */
      GridFunctionLoader(ObjectType& gf)
        : Parent(gf)
      {}

      /**
       * @brief Loads grid function from an input stream.
       * @param[in] is Input stream containing MFEM grid function data
       *
       * Parses the header and coefficient data, handling different data orderings.
       */
      void load(std::istream& is) override
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::blank;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::char_;

        MFEM::GridFunctionHeader header;
        const auto get_fec = [&](auto& ctx) { header.fec = _attr(ctx); };
        const auto get_vdim = [&](auto& ctx) { header.vdim = _attr(ctx); };
        const auto get_ordering = [&](auto& ctx) { header.ordering = static_cast<MFEM::Ordering>(_attr(ctx)); };

        std::string line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        auto it = line.begin();
        const auto pfes = boost::spirit::x3::string("FiniteElementSpace");
        const bool rfes = boost::spirit::x3::phrase_parse(it, line.end(), pfes, space);
        assert(it == line.end() && rfes);

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it = line.begin();
        const auto pfec = boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[get_fec];
        bool rfec = boost::spirit::x3::phrase_parse(it, line.end(), pfec, space);
        assert(it == line.end() && rfec);

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it = line.begin();
        const auto pvdim = boost::spirit::x3::string("VDim:") >> uint_[get_vdim];
        bool rvdim = boost::spirit::x3::phrase_parse(it, line.end(), pvdim, space);
        assert(it == line.end() && rvdim);

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it = line.begin();
        const auto pordering = boost::spirit::x3::string("Ordering:") >> uint_[get_ordering];
        bool rordering = boost::spirit::x3::phrase_parse(it, line.end(), pordering, space);
        assert(it == line.end() && rordering);

        auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        assert(header.vdim == fes.getVectorDimension());
        auto& data = gf.getData();
        if (data.size() > 0)
        {
          line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
          data.coeffRef(0) = std::stod(line);
          assert(data.size() >= 0);
          for (size_t i = 1; i < static_cast<size_t>(data.size()); i++)
            is >> data.coeffRef(i);
          if (header.ordering == MFEM::Ordering::VectorDimension)
            data.transposeInPlace();
        }
      }

    private:
      size_t m_dimension;
      size_t m_spaceDimension;
      size_t m_currentLineNumber;
  };

  /**
   * @brief Base class for printing P0 (piecewise constant) grid functions in MFEM format.
   *
   * Handles output of discontinuous finite element solutions on cells.
   *
   * @tparam Range Range type for the finite element space
   * @tparam Context Context type (e.g., Context::Local)
   * @tparam Data Data storage type
   *
   * @see GridFunctionPrinter
   */
  template <class Range, class Context, class Data>
  class GridFunctionPrinterBase<FileFormat::MFEM, Variational::P0<Range, Geometry::Mesh<Context>>, Data>
  : public Printer<Variational::GridFunction<Variational::P0<Range, Geometry::Mesh<Context>>, Data>>
  {
    public:
      using RangeType = Range;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using DataType = Data;

      using MeshType = Geometry::Mesh<Context>;

      using FESType = Variational::P0<Range, MeshType>;

      static constexpr FileFormat Format = FileFormat::MFEM;

      using ObjectType = Variational::GridFunction<FESType, DataType>;

      using Parent = Printer<ObjectType>;

      /**
       * @brief Constructs a grid function printer.
       * @param[in] gf Grid function to write to output
       */
      GridFunctionPrinterBase(const ObjectType& gf)
        : m_gf(gf)
      {}

      /**
       * @brief Prints grid function to an output stream.
       * @param[in,out] os Output stream
       *
       * Writes the finite element space header and coefficient data.
       */
      void print(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        os << "FiniteElementSpace\n"
           << "FiniteElementCollection: " << "L2_" << fes.getMesh().getDimension() << "D_P0\n"
           << "VDim: " << fes.getVectorDimension() << '\n'
           << "Ordering: " << MFEM::Ordering::Nodes
           << "\n\n";
        this->printData(os);
      }

      const ObjectType& getObject() const override
      {
        return m_gf.get();
      }

      /**
       * @brief Prints the coefficient data.
       * @param[in,out] os Output stream
       *
       * Must be implemented by derived classes.
       */
      virtual void printData(std::ostream& os) = 0;

    private:
      std::reference_wrapper<const ObjectType> m_gf;
  };


  /**
   * @brief Base class for printing P1 (continuous Lagrange) grid functions in MFEM format.
   *
   * Handles output of continuous finite element solutions on nodes.
   *
   * @tparam Range Range type for the finite element space
   * @tparam Context Context type (e.g., Context::Local)
   * @tparam Data Data storage type
   *
   * @see GridFunctionPrinter
   */
  template <class Range, class Context, class Data>
  class GridFunctionPrinterBase<FileFormat::MFEM, Variational::P1<Range, Geometry::Mesh<Context>>, Data>
  : public Printer<Variational::GridFunction<Variational::P1<Range, Geometry::Mesh<Context>>, Data>>
  {
    public:
      using RangeType = Range;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using DataType = Data;

      using MeshType = Geometry::Mesh<Context>;

      using FESType = Variational::P1<Range, MeshType>;

      static constexpr FileFormat Format = FileFormat::MFEM;

      using ObjectType = Variational::GridFunction<FESType, DataType>;

      using Parent = Printer<ObjectType>;

      /**
       * @brief Constructs a grid function printer.
       * @param[in] gf Grid function to write to output
       */
      GridFunctionPrinterBase(const ObjectType& gf)
        : m_gf(gf)
      {}

      /**
       * @brief Prints grid function to an output stream.
       * @param[in,out] os Output stream
       *
       * Writes the finite element space header and coefficient data.
       */
      void print(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();
        os << "FiniteElementSpace\n"
           << "FiniteElementCollection: " << "H1_" << fes.getMesh().getDimension() << "D_P1\n"
           << "VDim: " << fes.getVectorDimension() << '\n'
           << "Ordering: " << MFEM::Ordering::Nodes
           << "\n\n";
        this->printData(os);
      }

      const ObjectType& getObject() const override
      {
        return m_gf.get();
      }

      /**
       * @brief Prints the coefficient data.
       * @param[in,out] os Output stream
       *
       * Must be implemented by derived classes.
       */
      virtual void printData(std::ostream& os) = 0;

    private:
      std::reference_wrapper<const ObjectType> m_gf;
  };

  /**
   * @brief Final specialization for printing grid functions with vector data in MFEM format.
   *
   * Implements the complete printer for grid functions stored as Math::Vector objects.
   *
   * @tparam FES Finite element space type
   * @tparam Scalar Scalar type for the vector data
   *
   * ## Usage Example
   * ```cpp
   * P1<Real> Vh(mesh);
   * GridFunction<P1<Real>> u(Vh);
   * // ... compute solution ...
   * GridFunctionPrinter<FileFormat::MFEM, P1<Real>, Vector<Real>> printer(u);
   * std::ofstream file("solution.gf");
   * printer.print(file);
   * ```
   */
  template <class FES, class Scalar>
  class GridFunctionPrinter<FileFormat::MFEM, FES, Math::Vector<Scalar>> final
    : public GridFunctionPrinterBase<FileFormat::MFEM, FES, Math::Vector<Scalar>>
  {
    public:
      using DataType = Math::Vector<Scalar>;

      using ObjectType = Variational::GridFunction<FES, DataType>;

      using Parent = GridFunctionPrinterBase<FileFormat::MFEM, FES, DataType>;

      /**
       * @brief Constructs a grid function printer.
       * @param[in] gf Grid function to write to output
       */
      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      /**
       * @brief Prints the coefficient data values.
       * @param[in,out] os Output stream
       *
       * Writes all coefficient values, one per line.
       */
      void printData(std::ostream& os) override
      {
        const auto& gf = this->getObject();
        const auto& vec = gf.getData();
        const auto* data = vec.data();
        assert(vec.size() >= 0);
        for (size_t i = 0; i < static_cast<size_t>(vec.size()); i++)
          os << data[i] << '\n';
      }
  };
}

#endif
