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
#include "Rodin/Variational/H1/H1.h"

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
      size_t m_currentLineNumber = 0;
  };

  /**
   * @brief Specialization for loading H1 grid functions from MFEM format.
   *
   * Loads finite element solution data for H1-conforming Lagrange elements
   * of arbitrary degree from MFEM grid function files. Handles DOF reordering
   * between MFEM's ordering and Rodin's internal ordering.
   *
   * @tparam K Polynomial degree
   * @tparam Range Range type for the finite element space
   *
   * ## MFEM Grid Function Format
   * The format includes:
   * - FiniteElementSpace header
   * - FiniteElementCollection name (e.g., "H1_2D_P2" for degree 2 in 2D)
   * - VDim: vector dimension
   * - Ordering: data layout (Nodes=0 or VectorDimension=1)
   * - Coefficient data values (ordered by MFEM convention)
   *
   * MFEM DOF ordering: vertices -> edge interiors -> face interiors -> element interiors
   *
   * ## Usage Example
   * ```cpp
   * H1 fes(std::integral_constant<size_t, 2>{}, mesh);
   * GridFunction gf(fes);
   * GridFunctionLoader<FileFormat::MFEM, H1<2, Real>, Vector<Real>> loader(gf);
   * loader.load("solution.gf");
   * ```
   *
   * @see GridFunctionPrinter
   */
  template <size_t K, class Range>
  class GridFunctionLoader<
    FileFormat::MFEM,
    Variational::H1<K, Range, Geometry::Mesh<Context::Local>>,
    Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
    : public GridFunctionLoaderBase<
        Variational::H1<K, Range, Geometry::Mesh<Context::Local>>,
        Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
  {
    public:
      using FESType = Variational::H1<K, Range, Geometry::Mesh<Context::Local>>;

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
       * Parses the header and coefficient data, mapping from MFEM ordering to Rodin ordering.
       */
      void load(std::istream& is) override
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::char_;

        MFEM::GridFunctionHeader header;
        const auto get_fec      = [&](auto& ctx) { header.fec      = _attr(ctx); };
        const auto get_vdim     = [&](auto& ctx) { header.vdim     = _attr(ctx); };
        const auto get_ordering = [&](auto& ctx) {
          header.ordering = static_cast<MFEM::Ordering>(_attr(ctx));
        };

        // --- parse header ------------------------------------------------------
        std::string line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        auto it = line.begin();
        const auto pfes = boost::spirit::x3::string("FiniteElementSpace");
        const bool rfes = boost::spirit::x3::phrase_parse(it, line.end(), pfes, space);
        assert(it == line.end() && rfes);

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it = line.begin();
        const auto pfec =
          boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[get_fec];
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

        // --- basic info --------------------------------------------------------
        auto& gf   = this->getObject();
        auto& data = gf.getData();
        const auto& fes  = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        const size_t vdim = fes.getVectorDimension();
        const size_t D = mesh.getDimension();
        const size_t scalarSize = fes.getSize() / vdim;

        assert(header.vdim == vdim);
        if (data.size() == 0)
          return;

        assert(static_cast<size_t>(data.size()) == vdim * scalarSize);

        // --- build scalar dof -> MFEM scalar position map ----------------------
        // MFEM scalar ordering: vertices -> edges (oriented) -> faces -> elements.
        std::vector<size_t> mfem_pos(scalarSize, static_cast<size_t>(-1));
        size_t p = 0;

        // 0. vertex -> scalar dof map (as in printer)
        const size_t nVertices = mesh.getConnectivity().getCount(0);
        std::vector<Index> vertexScalarDof(nVertices);
        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
        {
          const auto& vdofs = fes.getDOFs(0, v);
          assert(vdofs.size() == 1); // H1 nodal: one scalar dof per vertex
          vertexScalarDof[v] = vdofs(0);
        }

        auto register_dof = [&](Index d)
        {
          const size_t idx = static_cast<size_t>(d);
          if (mfem_pos[idx] == static_cast<size_t>(-1))
          {
            mfem_pos[idx] = p++;
          }
        };

        // 1. vertices
        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
        {
          register_dof(vertexScalarDof[v]);
        }

        // 2. edges with MFEM orientation (same logic as printer)
        if (D >= 1)
        {
          const auto& conn10 = mesh.getConnectivity().getIncidence(1, 0);
          const size_t nEdges = mesh.getConnectivity().getCount(1);

          for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
          {
            const auto& edgeVerts = conn10[e];
            assert(edgeVerts.size() == 2);

            Index v0 = edgeVerts[0];
            Index v1 = edgeVerts[1];

            Index vmin = std::min(v0, v1);
            Index vmax = std::max(v0, v1);

            Index vminDof = vertexScalarDof[vmin];
            Index vmaxDof = vertexScalarDof[vmax];

            const auto& edofs = fes.getDOFs(1, e);

            std::vector<Index> interior;
            interior.reserve(edofs.size());
            for (Index k = 0; k < static_cast<Index>(edofs.size()); ++k)
            {
              Index d = edofs(k);
              if (d != vminDof && d != vmaxDof)
                interior.push_back(d);
            }

            // If local orientation is opposite to MFEM's (low->high), reverse
            if (v0 > v1)
              std::reverse(interior.begin(), interior.end());

            for (Index d : interior)
              register_dof(d);
          }
        }

        // 3. faces (we assume H1 face interior dofs are already stored in
        //    MFEM-consistent local order in fes.getDOFs(faceDim, f))
        if (D >= 2)
        {
          const size_t faceDim   = (D == 3) ? 2 : (D - 1);
          const size_t faceCount = mesh.getConnectivity().getCount(faceDim);
          for (Index f = 0; f < static_cast<Index>(faceCount); ++f)
          {
            const auto& fdofs = fes.getDOFs(faceDim, f);
            for (Index k = 0; k < static_cast<Index>(fdofs.size()); ++k)
              register_dof(fdofs(k));
          }
        }

        // 4. element interiors
        const size_t nCells = mesh.getConnectivity().getCount(D);
        for (Index c = 0; c < static_cast<Index>(nCells); ++c)
        {
          const auto& cdofs = fes.getDOFs(D, c);
          for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
            register_dof(cdofs(k));
        }

        assert(p == scalarSize);
        for (size_t i = 0; i < scalarSize; ++i)
        {
          assert(mfem_pos[i] != static_cast<size_t>(-1));
        }

        // --- read MFEM data (all components, all scalar dofs) -------------------
        std::vector<ScalarType> mfem_data(static_cast<size_t>(data.size()));

        // First value: line may contain only that number
        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        mfem_data[0] = static_cast<ScalarType>(std::stod(line));

        for (size_t i = 1; i < mfem_data.size(); ++i)
        {
          ScalarType val;
          is >> val;
          mfem_data[i] = val;
        }

        // --- map MFEM data -> Rodin internal storage (always Nodes ordering) ----
        // Internal layout: data[s + c*scalarSize] is component c of scalar dof s.
        // MFEM layout:
        //   Nodes:           for each scalar position p, components c = 0..vdim-1
        //   VectorDimension: for each component c, positions p = 0..scalarSize-1
        for (size_t c = 0; c < vdim; ++c)
        {
          for (size_t s = 0; s < scalarSize; ++s)
          {
            const size_t p_scalar = mfem_pos[s];

            size_t mfem_index = 0;
            if (header.ordering == MFEM::Ordering::Nodes)
            {
              // Node-major: [ (p0,c0..c_{vdim-1}), (p1,c0..), ... ]
              mfem_index = p_scalar * vdim + c;
            }
            else // MFEM::Ordering::VectorDimension
            {
              // Component-major: [ (c0,p0..p_{scalarSize-1}), (c1,p0..), ... ]
              mfem_index = c * scalarSize + p_scalar;
            }

            data.coeffRef(s + c * scalarSize) = mfem_data[mfem_index];
          }
        }
      }

    private:
      size_t m_currentLineNumber = 0;
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
   * @brief Base class for printing H1 (continuous Lagrange) grid functions in MFEM format.
   *
   * Handles output of continuous finite element solutions on nodes for arbitrary degree H1 spaces.
   *
   * @tparam K Polynomial degree
   * @tparam Range Range type for the finite element space
   * @tparam Context Context type (e.g., Context::Local)
   * @tparam Data Data storage type
   *
   * @see GridFunctionPrinter
   */
  template <size_t K, class Range, class Context, class Data>
  class GridFunctionPrinterBase<FileFormat::MFEM, Variational::H1<K, Range, Geometry::Mesh<Context>>, Data>
  : public Printer<Variational::GridFunction<Variational::H1<K, Range, Geometry::Mesh<Context>>, Data>>
  {
    public:
      using RangeType = Range;

      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      using DataType = Data;

      using MeshType = Geometry::Mesh<Context>;

      using FESType = Variational::H1<K, Range, MeshType>;

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
           << "FiniteElementCollection: " << "H1_" << fes.getMesh().getDimension() << "D_P" << K << "\n"
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

  /**
   * @brief Specialization for printing H1 grid functions with vector data in MFEM format.
   *
   * Handles DOF reordering between Rodin's internal ordering and MFEM's expected ordering.
   * MFEM orders DOFs as: vertices -> edge interiors -> face interiors -> element interiors
   *
   * @tparam K Polynomial degree
   * @tparam Range Range type for the finite element space
   * @tparam Scalar Scalar type for the vector data
   */
  template <size_t K, class Range, class Scalar>
  class GridFunctionPrinter<
      FileFormat::MFEM,
      Variational::H1<K, Range, Geometry::Mesh<Context::Local>>,
      Math::Vector<Scalar>> final
    : public GridFunctionPrinterBase<
          FileFormat::MFEM,
          Variational::H1<K, Range, Geometry::Mesh<Context::Local>>,
          Math::Vector<Scalar>>
  {
  public:
    using FESType   = Variational::H1<K, Range, Geometry::Mesh<Context::Local>>;
    using DataType  = Math::Vector<Scalar>;
    using ObjectType = Variational::GridFunction<FESType, DataType>;
    using Parent     = GridFunctionPrinterBase<FileFormat::MFEM, FESType, DataType>;

    GridFunctionPrinter(const ObjectType& gf)
      : Parent(gf)
    {}

    void printData(std::ostream& os) override
    {
      const auto& gf   = this->getObject();
      const auto& fes  = gf.getFiniteElementSpace();
      const auto& mesh = fes.getMesh();
      const auto& data = gf.getData();

      const size_t vdim       = fes.getVectorDimension();
      const size_t D          = mesh.getDimension();
      const size_t sdim       = mesh.getSpaceDimension();
      const size_t scalarSize = fes.getSize() / vdim;

      const std::streamsize old_prec = os.precision();
      const std::ios::fmtflags old_flags = os.flags();

      os << std::setprecision(std::numeric_limits<Scalar>::digits10 + 2);
      os.setf(std::ios::scientific, std::ios::floatfield);

      // Track which scalar DOFs are already written
      std::vector<uint8_t> written(scalarSize, uint8_t(0));

      // Vertex -> scalar DOF (H1: one per vertex)
      const size_t nVertices = mesh.getConnectivity().getCount(0);
      std::vector<Index> vertexScalarDof(nVertices);
      for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
      {
        const auto& vdofs = fes.getDOFs(0, v);
        assert(vdofs.size() == 1 && "H1 vertex should have exactly one scalar dof.");
        vertexScalarDof[v] = vdofs(0);
      }

      // Emit coefficient(s) for a *scalar DOF index* (original Rodin dof)
      auto emit_scalar_dof = [&](Index rodin_dof)
      {
        if (written[rodin_dof])
          return;
        for (size_t c = 0; c < vdim; ++c)
          os << data.coeffRef(rodin_dof + c * scalarSize) << '\n';
        written[rodin_dof] = uint8_t(1);
      };

      // Emit coefficient(s) for a value returned by gf(..)
      auto emit_value = [&](const auto& val)
      {
        if constexpr (std::is_same_v<Range, Scalar>)
        {
          os << val << '\n';
        }
        else
        {
          static_assert(std::is_same_v<Range, Math::Vector<Scalar>>,
                        "Range must be Scalar or Math::Vector<Scalar>.");
          for (size_t c = 0; c < vdim; ++c)
            os << val[c] << '\n';
        }
      };

      // -----------------------------------------------------------------------
      // 1. Vertices
      // -----------------------------------------------------------------------
      for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
        emit_scalar_dof(vertexScalarDof[v]);

      // -----------------------------------------------------------------------
      // 2. Edges: interior DOFs, oriented low-vertex -> high-vertex
      // -----------------------------------------------------------------------
      if (D >= 1)
      {
        const auto& conn10  = mesh.getConnectivity().getIncidence(1, 0);
        const size_t nEdges = mesh.getConnectivity().getCount(1);

        for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
        {
          const auto& edgeVerts = conn10[e];
          assert(edgeVerts.size() == 2 && "Segment should have 2 vertices.");

          const Index v0 = edgeVerts[0];
          const Index v1 = edgeVerts[1];

          const Index vmin = std::min(v0, v1);
          const Index vmax = std::max(v0, v1);

          const Index vminDof = vertexScalarDof[vmin];
          const Index vmaxDof = vertexScalarDof[vmax];

          const auto& edofs = fes.getDOFs(1, e);

          std::vector<Index> interior;
          interior.reserve(edofs.size());
          for (Index k = 0; k < static_cast<Index>(edofs.size()); ++k)
          {
            Index d = edofs(k);
            if (d != vminDof && d != vmaxDof)
              interior.push_back(d);
          }

          // Reverse if our local orientation is opposite to MFEM's
          if (v0 > v1)
            std::reverse(interior.begin(), interior.end());

          for (Index d : interior)
            emit_scalar_dof(d);
        }
      }

      // -----------------------------------------------------------------------
      // 3. Faces
      //    - D == 2: faces are 1D edges -> handled by edge block
      //    - D == 3: 2D faces; we special-case triangular faces to use MFEM
      //              H1_TriangleElement interior nodes.
      // -----------------------------------------------------------------------
      if (D >= 2)
      {
        const size_t faceDim   = (D == 3) ? 2 : (D - 1);
        const size_t faceCount = mesh.getConnectivity().getCount(faceDim);

        if (D == 3)
        {
          const auto& conn20 = mesh.getConnectivity().getIncidence(2, 0);
          const auto& cp     = Variational::GLL01<K>::getNodes();

          Math::SpatialPoint xref(2);
          Math::SpatialPoint xphys(sdim);

          for (Index f = 0; f < static_cast<Index>(faceCount); ++f)
          {
            const auto faceGeom = mesh.getGeometry(2, f);

            if (faceGeom == Geometry::Polytope::Type::Triangle)
            {
              const auto& fverts = conn20[f];
              assert(fverts.size() == 3 && "Triangle face must have 3 vertices.");

              const Index v0 = fverts[0];
              const Index v1 = fverts[1];
              const Index v2 = fverts[2];

              const auto X0 = mesh.getVertexCoordinates(v0);
              const auto X1 = mesh.getVertexCoordinates(v1);
              const auto X2 = mesh.getVertexCoordinates(v2);

              // MFEM H1_TriangleElement interior nodes on this face
              for (int j = 1; j < static_cast<int>(K); ++j)
              {
                for (int i = 1; i + j < static_cast<int>(K); ++i)
                {
                  const Real ci = cp[i];
                  const Real cj = cp[j];
                  const Real ck = cp[K - i - j];
                  const Real w  = ci + cj + ck;

                  const Real L2 = ci / w;
                  const Real L3 = cj / w;
                  const Real L1 = Real(1) - L2 - L3;

                  // reference (2D) coordinates on the face
                  xref.x() = L2;
                  xref.y() = L3;

                  // physical coordinates in R^sdim
                  xphys.x() = L1 * X0.x() + L2 * X1.x() + L3 * X2.x();
                  xphys.y() = L1 * X0.y() + L2 * X1.y() + L3 * X2.y();
                  if (sdim == 3)
                    xphys.z() = L1 * X0.z() + L2 * X1.z() + L3 * X2.z();

                  const auto  faceIt = mesh.getFace(f);
                  const auto& face   = *faceIt;
                  const Geometry::Point p(face, xref, xphys);

                  decltype(auto) val = gf(p);
                  emit_value(val);
                }
              }
            }
            else
            {
              // Non-triangular faces (e.g. quads): fall back to DOF-based ordering
              const auto& fdofs = fes.getDOFs(faceDim, f);
              for (Index k = 0; k < static_cast<Index>(fdofs.size()); ++k)
                emit_scalar_dof(fdofs(k));
            }
          }
        }
        else
        {
          // D == 2 or other: keep DOF-based behavior for faces
          for (Index f = 0; f < static_cast<Index>(faceCount); ++f)
          {
            const auto& fdofs = fes.getDOFs(faceDim, f);
            for (Index k = 0; k < static_cast<Index>(fdofs.size()); ++k)
              emit_scalar_dof(fdofs(k));
          }
        }
      }

      // -----------------------------------------------------------------------
      // 4. Element interiors
      // -----------------------------------------------------------------------
      const size_t nCells = mesh.getConnectivity().getCount(D);

      if (D == 2)
      {
        // 2D: triangles (possibly embedded in R^2 or R^3)
        const auto& conn20 = mesh.getConnectivity().getIncidence(2, 0);
        const auto& cp     = Variational::GLL01<K>::getNodes();

        Math::SpatialPoint xref(2);
        Math::SpatialPoint xphys(sdim);

        for (Index c = 0; c < static_cast<Index>(nCells); ++c)
        {
          const auto geom = mesh.getGeometry(2, c);
          if (geom != Geometry::Polytope::Type::Triangle)
          {
            // Fallback: non-triangle cells in 2D
            const auto& cdofs = fes.getDOFs(D, c);
            for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
              emit_scalar_dof(cdofs(k));
            continue;
          }

          const auto  cellIt = mesh.getCell(c);
          const auto& cell   = *cellIt;
          const auto& triVer = conn20[c];

          const Index v0 = triVer[0];
          const Index v1 = triVer[1];
          const Index v2 = triVer[2];

          const auto X0 = mesh.getVertexCoordinates(v0);
          const auto X1 = mesh.getVertexCoordinates(v1);
          const auto X2 = mesh.getVertexCoordinates(v2);

          // MFEM H1_TriangleElement interior nodes
          for (int j = 1; j < static_cast<int>(K); ++j)
          {
            for (int i = 1; i + j < static_cast<int>(K); ++i)
            {
              const Real ci = cp[i];
              const Real cj = cp[j];
              const Real ck = cp[K - i - j];
              const Real w  = ci + cj + ck;

              const Real L2 = ci / w;
              const Real L3 = cj / w;
              const Real L1 = Real(1) - L2 - L3;

              // reference coords (2D)
              xref.x() = L2;
              xref.y() = L3;

              // physical coords
              xphys.x() = L1 * X0.x() + L2 * X1.x() + L3 * X2.x();
              xphys.y() = L1 * X0.y() + L2 * X1.y() + L3 * X2.y();
              if (sdim == 3)
                xphys.z() = L1 * X0.z() + L2 * X1.z() + L3 * X2.z();

              const Geometry::Point p(cell, xref, xphys);
              decltype(auto) val = gf(p);
              emit_value(val);
            }
          }
        }
      }
      else if (D == 3)
      {
        const auto& conn30 = mesh.getConnectivity().getIncidence(3, 0);
        const auto& cp     = Variational::GLL01<K>::getNodes();

        Math::SpatialPoint xref(3);
        Math::SpatialPoint xphys(sdim);
        Math::SpatialPoint Xb(sdim);
        Math::SpatialPoint Xt(sdim);

        for (Index c = 0; c < static_cast<Index>(nCells); ++c)
        {
          const auto  cellIt = mesh.getCell(c);
          const auto& cell   = *cellIt;
          const auto& cellVer = conn30[c];

          const auto geom = mesh.getGeometry(3, c);

          // -------------------------------------------------------------
          // Tetrahedron: use MFEM H1_TetrahedronElement interior nodes
          // -------------------------------------------------------------
          if (geom == Geometry::Polytope::Type::Tetrahedron)
          {
            assert(cellVer.size() == 4 && "Tetrahedron must have 4 vertices.");

            const Index v0 = cellVer[0];
            const Index v1 = cellVer[1];
            const Index v2 = cellVer[2];
            const Index v3 = cellVer[3];

            const auto X0 = mesh.getVertexCoordinates(v0);
            const auto X1 = mesh.getVertexCoordinates(v1);
            const auto X2 = mesh.getVertexCoordinates(v2);
            const auto X3 = mesh.getVertexCoordinates(v3);

            // MFEM interior:
            // for (k=1..K-1) for (j=1..K-1-k) for (i=1..K-1-j-k)
            for (int k = 1; k < static_cast<int>(K); ++k)
            {
              for (int j = 1; j + k < static_cast<int>(K); ++j)
              {
                for (int i = 1; i + j + k < static_cast<int>(K); ++i)
                {
                  const Real ci = cp[i];
                  const Real cj = cp[j];
                  const Real ck = cp[k];
                  const Real cl = cp[K - i - j - k];
                  const Real w  = ci + cj + ck + cl;

                  const Real L2 = ci / w; // v1
                  const Real L3 = cj / w; // v2
                  const Real L4 = ck / w; // v3
                  const Real L1 = Real(1) - L2 - L3 - L4; // v0

                  // reference coordinates on the reference tetrahedron
                  xref.x() = L2;
                  xref.y() = L3;
                  xref.z() = L4;

                  // physical coordinates
                  xphys.x() = L1 * X0.x() + L2 * X1.x() + L3 * X2.x() + L4 * X3.x();
                  xphys.y() = L1 * X0.y() + L2 * X1.y() + L3 * X2.y() + L4 * X3.y();
                  xphys.z() = L1 * X0.z() + L2 * X1.z() + L3 * X2.z() + L4 * X3.z();

                  const Geometry::Point p(cell, xref, xphys);
                  decltype(auto) val = gf(p);
                  emit_value(val);
                }
              }
            }
          }
          // -------------------------------------------------------------
          // Wedge / prism: use MFEM H1_WedgeElement interior nodes
          // (triangle-interior × segment-interior)
          // -------------------------------------------------------------
          else if (geom == Geometry::Polytope::Type::Wedge)
          {
            assert(cellVer.size() == 6 && "Wedge must have 6 vertices.");

            // Assume bottom: (v0,v1,v2), top: (v3,v4,v5)
            const Index v0 = cellVer[0];
            const Index v1 = cellVer[1];
            const Index v2 = cellVer[2];
            const Index v3 = cellVer[3];
            const Index v4 = cellVer[4];
            const Index v5 = cellVer[5];

            const auto V0 = mesh.getVertexCoordinates(v0);
            const auto V1 = mesh.getVertexCoordinates(v1);
            const auto V2 = mesh.getVertexCoordinates(v2);
            const auto V3 = mesh.getVertexCoordinates(v3);
            const auto V4 = mesh.getVertexCoordinates(v4);
            const auto V5 = mesh.getVertexCoordinates(v5);

            // Precompute vertex coords into Xb/Xt via barycentric combination
            // Interior wedge nodes: tensor product of:
            //   - triangle-interior (i,j) on base tri (0,1,2)
            //   - segment-interior cp[k], k=1..K-1
            for (int k = 1; k < static_cast<int>(K); ++k)
            {
              const Real s           = cp[k];           // segment param in [0,1]
              const Real one_minus_s = Real(1) - s;

              for (int j = 1; j < static_cast<int>(K); ++j)
              {
                for (int i = 1; i + j < static_cast<int>(K); ++i)
                {
                  const Real ci = cp[i];
                  const Real cj = cp[j];
                  const Real ck = cp[K - i - j];
                  const Real w_tri = ci + cj + ck;

                  const Real L2 = ci / w_tri;
                  const Real L3 = cj / w_tri;
                  const Real L1 = Real(1) - L2 - L3;

                  // reference coordinates (r,s,t): (L2,L3,s)
                  xref.x() = L2;
                  xref.y() = L3;
                  xref.z() = s;

                  // bottom triangle physical point Xb
                  Xb.x() = L1 * V0.x() + L2 * V1.x() + L3 * V2.x();
                  Xb.y() = L1 * V0.y() + L2 * V1.y() + L3 * V2.y();
                  Xb.z() = L1 * V0.z() + L2 * V1.z() + L3 * V2.z();

                  // top triangle physical point Xt
                  Xt.x() = L1 * V3.x() + L2 * V4.x() + L3 * V5.x();
                  Xt.y() = L1 * V3.y() + L2 * V4.y() + L3 * V5.y();
                  Xt.z() = L1 * V3.z() + L2 * V4.z() + L3 * V5.z();

                  // interpolate along the segment
                  xphys.x() = one_minus_s * Xb.x() + s * Xt.x();
                  xphys.y() = one_minus_s * Xb.y() + s * Xt.y();
                  xphys.z() = one_minus_s * Xb.z() + s * Xt.z();

                  const Geometry::Point p(cell, xref, xphys);
                  decltype(auto) val = gf(p);
                  emit_value(val);
                }
              }
            }
          }
          // -------------------------------------------------------------
          // Other 3D cell types: fallback DOF-based
          // -------------------------------------------------------------
          else
          {
            const auto& cdofs = fes.getDOFs(D, c);
            for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
              emit_scalar_dof(cdofs(k));
          }
        }
      }
      else
      {
        // Other dimensions: fallback DOF-based
        for (Index c = 0; c < static_cast<Index>(nCells); ++c)
        {
          const auto& cdofs = fes.getDOFs(D, c);
          for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
            emit_scalar_dof(cdofs(k));
        }
      }

      os.precision(old_prec);
      os.flags(old_flags);
    }
  };
}

#endif
