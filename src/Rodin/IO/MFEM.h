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
#include <optional>

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

  // ----- MFEM-compatible nodal sets and Vandermonde matrices ------------------

  /**
   * @brief MFEM-style nodal set on the reference triangle (0,0)-(1,0)-(0,1).
   *
   * Reproduces the node placement of mfem::H1_TriangleElement using
   * poly1d.ClosedPoints(p, ...) = GLL01<K>::getNodes() in [0,1]:
   *
   *  - vertices:
   *      (cp[0], cp[0]), (cp[p], cp[0]), (cp[0], cp[p])
   *  - edge nodes: placed along each edge with cp[i]
   *  - interior: barycentric from cp[i]/w, cp[j]/w, cp[p-i-j]/w
   *
   * The ordering matches MFEM's construction: vertices, edges, interior.
   */
  template <size_t K>
  class TriangleNodes
  {
    public:
      static constexpr size_t Count = (K + 1) * (K + 2) / 2;

      static const std::vector<Math::SpatialPoint>& getNodes()
      {
        static thread_local std::vector<Math::SpatialPoint> s_nodes;

        if (s_nodes.empty())
        {
          s_nodes.reserve(Count);
          s_nodes.clear();

          const auto& cp = Variational::GLL01<K>::getNodes(); // assumed in [0,1], as in MFEM
          const int p = static_cast<int>(K);

          // Vertices
          s_nodes.emplace_back(Math::SpatialPoint{{ cp[0], cp[0]          }}); // (0,0)
          s_nodes.emplace_back(Math::SpatialPoint{{ cp[p], cp[0]          }}); // (1,0)
          s_nodes.emplace_back(Math::SpatialPoint{{ cp[0], cp[p]          }}); // (0,1)

          // Edge nodes (same loops/order as MFEM)
          // Edge (0,1): (cp[i], cp[0])
          for (int i = 1; i < p; ++i)
            s_nodes.emplace_back(Math::SpatialPoint{{ cp[i], cp[0] }});

          // Edge (1,2): (cp[p-i], cp[i])
          for (int i = 1; i < p; ++i)
            s_nodes.emplace_back(Math::SpatialPoint{{ cp[p - i], cp[i] }});

          // Edge (2,0): (cp[0], cp[p-i])
          for (int i = 1; i < p; ++i)
            s_nodes.emplace_back(Math::SpatialPoint{{ cp[0], cp[p - i] }});

          // Interior nodes
          for (int j = 1; j < p; ++j)
          {
            for (int i = 1; i + j < p; ++i)
            {
              const Real ci = cp[i];
              const Real cj = cp[j];
              const Real ck = cp[p - i - j];
              const Real w  = ci + cj + ck;

              const Real x = ci / w;
              const Real y = cj / w;

              s_nodes.emplace_back(Math::SpatialPoint{{ x, y }});
            }
          }

          assert(s_nodes.size() == Count);
        }

        return s_nodes;
      }
  };

  /**
   * @brief Vandermonde matrix on the triangle using MFEM's nodal set.
   *
   * V_MFEM(i,j) = ψ_j(x_i), where ψ_j are Dubiner modes and x_i are MFEM nodes.
   * This is the analogue of VandermondeTriangle<K> but with MfemTriangleNodes<K>.
   */
  template <size_t K>
  class VandermondeTriangle
  {
    public:
      static const Math::Matrix<Real>& getMatrix()
      {
        static thread_local Math::Matrix<Real> s_vandermonde;

        constexpr size_t N = TriangleNodes<K>::Count;

        if (s_vandermonde.size() == 0)
        {
          const auto& nodes = TriangleNodes<K>::getNodes();
          s_vandermonde.resize(N, N);

          size_t mode_idx = 0;
          Rodin::Utility::ForIndex<K + 1>(
            [&](auto p_idx)
            {
              constexpr size_t P = p_idx.value;
              Rodin::Utility::ForIndex<K + 1 - P>(
                [&](auto q_idx)
                {
                  constexpr size_t Q = q_idx.value;
                  for (size_t node_idx = 0; node_idx < N; ++node_idx)
                  {
                    const auto& pt = nodes[node_idx];
                    const Real x = pt.x();
                    const Real y = pt.y();

                    Real r, s;
                    Variational::DubinerTriangle<K>::getCollapsed(r, s, x, y);

                    Variational::DubinerTriangle<K>::template getBasis<P, Q>(
                      s_vandermonde(node_idx, mode_idx), r, s);
                  }
                  ++mode_idx;
                });
            });
        }

        return s_vandermonde;
      }

      static const Math::Matrix<Real>& getInverse()
      {
        static thread_local Math::Matrix<Real> s_inv;

        if (s_inv.size() == 0)
        {
          const auto& V = VandermondeTriangle<K>::getMatrix();
          Eigen::BDCSVD<Math::Matrix<Real>> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
          const Math::Matrix<Real> I = Math::Matrix<Real>::Identity(V.rows(), V.cols());
          s_inv = svd.solve(I);
        }

        return s_inv;
      }
  };

  // ---------------------------------------------------------------------------
  // MFEM-style nodal set and Vandermonde for the tetrahedron
  // ---------------------------------------------------------------------------

  /**
   * @brief MFEM-style nodal set on the reference tetrahedron
   * (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1).
   *
   * Reproduces node placement of mfem::H1_TetrahedronElement using
   * cp = poly1d.ClosedPoints(p, ...) = GLL01<K>::getNodes() in [0,1].
   *
   * Ordering matches MFEM's construction:
   *  - vertices
   *  - edges
   *  - faces
   *  - interior
   */
  template <size_t K>
  class TetrahedronNodes
  {
    public:
      static constexpr size_t Count = (K + 1) * (K + 2) * (K + 3) / 6;

      static const std::vector<Math::SpatialPoint>& getNodes()
      {
        static thread_local std::vector<Math::SpatialPoint> s_nodes;

        if (s_nodes.empty())
        {
          s_nodes.reserve(Count);
          s_nodes.clear();

          const auto& cp = Variational::GLL01<K>::getNodes(); // assumed in [0,1]
          const int p = static_cast<int>(K);

          // Vertices
          s_nodes.emplace_back(Math::SpatialPoint{{ cp[0], cp[0], cp[0] }}); // (0,0,0)
          s_nodes.emplace_back(Math::SpatialPoint{{ cp[p], cp[0], cp[0] }}); // (1,0,0)
          s_nodes.emplace_back(Math::SpatialPoint{{ cp[0], cp[p], cp[0] }}); // (0,1,0)
          s_nodes.emplace_back(Math::SpatialPoint{{ cp[0], cp[0], cp[p] }}); // (0,0,1)

          // Edges, in the same order as MFEM comments (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)

          // (0,1): (cp[i], cp[0], cp[0])
          for (int i = 1; i < p; ++i)
            s_nodes.emplace_back(Math::SpatialPoint{{ cp[i], cp[0], cp[0] }});

          // (0,2): (cp[0], cp[i], cp[0])
          for (int i = 1; i < p; ++i)
            s_nodes.emplace_back(Math::SpatialPoint{{ cp[0], cp[i], cp[0] }});

          // (0,3): (cp[0], cp[0], cp[i])
          for (int i = 1; i < p; ++i)
            s_nodes.emplace_back(Math::SpatialPoint{{ cp[0], cp[0], cp[i] }});

          // (1,2): (cp[p-i], cp[i], cp[0])
          for (int i = 1; i < p; ++i)
            s_nodes.emplace_back(Math::SpatialPoint{{ cp[p - i], cp[i], cp[0] }});

          // (1,3): (cp[p-i], cp[0], cp[i])
          for (int i = 1; i < p; ++i)
            s_nodes.emplace_back(Math::SpatialPoint{{ cp[p - i], cp[0], cp[i] }});

          // (2,3): (cp[0], cp[p-i], cp[i])
          for (int i = 1; i < p; ++i)
            s_nodes.emplace_back(Math::SpatialPoint{{ cp[0], cp[p - i], cp[i] }});

          // Faces
          // Face (1,2,3): barycentric (λ0, λ1, λ2, λ3) with λ0 = 0
          for (int j = 1; j < p; ++j)
          {
            for (int i = 1; i + j < p; ++i)
            {
              const Real ci = cp[i];
              const Real cj = cp[j];
              const Real ck = cp[p - i - j];
              const Real w  = ci + cj + ck;

              const Real x = cp[p - i - j] / w; // λ1
              const Real y = cp[i]             / w; // λ2
              const Real z = cp[j]             / w; // λ3

              s_nodes.emplace_back(Math::SpatialPoint{{ x, y, z }});
            }
          }

          // Face (0,3,2): λ1 = 0
          for (int j = 1; j < p; ++j)
          {
            for (int i = 1; i + j < p; ++i)
            {
              const Real ci = cp[i];
              const Real cj = cp[j];
              const Real ck = cp[p - i - j];
              const Real w  = ci + cj + ck;

              const Real x = cp[0];            // λ0
              const Real y = cp[j] / w;        // λ2
              const Real z = cp[i] / w;        // λ3

              s_nodes.emplace_back(Math::SpatialPoint{{ x, y, z }});
            }
          }

          // Face (0,1,3): λ2 = 0
          for (int j = 1; j < p; ++j)
          {
            for (int i = 1; i + j < p; ++i)
            {
              const Real ci = cp[i];
              const Real cj = cp[j];
              const Real ck = cp[p - i - j];
              const Real w  = ci + cj + ck;

              const Real x = cp[i] / w;        // λ0
              const Real y = cp[0];            // λ1
              const Real z = cp[j] / w;        // λ3

              s_nodes.emplace_back(Math::SpatialPoint{{ x, y, z }});
            }
          }

          // Face (0,2,1): λ3 = 0
          for (int j = 1; j < p; ++j)
          {
            for (int i = 1; i + j < p; ++i)
            {
              const Real ci = cp[i];
              const Real cj = cp[j];
              const Real ck = cp[p - i - j];
              const Real w  = ci + cj + ck;

              const Real x = cp[j] / w;        // λ0
              const Real y = cp[i] / w;        // λ1
              const Real z = cp[0];            // λ2

              s_nodes.emplace_back(Math::SpatialPoint{{ x, y, z }});
            }
          }

          // Interior
          for (int k = 1; k < p; ++k)
          {
            for (int j = 1; j + k < p; ++j)
            {
              for (int i = 1; i + j + k < p; ++i)
              {
                const Real ci = cp[i];
                const Real cj = cp[j];
                const Real ck = cp[k];
                const Real cl = cp[p - i - j - k];
                const Real w  = ci + cj + ck + cl;

                const Real x = ci / w;
                const Real y = cj / w;
                const Real z = ck / w;

                s_nodes.emplace_back(Math::SpatialPoint{{ x, y, z }});
              }
            }
          }

          assert(s_nodes.size() == Count);
        }

        return s_nodes;
      }
  };

  /**
   * @brief Vandermonde matrix on the tetrahedron using MFEM's nodal set.
   *
   * V_MFEM(i,j) = ψ_j(x_i), where ψ_j are Dubiner tetrahedral modes and
   * x_i are MfemTetrahedronNodes<K>.
   */
  template <size_t K>
  class VandermondeTetrahedron
  {
    public:
      static const Math::Matrix<Real>& getMatrix()
      {
        static thread_local Math::Matrix<Real> s_vandermonde;

        constexpr size_t N = TetrahedronNodes<K>::Count;

        if (s_vandermonde.size() == 0)
        {
          const auto& nodes = TetrahedronNodes<K>::getNodes();
          s_vandermonde.resize(N, N);

          size_t mode_idx = 0;
          Rodin::Utility::ForIndex<K + 1>(
            [&](auto p_idx)
            {
              constexpr size_t P = p_idx.value;
              Rodin::Utility::ForIndex<K + 1 - P>(
                [&](auto q_idx)
                {
                  constexpr size_t Q = q_idx.value;
                  Rodin::Utility::ForIndex<K + 1 - P - Q>(
                    [&](auto r_idx)
                    {
                      constexpr size_t R = r_idx.value;
                      for (size_t node_idx = 0; node_idx < N; ++node_idx)
                      {
                        const auto& pt = nodes[node_idx];
                        const Real x = pt.x();
                        const Real y = pt.y();
                        const Real z = pt.z();

                        Real a, b, c;
                        Variational::DubinerTetrahedron<K>::getCollapsed(a, b, c, x, y, z);

                        Variational::DubinerTetrahedron<K>::template getBasis<P, Q, R>(
                          s_vandermonde(node_idx, mode_idx), a, b, c);
                      }
                      ++mode_idx;
                    });
                });
            });
        }

        return s_vandermonde;
      }

      static const Math::Matrix<Real>& getInverse()
      {
        static thread_local Math::Matrix<Real> s_inv;

        if (s_inv.size() == 0)
        {
          const auto& V = VandermondeTetrahedron<K>::getMatrix();
          Eigen::BDCSVD<Math::Matrix<Real>> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
          const Math::Matrix<Real> I = Math::Matrix<Real>::Identity(V.rows(), V.cols());
          s_inv = svd.solve(I);
        }

        return s_inv;
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
      using FESType    = Variational::H1<K, Range, Geometry::Mesh<Context::Local>>;
      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
      using DataType   = Math::Vector<ScalarType>;
      using ObjectType = Variational::GridFunction<FESType, DataType>;
      using Parent     = GridFunctionLoaderBase<FESType, DataType>;

      GridFunctionLoader(ObjectType& gf)
        : Parent(gf)
      {}

      void load(std::istream& is) override
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::char_;

        // -------------------------------------------------------------
        // 1. Parse MFEM GridFunction header
        // -------------------------------------------------------------
        MFEM::GridFunctionHeader header;

        const auto get_fec      = [&](auto& ctx) { header.fec      = _attr(ctx); };
        const auto get_vdim     = [&](auto& ctx) { header.vdim     = _attr(ctx); };
        const auto get_ordering = [&](auto& ctx)
        {
          header.ordering = static_cast<MFEM::Ordering>(_attr(ctx));
        };

        std::string line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);

        auto it = line.begin();
        const auto pfes  = boost::spirit::x3::string("FiniteElementSpace");
        const bool rfes  = boost::spirit::x3::phrase_parse(it, line.end(), pfes, space);
        assert(rfes && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pfec  =
          boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[get_fec];
        const bool rfec  = boost::spirit::x3::phrase_parse(it, line.end(), pfec, space);
        assert(rfec && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pvdim = boost::spirit::x3::string("VDim:") >> uint_[get_vdim];
        const bool rvdim = boost::spirit::x3::phrase_parse(it, line.end(), pvdim, space);
        assert(rvdim && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pord =
          boost::spirit::x3::string("Ordering:") >> uint_[get_ordering];
        const bool rord = boost::spirit::x3::phrase_parse(it, line.end(), pord, space);
        assert(rord && it == line.end());

        // -------------------------------------------------------------
        // 2. Basic objects and sizes
        // -------------------------------------------------------------
        auto& gf        = this->getObject();
        auto& data      = gf.getData();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        const size_t vdim       = fes.getVectorDimension();
        const size_t D          = mesh.getDimension();
        const size_t scalarSize = fes.getSize() / vdim;

        assert(header.vdim == vdim);
        if (data.size() == 0)
          return;

        assert(static_cast<size_t>(data.size()) == vdim * scalarSize);

        // -------------------------------------------------------------
        // 3. Precompute triangle / tetra change-of-nodes blocks
        // -------------------------------------------------------------
        constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
        constexpr size_t TetN = MFEM::TetrahedronNodes<K>::Count;

        struct TriChangeBlocks
        {
          Math::Matrix<ScalarType> C_ii_inv; // (nTriInt x nTriInt)
          Math::Matrix<ScalarType> C_ib;     // (nTriInt x nBnd)
          std::vector<size_t>      col_bnd;  // Fekete boundary columns
          std::vector<size_t>      col_int;  // Fekete interior columns
          std::vector<size_t>      row_int;  // MFEM interior rows
          size_t nInt() const { return col_int.size(); }
        };

        struct TetChangeBlocks
        {
          Math::Matrix<ScalarType> C_ii_inv; // (nTetInt x nTetInt)
          Math::Matrix<ScalarType> C_ib;     // (nTetInt x nBnd)
          std::vector<size_t>      col_bnd;  // Fekete boundary columns
          std::vector<size_t>      col_int;  // Fekete interior columns
          std::vector<size_t>      row_int;  // MFEM interior rows
          size_t nInt() const { return col_int.size(); }
        };

        static thread_local TriChangeBlocks triBlocks;
        static thread_local TetChangeBlocks tetBlocks;

        const int p_int   = static_cast<int>(K);
        const int nTriInt = (p_int > 1) ? (p_int - 1) * (p_int - 2) / 2 : 0;
        const int nTetInt = (p_int > 2) ? (p_int - 1) * (p_int - 2) * (p_int - 3) / 6 : 0;

        if (nTriInt > 0 && triBlocks.C_ii_inv.size() == 0)
        {
          const auto& V_mfem   = MFEM::VandermondeTriangle<K>::getMatrix();   // TriN x TriN
          const auto& V_rodInv = Variational::VandermondeTriangle<K>::getInverse(); // TriN x TriN
          const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
          Math::Matrix<ScalarType> C = C_real.template cast<ScalarType>();

          const int nV = 3;
          const int nE = 3 * (p_int - 1);
          const int triInteriorOffset = nV + nE; // 3 + 3(p-1) = 3p

          // MFEM interior rows: those corresponding to MFEM interior triangle nodes
          triBlocks.row_int.clear();
          for (int r = triInteriorOffset; r < static_cast<int>(TriN); ++r)
            triBlocks.row_int.push_back(static_cast<size_t>(r));

          // Fekete interior vs boundary columns (based on barycentric coords)
          triBlocks.col_bnd.clear();
          triBlocks.col_int.clear();

          const auto& fk = Variational::FeketeTriangle<K>::getNodes();
          const Real tol = 1e-12;

          for (size_t j = 0; j < TriN; ++j)
          {
            assert(j < fk.size());
            const auto& pt = fk[j];
            const Real x = pt.x();
            const Real y = pt.y();
            const Real l0 = Real(1) - x - y;
            const Real l1 = x;
            const Real l2 = y;
            if (l0 > tol && l1 > tol && l2 > tol)
              triBlocks.col_int.push_back(j);
            else
              triBlocks.col_bnd.push_back(j);
          }

          assert(static_cast<int>(triBlocks.col_int.size()) == nTriInt);

          const int nInt = nTriInt;
          const int nBnd = static_cast<int>(triBlocks.col_bnd.size());

          triBlocks.C_ii_inv.resize(nInt, nInt);
          triBlocks.C_ib.resize(nInt, nBnd);

          for (int i = 0; i < nInt; ++i)
          {
            assert(i < static_cast<int>(triBlocks.row_int.size()));
            const size_t ri = triBlocks.row_int[static_cast<size_t>(i)];
            for (int k = 0; k < nInt; ++k)
            {
              assert(k < static_cast<int>(triBlocks.col_int.size()));
              const size_t cj = triBlocks.col_int[static_cast<size_t>(k)];
              triBlocks.C_ii_inv(i, k) = C(static_cast<Eigen::Index>(ri),
                                           static_cast<Eigen::Index>(cj));
            }
            for (int k = 0; k < nBnd; ++k)
            {
              assert(k < static_cast<int>(triBlocks.col_bnd.size()));
              const size_t cj = triBlocks.col_bnd[static_cast<size_t>(k)];
              triBlocks.C_ib(i, k) = C(static_cast<Eigen::Index>(ri),
                                       static_cast<Eigen::Index>(cj));
            }
          }

          Eigen::BDCSVD<Math::Matrix<ScalarType>> svd(
            triBlocks.C_ii_inv, Eigen::ComputeThinU | Eigen::ComputeThinV);
          Math::Matrix<ScalarType> I =
            Math::Matrix<ScalarType>::Identity(nInt, nInt);
          triBlocks.C_ii_inv = svd.solve(I);
        }

        if (nTetInt > 0 && tetBlocks.C_ii_inv.size() == 0)
        {
          const auto& V_mfem   = MFEM::VandermondeTetrahedron<K>::getMatrix();   // TetN x TetN
          const auto& V_rodInv = Variational::VandermondeTetrahedron<K>::getInverse(); // TetN x TetN
          const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
          Math::Matrix<ScalarType> C = C_real.template cast<ScalarType>();

          const int nV  = 4;
          const int nE  = 6 * (p_int - 1);
          const int nF  = 4 * (p_int - 1) * (p_int - 2) / 2;
          const int tetInteriorOffset = nV + nE + nF;

          tetBlocks.row_int.clear();
          for (int r = tetInteriorOffset; r < static_cast<int>(TetN); ++r)
            tetBlocks.row_int.push_back(static_cast<size_t>(r));

          tetBlocks.col_bnd.clear();
          tetBlocks.col_int.clear();

          const auto& fk = Variational::FeketeTetrahedron<K>::getNodes();
          const Real tol = 1e-12;
          for (size_t j = 0; j < TetN; ++j)
          {
            const auto& pt = fk[j];
            const Real x = pt.x();
            const Real y = pt.y();
            const Real z = pt.z();
            const Real l0 = Real(1) - x - y - z;
            const Real l1 = x;
            const Real l2 = y;
            const Real l3 = z;
            if (l0 > tol && l1 > tol && l2 > tol && l3 > tol)
              tetBlocks.col_int.push_back(j);
            else
              tetBlocks.col_bnd.push_back(j);
          }

          assert(static_cast<int>(tetBlocks.col_int.size()) == nTetInt);

          const int nInt = nTetInt;
          const int nBnd = static_cast<int>(tetBlocks.col_bnd.size());

          tetBlocks.C_ii_inv.resize(nInt, nInt);
          tetBlocks.C_ib.resize(nInt, nBnd);

          for (int i = 0; i < nInt; ++i)
          {
            assert(i < static_cast<int>(tetBlocks.row_int.size()));
            const size_t ri = tetBlocks.row_int[static_cast<size_t>(i)];
            for (int k = 0; k < nInt; ++k)
            {
              assert(k < static_cast<int>(tetBlocks.col_int.size()));
              const size_t cj = tetBlocks.col_int[static_cast<size_t>(k)];
              tetBlocks.C_ii_inv(i, k) = C(static_cast<Eigen::Index>(ri),
                                           static_cast<Eigen::Index>(cj));
            }
            for (int k = 0; k < nBnd; ++k)
            {
              assert(k < static_cast<int>(tetBlocks.col_bnd.size()));
              const size_t cj = tetBlocks.col_bnd[static_cast<size_t>(k)];
              tetBlocks.C_ib(i, k) = C(static_cast<Eigen::Index>(ri),
                                       static_cast<Eigen::Index>(cj));
            }
          }

          Eigen::BDCSVD<Math::Matrix<ScalarType>> svd(
            tetBlocks.C_ii_inv, Eigen::ComputeThinU | Eigen::ComputeThinV);
          Math::Matrix<ScalarType> I =
            Math::Matrix<ScalarType>::Identity(nInt, nInt);
          tetBlocks.C_ii_inv = svd.solve(I);
        }

        // -------------------------------------------------------------
        // 4. Build global MFEM scalar index → Rodin scalar DOF map.
        //    We mirror *exactly* the traversal in the printer:
        //
        //      1. vertices
        //      2. edges (interior, oriented)
        //      3. faces:
        //         - D==3: triangle faces → change-of-nodes (interior only)
        //                 quads → DOF-based
        //         - D==2: faces (edges) → DOF-based
        //      4. cell interiors:
        //         - D==2: triangles → change-of-nodes (interior only)
        //                 quads → DOF-based
        //         - D==3: tets → change-of-nodes (interior only)
        //                 wedges → change-of-nodes (triangle slices)
        //                 others → DOF-based
        // -------------------------------------------------------------
        const size_t nVertices = mesh.getConnectivity().getCount(0);
        const size_t nEdges    = mesh.getConnectivity().getCount(1);
        const size_t nFaces    = (D >= 2) ? mesh.getConnectivity().getCount(2) : 0;
        const size_t nCells    = mesh.getConnectivity().getCount(D);

        std::vector<std::optional<size_t>> dof2pos(scalarSize, std::nullopt);
        std::vector<uint8_t> seen(scalarSize, uint8_t(0));

        // Interior blocks for change-of-nodes entities
        std::vector<std::optional<size_t>> triFaceStart;
        std::vector<std::optional<size_t>> triCellStart;
        std::vector<std::optional<size_t>> tetCellStart;
        std::vector<std::optional<size_t>> wedgeCellStart;

        if (D == 3 && nFaces > 0)
          triFaceStart.assign(nFaces, std::nullopt);
        if (D == 2 && nCells > 0)
          triCellStart.assign(nCells, std::nullopt);
        if (D == 3 && nCells > 0)
        {
          tetCellStart.assign(nCells, std::nullopt);
          wedgeCellStart.assign(nCells, std::nullopt);
        }

        // Helper to convert potentially global DOF index to scalar DOF index
        // For vector-valued spaces, getDOFs() may return indices in [0, vdim*scalarSize)
        // We need to extract the scalar DOF in [0, scalarSize)
        auto toScalarDOF = [&](Index dof) -> Index
        {
          // If vdim > 1 and DOF seems to be global, extract scalar part
          if (vdim > 1 && dof >= scalarSize)
          {
            // Assuming interleaved or block ordering, try modulo
            return dof % scalarSize;
          }
          return dof;
        };

        // v → scalar DOF
        std::vector<Index> vertexScalarDof(nVertices);
        for (size_t v = 0; v < nVertices; ++v)
        {
          const auto& vdofs = fes.getDOFs(0, static_cast<Index>(v));
          assert(vdofs.size() >= 1);
          vertexScalarDof[v] = toScalarDOF(vdofs(0));
        }

        size_t pos = 0; // MFEM scalar index (0..scalarSize-1)

        auto visit_dof = [&](Index d)
        {
          const Index scalar_dof = toScalarDOF(d);
          const size_t s = static_cast<size_t>(scalar_dof);
          if (s >= scalarSize)
          {
            // Skip invalid DOF indices that are out of range
            return;
          }
          if (!seen[s])
          {
            seen[s]   = uint8_t(1);
            assert(s < dof2pos.size());
            dof2pos[s] = pos++;
          }
          // if already seen, printer's emit_scalar_dof would skip -> no new MFEM value
        };

        const int nWedgeInt = (p_int > 1 && nTriInt > 0) ? (p_int - 1) * nTriInt : 0;

        // 4.1 Vertices
        for (size_t v = 0; v < nVertices; ++v)
        {
          assert(v < vertexScalarDof.size());
          visit_dof(vertexScalarDof[v]);
        }

        // 4.2 Edges (interior DOFs, oriented)
        if (D >= 1)
        {
          const auto& conn10 = mesh.getConnectivity().getIncidence(1, 0);

          for (size_t e = 0; e < nEdges; ++e)
          {
            const auto& edgeVerts = conn10[static_cast<Index>(e)];
            assert(edgeVerts.size() == 2);

            Index v0 = edgeVerts[0];
            Index v1 = edgeVerts[1];

            const Index vmin    = std::min(v0, v1);
            const Index vmax    = std::max(v0, v1);
            assert(static_cast<size_t>(vmin) < vertexScalarDof.size());
            const Index vminDof = vertexScalarDof[static_cast<size_t>(vmin)];
            assert(static_cast<size_t>(vmax) < vertexScalarDof.size());
            const Index vmaxDof = vertexScalarDof[static_cast<size_t>(vmax)];

            const auto& edofs = fes.getDOFs(1, static_cast<Index>(e));

            std::vector<Index> interior;
            interior.reserve(edofs.size());

            for (size_t k = 0; k < edofs.size(); ++k)
            {
              const Index d = edofs(static_cast<Index>(k));
              if (d != vminDof && d != vmaxDof)
                interior.push_back(d);
            }

            if (v0 > v1)
              std::reverse(interior.begin(), interior.end());

            for (Index d : interior)
              visit_dof(d);
          }
        }

        // 4.3 Faces
        if (D >= 2)
        {
          const size_t faceDim = (D == 3) ? 2 : (D - 1);

          if (D == 3)
          {
            for (size_t f = 0; f < nFaces; ++f)
            {
              const auto faceGeom = mesh.getGeometry(2, static_cast<Index>(f));
              const auto& fdofs   = fes.getDOFs(faceDim, static_cast<Index>(f));

              if (faceGeom == Geometry::Polytope::Type::Triangle && nTriInt > 0)
              {
                assert(f < triFaceStart.size());
                triFaceStart[f] = pos;
                pos += static_cast<size_t>(nTriInt);
              }
              else
              {
                // Quads / others: DOF-based
                for (size_t k = 0; k < fdofs.size(); ++k)
                  visit_dof(fdofs(static_cast<Index>(k)));
              }
            }
          }
          else // D == 2
          {
            // Faces are edges; printer re-emits them, but emit_scalar_dof suppresses duplicates.
            for (size_t f = 0; f < nFaces; ++f)
            {
              const auto& fdofs = fes.getDOFs(faceDim, static_cast<Index>(f));
              for (size_t k = 0; k < fdofs.size(); ++k)
                visit_dof(fdofs(static_cast<Index>(k)));
            }
          }
        }

        // 4.4 Cells
        if (D == 2)
        {
          for (size_t c = 0; c < nCells; ++c)
          {
            const auto geom   = mesh.getGeometry(2, static_cast<Index>(c));
            const auto& cdofs = fes.getDOFs(2, static_cast<Index>(c));

            if (geom == Geometry::Polytope::Type::Triangle && nTriInt > 0)
            {
              assert(c < triCellStart.size());
              triCellStart[c] = pos;
              pos += static_cast<size_t>(nTriInt);
            }
            else
            {
              for (size_t k = 0; k < cdofs.size(); ++k)
                visit_dof(cdofs(static_cast<Index>(k)));
            }
          }
        }
        else if (D == 3)
        {
          for (size_t c = 0; c < nCells; ++c)
          {
            const auto geom   = mesh.getGeometry(3, static_cast<Index>(c));
            const auto& cdofs = fes.getDOFs(3, static_cast<Index>(c));

            if (geom == Geometry::Polytope::Type::Tetrahedron && nTetInt > 0)
            {
              assert(c < tetCellStart.size());
              tetCellStart[c] = pos;
              pos += static_cast<size_t>(nTetInt);
            }
            else if (geom == Geometry::Polytope::Type::Wedge && nWedgeInt > 0)
            {
              assert(c < wedgeCellStart.size());
              wedgeCellStart[c] = pos;
              pos += static_cast<size_t>(nWedgeInt);
            }
            else
            {
              // Hexahedra / others: DOF-based
              for (size_t k = 0; k < cdofs.size(); ++k)
                visit_dof(cdofs(static_cast<Index>(k)));
            }
          }
        }
        else
        {
          // D == 1: cells are segments / points → DOF-based
          for (size_t c = 0; c < nCells; ++c)
          {
            const auto& cdofs = fes.getDOFs(D, static_cast<Index>(c));
            for (size_t k = 0; k < cdofs.size(); ++k)
              visit_dof(cdofs(static_cast<Index>(k)));
          }
        }

        assert(pos == scalarSize);

        // -------------------------------------------------------------
        // 5. Read MFEM data (all components, all scalar positions)
        // -------------------------------------------------------------
        std::vector<ScalarType> mfem_data(static_cast<size_t>(data.size()));

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        assert(mfem_data.size() >= 1);
        mfem_data[0] = static_cast<ScalarType>(std::stod(line));

        for (size_t i = 1; i < mfem_data.size(); ++i)
        {
          ScalarType val;
          is >> val;
          mfem_data[i] = val;
        }

        auto mfem_value = [&](size_t p_s, size_t c) -> ScalarType
        {
          size_t idx = 0;
          if (header.ordering == MFEM::Ordering::Nodes)
            idx = p_s * vdim + c;
          else
            idx = c * scalarSize + p_s;
          assert(idx < mfem_data.size());
          return mfem_data[idx];
        };

        // -------------------------------------------------------------
        // 6. Fill DOFs written directly (no change-of-nodes)
        // -------------------------------------------------------------
        // DOFs with dof2pos[d].has_value() were printed via emit_scalar_dof.
        for (size_t s = 0; s < scalarSize; ++s)
        {
          assert(s < dof2pos.size());
          const auto& p_s_opt = dof2pos[s];
          if (!p_s_opt.has_value())
            continue; // will be handled via change-of-nodes
          const size_t p_s = p_s_opt.value();
          for (size_t c = 0; c < vdim; ++c)
          {
            assert(static_cast<Index>(s + c * scalarSize) < data.size());
            data.coeffRef(static_cast<Index>(s + c * scalarSize)) = mfem_value(p_s, c);
          }
        }

        // -------------------------------------------------------------
        // 7. Undo change-of-nodes for interiors
        // -------------------------------------------------------------

        // Helper for triangle patches (faces or cells)
        auto invert_triangle_patch = [&](const std::vector<Index>& local_dofs,
                                         size_t mfem_start, size_t comp)
        {
          const auto& blk = triBlocks;
          const size_t nIntLoc = blk.col_int.size();
          if (nIntLoc == 0)
            return;

          // Boundary Fekete DOFs
          Math::Vector<ScalarType> u_bnd(blk.col_bnd.size());
          for (size_t k = 0; k < blk.col_bnd.size(); ++k)
          {
            const size_t loc = blk.col_bnd[k];
            assert(loc < local_dofs.size());
            const Index d = local_dofs[loc];
            assert(static_cast<Index>(k) < u_bnd.size());
            assert(static_cast<Index>(d + comp * scalarSize) < data.size());
            u_bnd(static_cast<Index>(k)) =
              data.coeffRef(d + static_cast<Index>(comp * scalarSize));
          }

          // Interior MFEM values
          Math::Vector<ScalarType> y_int(blk.row_int.size());
          for (size_t i = 0; i < blk.row_int.size(); ++i)
          {
            const size_t p_s = mfem_start + i;
            assert(static_cast<Index>(i) < y_int.size());
            y_int(static_cast<Index>(i)) = mfem_value(p_s, comp);
          }

          Math::Vector<ScalarType> rhs = y_int - blk.C_ib * u_bnd;
          Math::Vector<ScalarType> u_int = blk.C_ii_inv * rhs;

          for (size_t k = 0; k < blk.col_int.size(); ++k)
          {
            const size_t loc = blk.col_int[k];
            assert(loc < local_dofs.size());
            const Index d = local_dofs[loc];
            assert(static_cast<Index>(k) < u_int.size());
            assert(static_cast<Index>(d + comp * scalarSize) < data.size());
            data.coeffRef(d + static_cast<Index>(comp * scalarSize)) =
              u_int(static_cast<Index>(k));
          }
        };

        // 7.1 2D triangles (cell interiors)
        if (D == 2 && nTriInt > 0)
        {
          for (size_t c = 0; c < nCells; ++c)
          {
            if (mesh.getGeometry(2, static_cast<Index>(c))
                != Geometry::Polytope::Type::Triangle)
              continue;

            assert(c < triCellStart.size());
            const auto& start_opt = triCellStart[c];
            if (!start_opt.has_value())
              continue;

            const size_t start = start_opt.value();
            const auto& cdofs = fes.getDOFs(2, static_cast<Index>(c));
            assert(cdofs.size() == TriN);

            std::vector<Index> local(TriN);
            for (size_t i = 0; i < TriN; ++i)
            {
              assert(i < local.size());
              assert(static_cast<Index>(i) < cdofs.size());
              local[i] = cdofs(static_cast<Index>(i));
            }

            for (size_t comp = 0; comp < vdim; ++comp)
              invert_triangle_patch(local, start, comp);
          }
        }

        // 7.2 3D triangle faces (face interiors)
        if (D == 3 && nTriInt > 0 && !triFaceStart.empty())
        {
          for (size_t f = 0; f < nFaces; ++f)
          {
            if (mesh.getGeometry(2, static_cast<Index>(f))
                != Geometry::Polytope::Type::Triangle)
              continue;

            assert(f < triFaceStart.size());
            const auto& start_opt = triFaceStart[f];
            if (!start_opt.has_value())
              continue;

            const size_t start = start_opt.value();
            const auto& fdofs = fes.getDOFs(2, static_cast<Index>(f));
            assert(fdofs.size() == TriN);

            std::vector<Index> local(TriN);
            for (size_t i = 0; i < TriN; ++i)
            {
              assert(i < local.size());
              assert(static_cast<Index>(i) < fdofs.size());
              local[i] = fdofs(static_cast<Index>(i));
            }

            for (size_t comp = 0; comp < vdim; ++comp)
              invert_triangle_patch(local, start, comp);
          }
        }

        // 7.3 3D tetrahedra (cell interiors)
        if (D == 3 && nTetInt > 0)
        {
          const auto& blk = tetBlocks;

          for (size_t c = 0; c < nCells; ++c)
          {
            if (mesh.getGeometry(3, static_cast<Index>(c))
                != Geometry::Polytope::Type::Tetrahedron)
              continue;

            assert(c < tetCellStart.size());
            const auto& start_opt = tetCellStart[c];
            if (!start_opt.has_value())
              continue;

            const size_t start = start_opt.value();
            const auto& cdofs = fes.getDOFs(3, static_cast<Index>(c));
            assert(cdofs.size() == TetN);

            std::vector<Index> local(TetN);
            for (size_t i = 0; i < TetN; ++i)
            {
              assert(i < local.size());
              assert(static_cast<Index>(i) < cdofs.size());
              local[i] = cdofs(static_cast<Index>(i));
            }

            for (size_t comp = 0; comp < vdim; ++comp)
            {
              // Fekete boundary DOFs
              Math::Vector<ScalarType> u_bnd(blk.col_bnd.size());
              for (size_t k = 0; k < blk.col_bnd.size(); ++k)
              {
                const size_t loc = blk.col_bnd[k];
                const Index d = local[loc];
                assert(static_cast<Index>(k) < u_bnd.size());
                assert(static_cast<Index>(d + comp * scalarSize) < data.size());
                u_bnd(static_cast<Index>(k)) =
                  data.coeffRef(d + static_cast<Index>(comp * scalarSize));
              }

              // MFEM interior values
              Math::Vector<ScalarType> y_int(blk.row_int.size());
              for (size_t i = 0; i < blk.row_int.size(); ++i)
              {
                const size_t p_s = start + i;
                assert(static_cast<Index>(i) < y_int.size());
                y_int(static_cast<Index>(i)) = mfem_value(p_s, comp);
              }

              Math::Vector<ScalarType> rhs = y_int - blk.C_ib * u_bnd;
              Math::Vector<ScalarType> u_int = blk.C_ii_inv * rhs;

              for (size_t k = 0; k < blk.col_int.size(); ++k)
              {
                assert(k < blk.col_int.size());
                const size_t loc = blk.col_int[k];
                assert(loc < local.size());
                const Index d = local[loc];
                assert(static_cast<Index>(k) < u_int.size());
                assert(static_cast<Index>(d + comp * scalarSize) < data.size());
                data.coeffRef(d + static_cast<Index>(comp * scalarSize)) =
                  u_int(static_cast<Index>(k));
              }
            }
          }
        }

        // 7.4 3D wedges (cell interiors: tri-interior × segment-interior)
        if (D == 3 && nWedgeInt > 0 && nTriInt > 0)
        {
          for (size_t c = 0; c < nCells; ++c)
          {
            if (mesh.getGeometry(3, static_cast<Index>(c))
                != Geometry::Polytope::Type::Wedge)
              continue;

            assert(c < wedgeCellStart.size());
            const auto& start_opt = wedgeCellStart[c];
            if (!start_opt.has_value())
              continue;

            const size_t start = start_opt.value();
            const auto& cdofs = fes.getDOFs(3, static_cast<Index>(c));
            const size_t wedgeDofs = static_cast<size_t>((p_int + 1) * TriN);
            assert(cdofs.size() == wedgeDofs);

            // Local wedge DOFs: index = kseg * TriN + tri_idx, kseg = 0..p_int
            std::vector<Index> local(wedgeDofs);
            for (size_t i = 0; i < wedgeDofs; ++i)
            {
              assert(i < local.size());
              assert(static_cast<Index>(i) < cdofs.size());
              local[i] = cdofs(static_cast<Index>(i));
            }

            for (size_t comp = 0; comp < vdim; ++comp)
            {
              for (int kseg = 1; kseg < p_int; ++kseg) // segment interior: 1..p-1
              {
                // Triangle slice at fixed kseg
                std::vector<Index> local_tri(TriN);
                for (size_t it = 0; it < TriN; ++it)
                {
                  const size_t loc = static_cast<size_t>(kseg) * TriN + it;
                  assert(loc < local.size());
                  local_tri[it] = local[loc];
                }

                const size_t slice_start =
                  start + static_cast<size_t>(kseg - 1) * static_cast<size_t>(nTriInt);

                invert_triangle_patch(local_tri, slice_start, comp);
              }
            }
          }
        }
      }

    private:
      size_t m_currentLineNumber = 0;
  };

  /**
   * @brief Loader for P0 (piecewise constant) grid functions in MFEM format.
   *
   * Loads discontinuous finite element solutions from MFEM format.
   * P0 spaces have one DOF per element, stored in element order.
   *
   * @tparam Range Range type for the finite element space
   *
   * ## Expected MFEM Format
   * - Header: FiniteElementSpace
   * - FiniteElementCollection: L2_<dim>D_P0
   * - VDim: vector dimension
   * - Ordering: data layout (Nodes=0 or VectorDimension=1)
   * - Coefficient data values (one per element, in element order)
   *
   * ## Usage Example
   * ```cpp
   * P0 fes(mesh);
   * GridFunction gf(fes);
   * GridFunctionLoader<FileFormat::MFEM, P0<Real>, Vector<Real>> loader(gf);
   * loader.load("solution.gf");
   * ```
   *
   * @see GridFunctionPrinter
   */
  template <class Range>
  class GridFunctionLoader<
    FileFormat::MFEM,
    Variational::P0<Range, Geometry::Mesh<Context::Local>>,
    Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
    : public GridFunctionLoaderBase<
        Variational::P0<Range, Geometry::Mesh<Context::Local>>,
        Math::Vector<typename FormLanguage::Traits<Range>::ScalarType>>
  {
    public:
      using FESType = Variational::P0<Range, Geometry::Mesh<Context::Local>>;

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
       * Parses the header and coefficient data. P0 elements are discontinuous,
       * so data is simply ordered by element index.
       */
      void load(std::istream& is) override
      {
        using boost::spirit::x3::space;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::_attr;
        using boost::spirit::x3::char_;

        // -------------------------------------------------------------
        // 1. Parse MFEM GridFunction header
        // -------------------------------------------------------------
        MFEM::GridFunctionHeader header;

        const auto get_fec      = [&](auto& ctx) { header.fec      = _attr(ctx); };
        const auto get_vdim     = [&](auto& ctx) { header.vdim     = _attr(ctx); };
        const auto get_ordering = [&](auto& ctx)
        {
          header.ordering = static_cast<MFEM::Ordering>(_attr(ctx));
        };

        std::string line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);

        auto it = line.begin();
        const auto pfes  = boost::spirit::x3::string("FiniteElementSpace");
        const bool rfes  = boost::spirit::x3::phrase_parse(it, line.end(), pfes, space);
        assert(rfes && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pfec  =
          boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[get_fec];
        const bool rfec  = boost::spirit::x3::phrase_parse(it, line.end(), pfec, space);
        assert(rfec && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pvdim = boost::spirit::x3::string("VDim:") >> uint_[get_vdim];
        const bool rvdim = boost::spirit::x3::phrase_parse(it, line.end(), pvdim, space);
        assert(rvdim && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pord =
          boost::spirit::x3::string("Ordering:") >> uint_[get_ordering];
        const bool rord = boost::spirit::x3::phrase_parse(it, line.end(), pord, space);
        assert(rord && it == line.end());

        // -------------------------------------------------------------
        // 2. Read coefficient data
        // -------------------------------------------------------------
        auto& gf  = this->getObject();
        auto& fes = gf.getFiniteElementSpace();
        auto& data = gf.getData();

        const size_t vdim = fes.getVectorDimension();
        const size_t dofCount = fes.getSize();

        // Resize data vector
        data.resize(dofCount);

        // P0 spaces are discontinuous: one DOF per element
        // MFEM orders by elements, and Rodin also orders by elements for P0
        // So we can read directly in order
        for (Index i = 0; i < static_cast<Index>(dofCount); i++)
        {
          line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
          it = line.begin();

          ScalarType value;
          const auto get_value = [&](auto& ctx) { value = _attr(ctx); };

          using boost::spirit::x3::double_;
          const auto pvalue = double_[get_value];
          const bool rvalue = boost::spirit::x3::phrase_parse(it, line.end(), pvalue, space);
          assert(rvalue && it == line.end());

          data(i) = value;
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
      using FESType    = Variational::H1<K, Range, Geometry::Mesh<Context::Local>>;
      using DataType   = Math::Vector<Scalar>;
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

        const std::streamsize old_prec  = os.precision();
        const std::ios::fmtflags old_fl = os.flags();

        os << std::setprecision(std::numeric_limits<Scalar>::digits10 + 2);
        os.setf(std::ios::scientific, std::ios::floatfield);

        // Track which "Rodin scalar DOFs" have already been emitted in the
        // fallback DOF-based paths (vertices / edges / non-tri/tet cells).
        std::vector<uint8_t> written(scalarSize, uint8_t(0));

        // Helper to convert potentially global DOF index to scalar DOF index
        // For vector-valued spaces, getDOFs() may return indices in [0, vdim*scalarSize)
        // We need to extract the scalar DOF in [0, scalarSize)
        auto toScalarDOF = [&](Index dof) -> Index
        {
          // If vdim > 1 and DOF seems to be global, extract scalar part
          if (vdim > 1 && dof >= scalarSize)
          {
            // Assuming interleaved or block ordering, try modulo
            return dof % scalarSize;
          }
          return dof;
        };

        // Vertex -> scalar DOF (H1: exactly one scalar DOF per vertex)
        const size_t nVertices = mesh.getConnectivity().getCount(0);
        std::vector<Index> vertexScalarDof(nVertices);
        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
        {
          const auto& vdofs = fes.getDOFs(0, v);
          assert(vdofs.size() >= 1 && "H1 vertex should have at least one DOF.");
          vertexScalarDof[v] = toScalarDOF(vdofs(0));
        }

        // Emit coefficient(s) for a *Rodin scalar DOF index* in Nodes ordering
        auto emit_scalar_dof = [&](Index rodin_dof)
        {
          const Index scalar_dof = toScalarDOF(rodin_dof);
          const size_t s = static_cast<size_t>(scalar_dof);
          if (s >= scalarSize)
          {
            // Skip invalid DOF indices that are out of range
            return;
          }
          if (written[s])
            return;
          for (size_t c = 0; c < vdim; ++c)
            os << data.coeffRef(scalar_dof + c * scalarSize) << '\n';
          written[s] = uint8_t(1);
        };

        // Helper for range-valued gf(p) (only used in wedge / fallback)
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

        //--------------------------------------------------------------------------
        // Precomputed change-of-nodes matrices (Rodin Fekete -> MFEM nodes)
        //--------------------------------------------------------------------------

        auto& triChangeScalar = []() -> const Math::Matrix<Scalar>&
        {
          static thread_local Math::Matrix<Scalar> C;
          if (C.size() == 0)
          {
            const auto& V_mfem  = MFEM::VandermondeTriangle<K>::getMatrix();
            const auto& V_rodInv = Variational::VandermondeTriangle<K>::getInverse();
            const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
            C = C_real.template cast<Scalar>();
          }
          return C;
        }();

        auto& tetChangeScalar = []() -> const Math::Matrix<Scalar>&
        {
          static thread_local Math::Matrix<Scalar> C;
          if (C.size() == 0)
          {
            const auto& V_mfem  = MFEM::VandermondeTetrahedron<K>::getMatrix();
            const auto& V_rodInv = Variational::VandermondeTetrahedron<K>::getInverse();
            const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
            C = C_real.template cast<Scalar>();
          }
          return C;
        }();

        //--------------------------------------------------------------------------
        // 1. Vertices (shared between Rodin and MFEM)
        //--------------------------------------------------------------------------

        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
          emit_scalar_dof(vertexScalarDof[v]);

        //--------------------------------------------------------------------------
        // 2. Edges: interior DOFs, oriented low-vertex -> high-vertex
        //    (Rodin and MFEM use the same segment nodal set)
        //--------------------------------------------------------------------------

        if (D >= 1)
        {
          const auto& conn10  = mesh.getConnectivity().getIncidence(1, 0);
          const size_t nEdges = mesh.getConnectivity().getCount(1);

          std::vector<Index> interior;
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

            interior.clear();
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

        //--------------------------------------------------------------------------
        // 3. Faces
        //    - D == 2: faces are 1D edges -> already handled in edge block.
        //    - D == 3: 2D faces.
        //
        //      For triangular faces, we use the change-of-nodes operator
        //      on the face (Rodin Fekete -> MFEM triangle nodes).
        //      We emit only *interior* face nodes here; vertex and edge
        //      nodes were already emitted above.
        //--------------------------------------------------------------------------

        if (D >= 2)
        {
          const size_t faceDim   = (D == 3) ? 2 : (D - 1);
          const size_t faceCount = mesh.getConnectivity().getCount(faceDim);

          if (D == 3)
          {
            const auto& conn20 = mesh.getConnectivity().getIncidence(2, 0);

            // Local triangle size and interior offset
            constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
            const int p  = static_cast<int>(K);
            const int nV = 3;
            const int nE = 3 * (p - 1);
            const int triInteriorOffset = nV + nE; // first interior node index

            // Temporary storage for local transforms on faces
            Math::Vector<Scalar> uR_face(TriN);
            std::vector<Math::Vector<Scalar>> uM_face(vdim, Math::Vector<Scalar>(TriN));

            for (Index f = 0; f < static_cast<Index>(faceCount); ++f)
            {
              const auto faceGeom = mesh.getGeometry(2, f);

              if (faceGeom == Geometry::Polytope::Type::Triangle)
              {
                const auto& fdofs = fes.getDOFs(2, f);
                assert(static_cast<size_t>(fdofs.size()) == TriN
                       && "Triangle face must have (K+1)(K+2)/2 DOFs.");

                // Local change-of-nodes on the face, per component
                for (size_t c = 0; c < vdim; ++c)
                {
                  for (size_t k = 0; k < TriN; ++k)
                  {
                    const Index d = fdofs(static_cast<Index>(k));
                    uR_face(static_cast<Index>(k)) =
                      data.coeffRef(d + c * scalarSize);
                  }
                  uM_face[c] = triChangeScalar * uR_face;
                }

                // Emit MFEM face interior DOFs in MFEM's local order
                int loc = 0;
                for (int j = 1; j < p; ++j)
                {
                  for (int i = 1; i + j < p; ++i)
                  {
                    const int idx = triInteriorOffset + loc++;
                    for (size_t c = 0; c < vdim; ++c)
                      os << uM_face[c](static_cast<Index>(idx)) << '\n';
                  }
                }
              }
              else
              {
                // Non-triangular faces (e.g. quads): Rodin and MFEM share
                // the same nodal set; use DOF-based ordering as before.
                const auto& fdofs = fes.getDOFs(faceDim, f);
                for (Index k = 0; k < static_cast<Index>(fdofs.size()); ++k)
                  emit_scalar_dof(fdofs(k));
              }
            }
          }
          else
          {
            // D == 2 or other: faces are edges (already covered above) or
            // lower-dimensional; keep DOF-based behavior.
            for (Index f = 0; f < static_cast<Index>(faceCount); ++f)
            {
              const auto& fdofs = fes.getDOFs(faceDim, f);
              for (Index k = 0; k < static_cast<Index>(fdofs.size()); ++k)
                emit_scalar_dof(fdofs(k));
            }
          }
        }

        //--------------------------------------------------------------------------
        // 4. Element interiors
        //
        //    2D triangles and 3D tetrahedra use change-of-nodes operators.
        //    Other cells fall back to DOF-based ordering.
        //--------------------------------------------------------------------------

        const size_t nCells = mesh.getConnectivity().getCount(D);

        if (D == 2)
        {
          // 2D: triangles (possibly embedded in R^2 or R^3)
          constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
          const int p  = static_cast<int>(K);
          const int nV = 3;
          const int nE = 3 * (p - 1);
          const int triInteriorOffset = nV + nE;

          Math::Vector<Scalar> uR_elem(TriN);
          std::vector<Math::Vector<Scalar>> uM_elem(vdim, Math::Vector<Scalar>(TriN));

          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto geom = mesh.getGeometry(2, c);
            if (geom != Geometry::Polytope::Type::Triangle)
            {
              // Fallback: non-triangle cells in 2D (e.g. quads)
              const auto& cdofs = fes.getDOFs(D, c);
              for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
                emit_scalar_dof(cdofs(k));
              continue;
            }

            const auto& cdofs = fes.getDOFs(2, c);
            assert(static_cast<size_t>(cdofs.size()) == TriN
                   && "Triangle cell must have (K+1)(K+2)/2 DOFs.");

            // Local change-of-nodes on the triangle, per component
            for (size_t comp = 0; comp < vdim; ++comp)
            {
              for (size_t k = 0; k < TriN; ++k)
              {
                const Index d = cdofs(static_cast<Index>(k));
                uR_elem(static_cast<Index>(k)) =
                  data.coeffRef(d + comp * scalarSize);
              }
              uM_elem[comp] = triChangeScalar * uR_elem;
            }

            // Emit MFEM triangle interior DOFs in MFEM's local order
            int loc = 0;
            for (int j = 1; j < p; ++j)
            {
              for (int i = 1; i + j < p; ++i)
              {
                const int idx = triInteriorOffset + loc++;
                for (size_t comp = 0; comp < vdim; ++comp)
                  os << uM_elem[comp](static_cast<Index>(idx)) << '\n';
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

          // Tetra local sizes / offsets
          constexpr size_t TetN = MFEM::TetrahedronNodes<K>::Count;
          const int p = static_cast<int>(K);

          const int nV    = 4;
          const int nE    = 6 * (p - 1);
          const int nF    = 2 * (p - 1) * (p - 2); // 4 faces × (p-1)(p-2)/2
          const int tetInteriorOffset = nV + nE + nF;

          Math::Vector<Scalar> uR_elem(TetN);
          std::vector<Math::Vector<Scalar>> uM_elem(vdim, Math::Vector<Scalar>(TetN));

          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto  cellGeom = mesh.getGeometry(3, c);
            const auto& cellVer  = conn30[c];

            const auto  cellIt   = mesh.getCell(c);
            const auto& cell     = *cellIt;

            // ---------------------------------------------------------------
            // Tetrahedron: use change-of-nodes for element interior DOFs
            // ---------------------------------------------------------------
            if (cellGeom == Geometry::Polytope::Type::Tetrahedron)
            {
              assert(cellVer.size() == 4 && "Tetrahedron must have 4 vertices.");

              const auto& cdofs = fes.getDOFs(3, c);
              assert(static_cast<size_t>(cdofs.size()) == TetN
                     && "Tetrahedron must have (K+1)(K+2)(K+3)/6 DOFs.");

              // Local change-of-nodes, per component
              for (size_t comp = 0; comp < vdim; ++comp)
              {
                for (size_t k = 0; k < TetN; ++k)
                {
                  const Index d = cdofs(static_cast<Index>(k));
                  uR_elem(static_cast<Index>(k)) =
                    data.coeffRef(d + comp * scalarSize);
                }
                uM_elem[comp] = tetChangeScalar * uR_elem;
              }

              // Emit only element-interior nodes in MFEM's local order
              for (int idx = tetInteriorOffset; idx < static_cast<int>(TetN); ++idx)
              {
                for (size_t comp = 0; comp < vdim; ++comp)
                  os << uM_elem[comp](static_cast<Index>(idx)) << '\n';
              }
            }
            // ---------------------------------------------------------------
            // Wedge / prism: MFEM and Rodin share the same nodal set
            // on the wedge (segment × triangle), but we must be careful
            // with orientation. We keep the evaluation-based approach
            // you already had, which is exact and consistent.
            // ---------------------------------------------------------------
            else if (cellGeom == Geometry::Polytope::Type::Wedge)
            {
              // Rodin wedge local DOFs: product ordering
              //   index = k * TriN + tri_idx, k = 0..p, tri_idx = 0..TriN-1
              //
              // MFEM wedge interior DOFs: triangle-interior × segment-interior, with
              //   tri interior indices = 3*p + l, l = 0..nt-1
              //   seg interior indices = 1..p-1
              //
              // We map each triangle slice via triChangeScalar and then read
              // MFEM triangle interior entries for k = 1..p-1.

              constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
              const int p  = static_cast<int>(K);
              const int nV = 3;
              const int nE = 3 * (p - 1);
              const int triInteriorOffset = nV + nE;              // = 3*p
              const int nTriInt = (p - 1) * (p - 2) / 2;

              // Wedge has (p+1)*TriN DOFs per scalar component
              const auto& cdofs = fes.getDOFs(3, c);
              const size_t wedgeDofs = static_cast<size_t>((p + 1) * TriN);
              assert(static_cast<size_t>(cdofs.size()) == wedgeDofs &&
                     "Wedge element must have (p+1)*(p+1)*(p+2)/2 DOFs.");

              // uR_elem: Rodin local wedge DOFs for one component
              Math::Vector<Scalar> uR_elem(wedgeDofs);

              // uR_tri[k]: Rodin local triangle slice at segment index k
              std::vector<Math::Vector<Scalar>> uR_tri(p + 1, Math::Vector<Scalar>(TriN));

              // uM_tri[k]: MFEM local triangle slice at segment index k
              std::vector<Math::Vector<Scalar>> uM_tri(p + 1, Math::Vector<Scalar>(TriN));

              for (size_t comp = 0; comp < vdim; ++comp)
              {
                // 1) Collect Rodin wedge DOFs for this component
                for (size_t loc = 0; loc < wedgeDofs; ++loc)
                {
                  const Index d = cdofs(static_cast<Index>(loc));
                  uR_elem(static_cast<Index>(loc)) =
                    data.coeffRef(d + comp * scalarSize);
                }

                // 2) Split into triangle slices and apply triangle change-of-nodes
                for (int kseg = 0; kseg <= p; ++kseg)
                {
                  const size_t offset = static_cast<size_t>(kseg) * TriN;
                  for (size_t itri = 0; itri < TriN; ++itri)
                    uR_tri[kseg](static_cast<Index>(itri)) =
                      uR_elem(static_cast<Index>(offset + itri));

                  uM_tri[kseg] = triChangeScalar * uR_tri[kseg];
                }

                // 3) Emit wedge *interior* DOFs in MFEM's ordering:
                //    for k = 1..p-1 (segment interior) and triangle interior
                //    nodes in (j,i) loops (same as H1_WedgeElement).
                for (int kseg = 1; kseg < p; ++kseg)
                {
                  int l = 0; // interior triangle counter
                  for (int j = 1; j < p; ++j)
                  {
                    for (int i = 1; i + j < p; ++i)
                    {
                      const int triIdx = triInteriorOffset + l++;
                      os << uM_tri[kseg](static_cast<Index>(triIdx)) << '\n';
                    }
                  }
                  assert(l == nTriInt);
                }
              }
            }

            // ---------------------------------------------------------------
            // Other 3D cell types: fallback DOF-based ordering
            // ---------------------------------------------------------------
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
          // Other dimensions: fallback DOF-based ordering
          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto& cdofs = fes.getDOFs(D, c);
            for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
              emit_scalar_dof(cdofs(k));
          }
        }

        os.precision(old_prec);
        os.flags(old_fl);
      }
  };
}

#endif
