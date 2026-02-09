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
#include <iomanip>
#include <optional>
#include <limits>

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
      case GeometryType::CUBE:
      {
        return Rodin::Geometry::Polytope::Type::Hexahedron;
      }
      case GeometryType::PYRAMID:
      {
        return {};
      }
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
      case Geometry::Polytope::Type::Hexahedron:
        return GeometryType::CUBE;
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
        const auto get_fec      = [&](auto& ctx) { header.fec = _attr(ctx); };
        const auto get_vdim     = [&](auto& ctx) { header.vdim = _attr(ctx); };
        const auto get_ordering = [&](auto& ctx)
        {
          header.ordering = static_cast<MFEM::Ordering>(_attr(ctx));
        };

        std::string line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        auto it = line.begin();
        const auto pfes = boost::spirit::x3::string("FiniteElementSpace");
        const bool rfes = boost::spirit::x3::phrase_parse(it, line.end(), pfes, space);
        (void) rfes;
        assert(it == line.end() && rfes);

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it = line.begin();
        const auto pfec = boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[get_fec];
        const bool rfec = boost::spirit::x3::phrase_parse(it, line.end(), pfec, space);
        (void) rfec;
        assert(it == line.end() && rfec);

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it = line.begin();
        const auto pvdim = boost::spirit::x3::string("VDim:") >> uint_[get_vdim];
        const bool rvdim = boost::spirit::x3::phrase_parse(it, line.end(), pvdim, space);
        (void) rvdim;
        assert(it == line.end() && rvdim);

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it = line.begin();
        const auto pordering = boost::spirit::x3::string("Ordering:") >> uint_[get_ordering];
        const bool rordering = boost::spirit::x3::phrase_parse(it, line.end(), pordering, space);
        (void) rordering;
        assert(it == line.end() && rordering);

        auto& gf  = this->getObject();
        const auto& fes = gf.getFiniteElementSpace();

        const size_t vdim = fes.getVectorDimension();
        assert(header.vdim == vdim);

        auto& data = gf.getData();
        const size_t n = static_cast<size_t>(data.size());
        if (n == 0)
          return;

        // P1 scalar DOFs = vertices
        const size_t vn = fes.getMesh().getVertexCount();
        assert(n == vn * vdim && "P1 GridFunction size must be vertexCount * vdim");

        // Read all coefficients as they appear in the file
        std::vector<ScalarType> tmp(n);

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        tmp[0] = static_cast<ScalarType>(std::stod(line));
        for (size_t i = 1; i < n; ++i)
          is >> tmp[i];

        // Convert from MFEM file ordering -> Rodin internal (block-by-component):
        // internal index = v + c*vn
        if (header.ordering == MFEM::Ordering::Nodes) // Ordering: 0
        {
          for (size_t c = 0; c < vdim; ++c)
            for (size_t v = 0; v < vn; ++v)
              data.coeffRef(v + c * vn) = tmp[v + c * vn];
        }
        else // Ordering: 1 (VectorDimension)
        {
          for (size_t v = 0; v < vn; ++v)
            for (size_t c = 0; c < vdim; ++c)
              data.coeffRef(v + c * vn) = tmp[vdim * v + c];
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
        using boost::spirit::x3::blank;
        using boost::spirit::x3::uint_;
        using boost::spirit::x3::double_;
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
        (void) rfes;
        assert(rfes && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pfec =
          boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[get_fec];
        const bool rfec = boost::spirit::x3::phrase_parse(it, line.end(), pfec, space);
        (void) rfec;
        assert(rfec && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pvdim = boost::spirit::x3::string("VDim:") >> uint_[get_vdim];
        const bool rvdim = boost::spirit::x3::phrase_parse(it, line.end(), pvdim, space);
        (void) rvdim;
        assert(rvdim && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pord =
          boost::spirit::x3::string("Ordering:") >> uint_[get_ordering];
        const bool rord = boost::spirit::x3::phrase_parse(it, line.end(), pord, space);
        (void) rord;
        assert(rord && it == line.end());

        // -------------------------------------------------------------
        // 2. Read coefficient data in MFEM order and map to Rodin DOFs
        // -------------------------------------------------------------
        auto& gf   = this->getObject();
        const auto& fes  = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        auto& data = gf.getData();

        const size_t vdim       = fes.getVectorDimension();
        const size_t D          = mesh.getDimension();
        const size_t scalarSize = fes.getSize() / vdim;

        assert(header.vdim == vdim);

        data.resize(fes.getSize());

        // Read all values exactly as in the file
        std::vector<ScalarType> mfem_values;
        mfem_values.reserve(fes.getSize());

        while (true)
        {
          line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
          if (line.empty() || is.eof())
            break;

          it = line.begin();
          ScalarType value;
          const auto get_value = [&](auto& ctx) { value = _attr(ctx); };
          const auto pvalue = double_[get_value];
          const bool rvalue = boost::spirit::x3::phrase_parse(it, line.end(), pvalue, space);

          if (!rvalue || it != line.end())
            break;

          mfem_values.push_back(value);
        }

        assert(mfem_values.size() == fes.getSize() && "Mismatch in number of coefficients read");

        // -------------------------------------------------------------
        // 2.a Normalize MFEM ordering into component-blocked storage:
        //      mfem_block[c][i] = value of component c at MFEM scalar position i
        // -------------------------------------------------------------
        std::vector<std::vector<ScalarType>> mfem_block(vdim, std::vector<ScalarType>(scalarSize));

        if (header.ordering == MFEM::Ordering::Nodes) // 0: XXX..YYY..ZZZ..
        {
          for (size_t c = 0; c < vdim; ++c)
            for (size_t i = 0; i < scalarSize; ++i)
              mfem_block[c][i] = mfem_values[c * scalarSize + i];
        }
        else // 1: XYZ,XYZ,XYZ...
        {
          for (size_t i = 0; i < scalarSize; ++i)
            for (size_t c = 0; c < vdim; ++c)
              mfem_block[c][i] = mfem_values[i * vdim + c];
        }

        // Scalar-position cursor in the MFEM traversal (increments once per scalar DOF consumed)
        size_t pos = 0;

        // Helper: Rodin DOF -> scalar DOF index
        const auto to_scalar_dof = [&](Index dof) -> Index
        {
          if (vdim > 1 && dof >= static_cast<Index>(scalarSize))
            return dof % static_cast<Index>(scalarSize);
          return dof;
        };

        // Helper: assign (all components) for a Rodin scalar DOF from current MFEM scalar position
        auto set_scalar_dof_from_pos = [&](Index scalar_dof)
        {
          assert(pos < scalarSize);
          for (size_t c = 0; c < vdim; ++c)
            data.coeffRef(scalar_dof + static_cast<Index>(c * scalarSize)) = mfem_block[c][pos];
          ++pos;
        };

        // Track which Rodin scalar DOFs have been set
        std::vector<uint8_t> written(scalarSize, false);

        //--------------------------------------------------------------------
        // Precomputed change-of-nodes matrices (Rodin Fekete -> MFEM nodes)
        // and inverse matrices (MFEM nodes -> Rodin Fekete)
        //--------------------------------------------------------------------
        auto& s_tri_change_scalar = []() -> const Math::Matrix<ScalarType>&
        {
          static thread_local Math::Matrix<ScalarType> C;
          if (C.size() == 0)
          {
            const auto& V_mfem   = MFEM::VandermondeTriangle<K>::getMatrix();
            const auto& V_rodInv = Variational::VandermondeTriangle<K>::getInverse();
            const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
            C = C_real.template cast<ScalarType>();
          }
          return C;
        }();

        auto& s_tri_inv_change_scalar = []() -> const Math::Matrix<ScalarType>&
        {
          static thread_local Math::Matrix<ScalarType> Cinv;
          if (Cinv.size() == 0)
          {
            const auto& V_mfem   = MFEM::VandermondeTriangle<K>::getMatrix();
            const auto& V_rodInv = Variational::VandermondeTriangle<K>::getInverse();
            const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
            const Math::Matrix<Real> Cinv_real = C_real.inverse();
            Cinv = Cinv_real.template cast<ScalarType>();
          }
          return Cinv;
        }();

        auto& s_tet_change_scalar = []() -> const Math::Matrix<ScalarType>&
        {
          static thread_local Math::Matrix<ScalarType> C;
          if (C.size() == 0)
          {
            const auto& V_mfem   = MFEM::VandermondeTetrahedron<K>::getMatrix();
            const auto& V_rodInv = Variational::VandermondeTetrahedron<K>::getInverse();
            const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
            C = C_real.template cast<ScalarType>();
          }
          return C;
        }();

        auto& s_tet_inv_change_scalar = []() -> const Math::Matrix<ScalarType>&
        {
          static thread_local Math::Matrix<ScalarType> Cinv;
          if (Cinv.size() == 0)
          {
            const auto& V_mfem   = MFEM::VandermondeTetrahedron<K>::getMatrix();
            const auto& V_rodInv = Variational::VandermondeTetrahedron<K>::getInverse();
            const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
            const Math::Matrix<Real> Cinv_real = C_real.inverse();
            Cinv = Cinv_real.template cast<ScalarType>();
          }
          return Cinv;
        }();

        //--------------------------------------------------------------------
        // 1. Read vertices
        //--------------------------------------------------------------------
        const size_t nVertices = mesh.getConnectivity().getCount(0);
        std::vector<Index> vertexScalarDof(nVertices);
        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
        {
          const auto& vdofs = fes.getDOFs(0, v);
          assert(vdofs.size() >= 1 && "H1 vertex should have at least one DOF.");
          vertexScalarDof[static_cast<size_t>(v)] = to_scalar_dof(vdofs(0));
        }

        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
        {
          const Index sdof = vertexScalarDof[static_cast<size_t>(v)];
          const size_t s = static_cast<size_t>(sdof);
          if (s < scalarSize && !written[s])
          {
            set_scalar_dof_from_pos(sdof);
            written[s] = true;
          }
        }

        //--------------------------------------------------------------------
        // 2. Read edges: interior DOFs, oriented vmin -> vmax
        //--------------------------------------------------------------------
        if (D >= 1)
        {
          const auto& conn10  = mesh.getConnectivity().getIncidence(1, 0);
          const size_t nEdges = mesh.getConnectivity().getCount(1);

          std::vector<Index> interior;
          for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
          {
            const auto& edgeVerts = conn10[e];
            assert(edgeVerts.size() == 2);

            const Index v0   = edgeVerts[0];
            const Index v1   = edgeVerts[1];
            const Index vmin = std::min(v0, v1);
            const Index vmax = std::max(v0, v1);

            const Index vminDof = vertexScalarDof[static_cast<size_t>(vmin)];
            const Index vmaxDof = vertexScalarDof[static_cast<size_t>(vmax)];

            const auto& edofs = fes.getDOFs(1, e);

            interior.clear();
            for (Index k = 0; k < static_cast<Index>(edofs.size()); ++k)
            {
              const Index sd = to_scalar_dof(edofs(k));
              if (sd != vminDof && sd != vmaxDof)
                interior.push_back(sd);
            }

            if (v0 > v1)
              std::reverse(interior.begin(), interior.end());

            for (Index sdof : interior)
            {
              const size_t s = static_cast<size_t>(sdof);
              if (s >= scalarSize || written[s])
                continue;

              set_scalar_dof_from_pos(sdof);
              written[s] = true;
            }
          }
        }

        //--------------------------------------------------------------------
        // 3. Read faces (D >= 2)
        //--------------------------------------------------------------------
        if (D >= 2)
        {
          const size_t faceDim   = (D == 3) ? 2 : (D - 1);
          const size_t faceCount = mesh.getConnectivity().getCount(faceDim);

          if (D == 3)
          {
            // Triangle face parameters
            constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
            const int p  = static_cast<int>(K);
            const int nV = 3;
            const int nE = 3 * (p - 1);
            const int triInteriorOffset = nV + nE;

            Math::Vector<ScalarType> uM_face(TriN);
            std::vector<Math::Vector<ScalarType>> uR_face(vdim, Math::Vector<ScalarType>(TriN));

            for (Index f = 0; f < static_cast<Index>(faceCount); ++f)
            {
              const auto faceGeom = mesh.getGeometry(2, f);
              const auto& fdofs   = fes.getDOFs(2, f);

              switch (faceGeom)
              {
                case Geometry::Polytope::Type::Triangle:
                {
                  assert(static_cast<size_t>(fdofs.size()) == TriN);
                  const int numInterior = static_cast<int>(TriN) - triInteriorOffset;

                  if (numInterior > 0)
                  {
                    for (size_t comp = 0; comp < vdim; ++comp)
                    {
                      uM_face.setZero();

                      Math::Vector<ScalarType> temp_uR(TriN);
                      for (size_t k = 0; k < TriN; ++k)
                      {
                        const Index sd = to_scalar_dof(fdofs(static_cast<Index>(k)));
                        temp_uR(static_cast<Index>(k)) =
                          (static_cast<size_t>(sd) < scalarSize)
                            ? data.coeffRef(sd + static_cast<Index>(comp * scalarSize))
                            : ScalarType(0);
                      }

                      Math::Vector<ScalarType> temp_uM = s_tri_change_scalar * temp_uR;

                      for (int k = 0; k < triInteriorOffset; ++k)
                        uM_face(k) = temp_uM(k);

                      // IMPORTANT: consume from MFEM scalar stream position (pos),
                      // not from mfem_values directly
                      for (int k = 0; k < numInterior; ++k)
                      {
                        assert(pos < scalarSize);
                        uM_face(triInteriorOffset + k) = mfem_block[comp][pos];
                        ++pos;
                      }

                      uR_face[comp] = s_tri_inv_change_scalar * uM_face;
                    }

                    for (size_t k = 0; k < TriN; ++k)
                    {
                      const Index sd = to_scalar_dof(fdofs(static_cast<Index>(k)));
                      const size_t s = static_cast<size_t>(sd);
                      if (s >= scalarSize)
                        continue;

                      for (size_t comp = 0; comp < vdim; ++comp)
                        data.coeffRef(sd + static_cast<Index>(comp * scalarSize)) = uR_face[comp](static_cast<Index>(k));
                      written[s] = true;
                    }
                  }
                  else
                  {
                    for (size_t k = 0; k < TriN; ++k)
                    {
                      const Index sd = to_scalar_dof(fdofs(static_cast<Index>(k)));
                      const size_t s = static_cast<size_t>(sd);
                      if (s < scalarSize)
                        written[s] = true;
                    }
                  }
                  break;
                }

                default:
                {
                  for (Index k = 0; k < static_cast<Index>(fdofs.size()); ++k)
                  {
                    const Index sdof = to_scalar_dof(fdofs(k));
                    const size_t s = static_cast<size_t>(sdof);
                    if (s >= scalarSize || written[s])
                      continue;

                    set_scalar_dof_from_pos(sdof);
                    written[s] = true;
                  }
                  break;
                }
              }
            }
          }
        }

        //--------------------------------------------------------------------
        // 4. Read element interiors
        //--------------------------------------------------------------------
        const size_t nCells = mesh.getConnectivity().getCount(D);

        if (D == 2)
        {
          constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
          const int p  = static_cast<int>(K);
          const int nV = 3;
          const int nE = 3 * (p - 1);
          const int triInteriorOffset = nV + nE;

          Math::Vector<ScalarType> uM_elem(TriN);
          std::vector<Math::Vector<ScalarType>> uR_elem(vdim, Math::Vector<ScalarType>(TriN));

          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto geom = mesh.getGeometry(2, c);
            const auto& cdofs = fes.getDOFs(2, c);

            switch (geom)
            {
              case Geometry::Polytope::Type::Triangle:
              {
                assert(static_cast<size_t>(cdofs.size()) == TriN);
                const int numInterior = static_cast<int>(TriN) - triInteriorOffset;

                if (numInterior > 0)
                {
                  for (size_t comp = 0; comp < vdim; ++comp)
                  {
                    uM_elem.setZero();

                    Math::Vector<ScalarType> temp_uR(TriN);
                    for (size_t k = 0; k < TriN; ++k)
                    {
                      const Index sd = to_scalar_dof(cdofs(static_cast<Index>(k)));
                      temp_uR(static_cast<Index>(k)) =
                        (static_cast<size_t>(sd) < scalarSize)
                          ? data.coeffRef(sd + static_cast<Index>(comp * scalarSize))
                          : ScalarType(0);
                    }

                    Math::Vector<ScalarType> temp_uM = s_tri_change_scalar * temp_uR;

                    for (int k = 0; k < triInteriorOffset; ++k)
                      uM_elem(k) = temp_uM(k);

                    for (int k = 0; k < numInterior; ++k)
                    {
                      assert(pos < scalarSize);
                      uM_elem(triInteriorOffset + k) = mfem_block[comp][pos];
                      ++pos;
                    }

                    uR_elem[comp] = s_tri_inv_change_scalar * uM_elem;
                  }

                  for (size_t k = 0; k < TriN; ++k)
                  {
                    const Index sd = to_scalar_dof(cdofs(static_cast<Index>(k)));
                    const size_t s = static_cast<size_t>(sd);
                    if (s >= scalarSize)
                      continue;

                    for (size_t comp = 0; comp < vdim; ++comp)
                      data.coeffRef(sd + static_cast<Index>(comp * scalarSize)) = uR_elem[comp](static_cast<Index>(k));
                    written[s] = true;
                  }
                }
                else
                {
                  for (size_t k = 0; k < TriN; ++k)
                  {
                    const Index sd = to_scalar_dof(cdofs(static_cast<Index>(k)));
                    const size_t s = static_cast<size_t>(sd);
                    if (s < scalarSize)
                      written[s] = true;
                  }
                }
                break;
              }

              default:
              {
                for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
                {
                  const Index sdof = to_scalar_dof(cdofs(k));
                  const size_t s = static_cast<size_t>(sdof);
                  if (s >= scalarSize || written[s])
                    continue;

                  set_scalar_dof_from_pos(sdof);
                  written[s] = true;
                }
                break;
              }
            }
          }
        }
        else if (D == 3)
        {
          const auto& conn30 = mesh.getConnectivity().getIncidence(3, 0);

          constexpr size_t TetN = MFEM::TetrahedronNodes<K>::Count;
          const int p = static_cast<int>(K);

          const int nV    = 4;
          const int nE    = 6 * (p - 1);
          const int nF    = 2 * (p - 1) * (p - 2);
          const int tetInteriorOffset = nV + nE + nF;

          Math::Vector<ScalarType> uM_elem(TetN);
          std::vector<Math::Vector<ScalarType>> uR_elem(vdim, Math::Vector<ScalarType>(TetN));

          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto  cellGeom = mesh.getGeometry(3, c);
            const auto& cellVer  = conn30[c];
            const auto& cdofs    = fes.getDOFs(3, c);

            switch (cellGeom)
            {
              case Geometry::Polytope::Type::Tetrahedron:
              {
                assert(cellVer.size() == 4);
                assert(static_cast<size_t>(cdofs.size()) == TetN);

                const int numInterior = static_cast<int>(TetN) - tetInteriorOffset;

                if (numInterior > 0)
                {
                  for (size_t comp = 0; comp < vdim; ++comp)
                  {
                    uM_elem.setZero();

                    Math::Vector<ScalarType> temp_uR(TetN);
                    for (size_t k = 0; k < TetN; ++k)
                    {
                      const Index sd = to_scalar_dof(cdofs(static_cast<Index>(k)));
                      temp_uR(static_cast<Index>(k)) =
                        (static_cast<size_t>(sd) < scalarSize)
                          ? data.coeffRef(sd + static_cast<Index>(comp * scalarSize))
                          : ScalarType(0);
                    }

                    Math::Vector<ScalarType> temp_uM = s_tet_change_scalar * temp_uR;

                    for (int k = 0; k < tetInteriorOffset; ++k)
                      uM_elem(k) = temp_uM(k);

                    for (int k = 0; k < numInterior; ++k)
                    {
                      assert(pos < scalarSize);
                      uM_elem(tetInteriorOffset + k) = mfem_block[comp][pos];
                      ++pos;
                    }

                    uR_elem[comp] = s_tet_inv_change_scalar * uM_elem;
                  }

                  for (size_t k = 0; k < TetN; ++k)
                  {
                    const Index sd = to_scalar_dof(cdofs(static_cast<Index>(k)));
                    const size_t s = static_cast<size_t>(sd);
                    if (s >= scalarSize)
                      continue;

                    for (size_t comp = 0; comp < vdim; ++comp)
                      data.coeffRef(sd + static_cast<Index>(comp * scalarSize)) = uR_elem[comp](static_cast<Index>(k));
                    written[s] = true;
                  }
                }
                else
                {
                  for (size_t k = 0; k < TetN; ++k)
                  {
                    const Index sd = to_scalar_dof(cdofs(static_cast<Index>(k)));
                    const size_t s = static_cast<size_t>(sd);
                    if (s < scalarSize)
                      written[s] = true;
                  }
                }
                break;
              }

              case Geometry::Polytope::Type::Wedge:
              {
                assert(false && "Wedge elements not supported in MFEM I/O");
                break;
              }

              default:
              {
                for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
                {
                  const Index sdof = to_scalar_dof(cdofs(k));
                  const size_t s = static_cast<size_t>(sdof);
                  if (s >= scalarSize || written[s])
                    continue;

                  set_scalar_dof_from_pos(sdof);
                  written[s] = true;
                }
                break;
              }
            }
          }
        }
        else
        {
          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto& cdofs = fes.getDOFs(D, c);
            for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
            {
              const Index sdof = to_scalar_dof(cdofs(k));
              const size_t s = static_cast<size_t>(sdof);
              if (s >= scalarSize || written[s])
                continue;

              set_scalar_dof_from_pos(sdof);
              written[s] = true;
            }
          }
        }

        // We must consume exactly scalarSize scalar positions (each position includes all components)
        assert(pos == scalarSize && "Not all scalar coefficients were consumed");
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
        (void) rfes;
        assert(rfes && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pfec  =
          boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[get_fec];
        const bool rfec  = boost::spirit::x3::phrase_parse(it, line.end(), pfec, space);
        (void) rfec;
        assert(rfec && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pvdim = boost::spirit::x3::string("VDim:") >> uint_[get_vdim];
        const bool rvdim = boost::spirit::x3::phrase_parse(it, line.end(), pvdim, space);
        (void) rvdim;
        assert(rvdim && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pord =
          boost::spirit::x3::string("Ordering:") >> uint_[get_ordering];
        const bool rord = boost::spirit::x3::phrase_parse(it, line.end(), pord, space);
        (void) rord;
        assert(rord && it == line.end());

        // -------------------------------------------------------------
        // 2. Read coefficient data
        // -------------------------------------------------------------
        auto& gf  = this->getObject();
        auto& fes = gf.getFiniteElementSpace();
        auto& data = gf.getData();

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
          (void) rvalue;
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
           << "Ordering: " << MFEM::Ordering::VectorDimension
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
        // Set maximum precision for floating-point output to avoid precision loss
        os << std::setprecision(std::numeric_limits<Scalar>::max_digits10);

        const auto& gf   = this->getObject();
        const auto& fes  = gf.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto& data = gf.getData();

        const size_t vdim       = fes.getVectorDimension();
        const size_t D          = mesh.getDimension();
        const size_t scalarSize = fes.getSize() / vdim;

        // Track which "Rodin scalar DOFs" have already been emitted.
        std::vector<uint8_t> written(scalarSize, false);

        const auto to_scalar_dof = [&](Index dof) -> Index
        {
          if (vdim > 1 && dof >= scalarSize)
            return dof % scalarSize;
          return dof;
        };

        const auto emit_scalar_dof = [&](Index rodin_dof)
        {
          const Index scalar_dof = to_scalar_dof(rodin_dof);
          const size_t s = static_cast<size_t>(scalar_dof);
          if (s >= scalarSize)
            return;
          if (written[s])
            return;
          for (size_t c = 0; c < vdim; ++c)
            os << data.coeffRef(scalar_dof + c * scalarSize) << '\n';
          written[s] = true;
        };

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

        //--------------------------------------------------------------------
        // Precomputed change-of-nodes matrices (Rodin Fekete -> MFEM nodes)
        //--------------------------------------------------------------------

        auto& s_tri_change_scalar = []() -> const Math::Matrix<Scalar>&
        {
          static thread_local Math::Matrix<Scalar> C;
          if (C.size() == 0)
          {
            const auto& V_mfem   = MFEM::VandermondeTriangle<K>::getMatrix();
            const auto& V_rodInv = Variational::VandermondeTriangle<K>::getInverse();
            const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
            C = C_real.template cast<Scalar>();
          }
          return C;
        }();

        auto& s_tet_change_scalar = []() -> const Math::Matrix<Scalar>&
        {
          static thread_local Math::Matrix<Scalar> C;
          if (C.size() == 0)
          {
            const auto& V_mfem   = MFEM::VandermondeTetrahedron<K>::getMatrix();
            const auto& V_rodInv = Variational::VandermondeTetrahedron<K>::getInverse();
            const Math::Matrix<Real> C_real = V_mfem * V_rodInv;
            C = C_real.template cast<Scalar>();
          }
          return C;
        }();

        //--------------------------------------------------------------------
        // 1. Vertices
        //--------------------------------------------------------------------

        const size_t nVertices = mesh.getConnectivity().getCount(0);
        std::vector<Index> vertexScalarDof(nVertices);
        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
        {
          const auto& vdofs = fes.getDOFs(0, v);
          assert(vdofs.size() >= 1 && "H1 vertex should have at least one DOF.");
          vertexScalarDof[v] = to_scalar_dof(vdofs(0));
        }

        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
          emit_scalar_dof(vertexScalarDof[v]);

        //--------------------------------------------------------------------
        // 2. Edges: interior DOFs, oriented vmin -> vmax
        //--------------------------------------------------------------------

        if (D >= 1)
        {
          const auto& conn10  = mesh.getConnectivity().getIncidence(1, 0);
          const size_t nEdges = mesh.getConnectivity().getCount(1);

          std::vector<Index> interior;
          for (Index e = 0; e < static_cast<Index>(nEdges); ++e)
          {
            const auto& edgeVerts = conn10[e];
            assert(edgeVerts.size() == 2);

            const Index v0   = edgeVerts[0];
            const Index v1   = edgeVerts[1];
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

            if (v0 > v1)
              std::reverse(interior.begin(), interior.end());

            for (Index d : interior)
              emit_scalar_dof(d);
          }
        }

        //--------------------------------------------------------------------
        // 3. Faces (D >= 2)
        //--------------------------------------------------------------------

        if (D >= 2)
        {
          const size_t faceDim   = (D == 3) ? 2 : (D - 1);
          const size_t faceCount = mesh.getConnectivity().getCount(faceDim);

          if (D == 3)
          {
            const auto& conn20 = mesh.getConnectivity().getIncidence(2, 0);
            const auto& conn32 = mesh.getConnectivity().getIncidence(3, 2);

            // Triangle face size / offsets in MFEM ordering
            constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
            const int p  = static_cast<int>(K);
            const int nV = 3;
            const int nE = 3 * (p - 1);
            const int triInteriorOffset = nV + nE;

            Math::Vector<Scalar> uR_face(TriN);
            std::vector<Math::Vector<Scalar>> uM_face(
              vdim, Math::Vector<Scalar>(TriN));

            for (Index f = 0; f < static_cast<Index>(faceCount); ++f)
            {
              const auto faceGeom = mesh.getGeometry(2, f);
              const auto& fdofs   = fes.getDOFs(2, f);

              switch (faceGeom)
              {
                case Geometry::Polytope::Type::Triangle:
                {
                  assert(static_cast<size_t>(fdofs.size()) == TriN);

                  for (size_t c = 0; c < vdim; ++c)
                  {
                    for (size_t k = 0; k < TriN; ++k)
                    {
                      const Index d = fdofs(static_cast<Index>(k));
                      uR_face(static_cast<Index>(k)) =
                        data.coeffRef(d + c * scalarSize);
                    }
                    uM_face[c] = s_tri_change_scalar * uR_face;
                  }

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
                  break;
                }

                default:
                {
                  for (Index k = 0; k < static_cast<Index>(fdofs.size()); ++k)
                    emit_scalar_dof(fdofs(k));
                  break;
                }
              }
            }
          }
          else
          {
            // D == 2: faces are edges (already oriented above); keep DOF-based.
            for (Index f = 0; f < static_cast<Index>(faceCount); ++f)
            {
              const auto& fdofs = fes.getDOFs(faceDim, f);
              for (Index k = 0; k < static_cast<Index>(fdofs.size()); ++k)
                emit_scalar_dof(fdofs(k));
            }
          }
        }

        //--------------------------------------------------------------------
        // 4. Element interiors
        //--------------------------------------------------------------------

        const size_t nCells = mesh.getConnectivity().getCount(D);

        if (D == 2)
        {
          constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
          const int p  = static_cast<int>(K);
          const int nV = 3;
          const int nE = 3 * (p - 1);
          const int triInteriorOffset = nV + nE;

          Math::Vector<Scalar> uR_elem(TriN);
          std::vector<Math::Vector<Scalar>> uM_elem(
            vdim, Math::Vector<Scalar>(TriN));

          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto geom = mesh.getGeometry(2, c);
            const auto& cdofs = fes.getDOFs(2, c);

            switch (geom)
            {
              case Geometry::Polytope::Type::Triangle:
              {
                assert(static_cast<size_t>(cdofs.size()) == TriN);

                for (size_t comp = 0; comp < vdim; ++comp)
                {
                  for (size_t k = 0; k < TriN; ++k)
                  {
                    const Index d = cdofs(static_cast<Index>(k));
                    uR_elem(static_cast<Index>(k)) =
                      data.coeffRef(d + comp * scalarSize);
                  }
                  uM_elem[comp] = s_tri_change_scalar * uR_elem;
                }

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
                break;
              }

              default:
              {
                for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
                  emit_scalar_dof(cdofs(k));
                break;
              }
            }
          }
        }
        else if (D == 3)
        {
          const auto& conn30 = mesh.getConnectivity().getIncidence(3, 0);
          const auto& cp     = Variational::GLL01<K>::getNodes();

          constexpr size_t TetN = MFEM::TetrahedronNodes<K>::Count;
          const int p = static_cast<int>(K);

          const int nV    = 4;
          const int nE    = 6 * (p - 1);
          const int nF    = 2 * (p - 1) * (p - 2); // 4 faces × (p-1)(p-2)/2
          const int tetInteriorOffset = nV + nE + nF;

          Math::Vector<Scalar> uR_elem(TetN);
          std::vector<Math::Vector<Scalar>> uM_elem(
            vdim, Math::Vector<Scalar>(TetN));

          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto  cellGeom = mesh.getGeometry(3, c);
            const auto& cellVer  = conn30[c];
            const auto  cellIt   = mesh.getCell(c);
            const auto& cell     = *cellIt;
            const auto& cdofs    = fes.getDOFs(3, c);

            switch (cellGeom)
            {
              //---------------------------------------------------------------
              // Tetrahedron: change-of-nodes for interior DOFs
              //---------------------------------------------------------------
              case Geometry::Polytope::Type::Tetrahedron:
              {
                assert(cellVer.size() == 4);
                assert(static_cast<size_t>(cdofs.size()) == TetN);

                for (size_t comp = 0; comp < vdim; ++comp)
                {
                  for (size_t k = 0; k < TetN; ++k)
                  {
                    const Index d = cdofs(static_cast<Index>(k));
                    uR_elem(static_cast<Index>(k)) =
                      data.coeffRef(d + comp * scalarSize);
                  }
                  uM_elem[comp] = s_tet_change_scalar * uR_elem;
                }

                for (int idx = tetInteriorOffset; idx < static_cast<int>(TetN); ++idx)
                {
                  for (size_t comp = 0; comp < vdim; ++comp)
                    os << uM_elem[comp](static_cast<Index>(idx)) << '\n';
                }
                break;
              }

              //---------------------------------------------------------------
              // Wedge / Prism: TODO
              //---------------------------------------------------------------
              case Geometry::Polytope::Type::Wedge:
              {
                assert(false); // Unsupported
                break;
              }

              //---------------------------------------------------------------
              // Other 3D cell types (including Hexahedron): DOF-based interiors
              //---------------------------------------------------------------
              default:
              {
                for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
                  emit_scalar_dof(cdofs(k));
                break;
              }
            }
          }
        }
        else
        {
          // Other dimensions: fallback DOF-based interiors
          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto& cdofs = fes.getDOFs(D, c);
            for (Index k = 0; k < static_cast<Index>(cdofs.size()); ++k)
              emit_scalar_dof(cdofs(k));
          }
        }
      }
  };
}

#endif
