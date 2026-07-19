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
#include <algorithm>
#include <array>
#include <map>
#include <ostream>
#include <iomanip>
#include <optional>
#include <limits>
#include <vector>

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

#define RODIN_IO_MFEM_DEFAULT_POLYTOPE_ATTRIBUTE 1

namespace Rodin::IO::MFEM
{
  /**
   * @brief Keywords used in MFEM mesh file format.
   *
   * These keywords identify different sections of an MFEM mesh file.
   */
  enum class Keyword
  {
    Dimension, ///< Dimension section keyword
    Elements,  ///< Elements section keyword
    Boundary,  ///< Boundary section keyword
    Vertices   ///< Vertices section keyword
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
      case Keyword::Dimension:
        return "dimension";
      case Keyword::Elements:
        return "elements";
      case Keyword::Boundary:
        return "boundary";
      case Keyword::Vertices:
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
    if (str == Keyword::Boundary)
      res = Keyword::Boundary;
    else if (str == Keyword::Dimension)
      res = Keyword::Dimension;
    else if (str == Keyword::Elements)
      res = Keyword::Elements;
    else if (str == Keyword::Vertices)
      res = Keyword::Vertices;
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

  /// @brief Face index and vertex ordering used when emitting MFEM face DOFs.
  struct FaceOrderEntry
  {
      Index index; ///< Face index in Rodin connectivity ordering.
      std::vector<Index> vertices; ///< Face vertices in MFEM local ordering.
  };

  inline
  std::array<Index, 2> getSortedEdgeKey(Index v0, Index v1)
  {
    if (v1 < v0)
      std::swap(v0, v1);
    return { v0, v1 };
  }

  inline
  std::array<Index, 3> getTriangleFaceKey(Index v0, Index v1, Index v2)
  {
    std::array<Index, 3> key{ v0, v1, v2 };
    std::sort(key.begin(), key.end());
    return key;
  }

  inline
  std::array<Index, 3> getQuadrilateralFaceKey(Index v0, Index v1, Index v2, Index v3)
  {
    std::array<Index, 4> key{ v0, v1, v2, v3 };
    std::sort(key.begin(), key.end());
    return { key[0], key[1], key[2] };
  }

  template <class Callback>
  inline
  void forEachMFEMLocalEdge(Geometry::Polytope::Type geometry, Callback&& callback)
  {
    switch (geometry)
    {
      case Geometry::Polytope::Type::Segment:
      {
        callback(0, 1);
        break;
      }
      case Geometry::Polytope::Type::Triangle:
      {
        callback(0, 1);
        callback(1, 2);
        callback(2, 0);
        break;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        callback(0, 1);
        callback(1, 2);
        callback(2, 3);
        callback(3, 0);
        break;
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        callback(0, 1);
        callback(0, 2);
        callback(0, 3);
        callback(1, 2);
        callback(1, 3);
        callback(2, 3);
        break;
      }
      case Geometry::Polytope::Type::Pyramid:
      {
        callback(0, 1);
        callback(1, 2);
        callback(2, 3);
        callback(3, 0);
        callback(0, 4);
        callback(1, 4);
        callback(2, 4);
        callback(3, 4);
        break;
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        callback(0, 1);
        callback(1, 2);
        callback(3, 2);
        callback(0, 3);
        callback(4, 5);
        callback(5, 6);
        callback(7, 6);
        callback(4, 7);
        callback(0, 4);
        callback(1, 5);
        callback(2, 6);
        callback(3, 7);
        break;
      }
      case Geometry::Polytope::Type::Wedge:
      {
        callback(0, 1);
        callback(1, 2);
        callback(2, 0);
        callback(3, 4);
        callback(4, 5);
        callback(5, 3);
        callback(0, 3);
        callback(1, 4);
        callback(2, 5);
        break;
      }
      case Geometry::Polytope::Type::Point:
      {
        break;
      }
    }
  }

  template <class Callback>
  inline
  void forEachMFEMLocalFace(Geometry::Polytope::Type geometry, Callback&& callback)
  {
    switch (geometry)
    {
      case Geometry::Polytope::Type::Tetrahedron:
      {
        callback(std::array<int, 3>{ 1, 2, 3 });
        callback(std::array<int, 3>{ 0, 3, 2 });
        callback(std::array<int, 3>{ 0, 1, 3 });
        callback(std::array<int, 3>{ 0, 2, 1 });
        break;
      }
      case Geometry::Polytope::Type::Wedge:
      {
        callback(std::array<int, 3>{ 0, 2, 1 });
        callback(std::array<int, 3>{ 3, 4, 5 });
        callback(std::array<int, 4>{ 0, 1, 4, 3 });
        callback(std::array<int, 4>{ 1, 2, 5, 4 });
        callback(std::array<int, 4>{ 2, 0, 3, 5 });
        break;
      }
      case Geometry::Polytope::Type::Pyramid:
      {
        callback(std::array<int, 4>{ 0, 1, 2, 3 });
        callback(std::array<int, 3>{ 0, 1, 4 });
        callback(std::array<int, 3>{ 1, 2, 4 });
        callback(std::array<int, 3>{ 2, 3, 4 });
        callback(std::array<int, 3>{ 3, 0, 4 });
        break;
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        callback(std::array<int, 4>{ 3, 2, 1, 0 });
        callback(std::array<int, 4>{ 0, 1, 5, 4 });
        callback(std::array<int, 4>{ 1, 2, 6, 5 });
        callback(std::array<int, 4>{ 2, 3, 7, 6 });
        callback(std::array<int, 4>{ 3, 0, 4, 7 });
        callback(std::array<int, 4>{ 4, 5, 6, 7 });
        break;
      }
      default:
      {
        break;
      }
    }
  }

  inline
  std::vector<Index> getMFEMEdgeOrder(const Geometry::Mesh<Context::Local>& mesh)
  {
    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();
    const size_t edgeCount = conn.getCount(1);

    std::vector<Index> out;
    out.reserve(edgeCount);

    if (D == 1)
    {
      for (Index e = 0; e < static_cast<Index>(edgeCount); ++e)
        out.push_back(e);
      return out;
    }

    std::map<std::array<Index, 2>, Index> edgeIndex;
    const auto& conn10 = conn.getIncidence(1, 0);
    for (Index e = 0; e < static_cast<Index>(edgeCount); ++e)
    {
      const auto& vertices = conn10[e];
      assert(vertices.size() == 2);
      edgeIndex.emplace(getSortedEdgeKey(vertices[0], vertices[1]), e);
    }

    std::map<std::array<Index, 2>, bool> seen;
    const size_t cellCount = conn.getCount(D);
    for (Index c = 0; c < static_cast<Index>(cellCount); ++c)
    {
      const auto geometry = conn.getGeometry(D, c);
      const auto& cell = conn.getPolytope(D, c);
      forEachMFEMLocalEdge(geometry,
        [&](int a, int b)
        {
          assert(a < static_cast<int>(cell.size()));
          assert(b < static_cast<int>(cell.size()));
          const auto key = getSortedEdgeKey(cell[a], cell[b]);
          const auto [seenIt, inserted] = seen.emplace(key, true);
          (void) seenIt;
          if (!inserted)
            return;
          const auto edgeIt = edgeIndex.find(key);
          assert(edgeIt != edgeIndex.end());
          out.push_back(edgeIt->second);
        });
    }

    assert(out.size() == edgeCount);
    return out;
  }

  inline
  std::vector<FaceOrderEntry> getMFEMFaceOrder(const Geometry::Mesh<Context::Local>& mesh)
  {
    const auto& conn = mesh.getConnectivity();
    const size_t D = mesh.getDimension();
    assert(D == 3);

    const size_t faceCount = conn.getCount(2);
    std::vector<FaceOrderEntry> out;
    out.reserve(faceCount);

    std::map<std::array<Index, 3>, Index> faceIndex;
    for (Index f = 0; f < static_cast<Index>(faceCount); ++f)
    {
      const auto& face = conn.getPolytope(2, f);
      if (face.size() == 3)
        faceIndex.emplace(getTriangleFaceKey(face[0], face[1], face[2]), f);
      else if (face.size() == 4)
        faceIndex.emplace(getQuadrilateralFaceKey(face[0], face[1], face[2], face[3]), f);
      else
        assert(false && "Unsupported face arity in MFEM face order.");
    }

    std::map<std::array<Index, 3>, bool> seen;
    const size_t cellCount = conn.getCount(D);
    for (Index c = 0; c < static_cast<Index>(cellCount); ++c)
    {
      const auto geometry = conn.getGeometry(D, c);
      const auto& cell = conn.getPolytope(D, c);

      forEachMFEMLocalFace(geometry,
        [&](const auto& localFace)
        {
          std::vector<Index> vertices;
          vertices.reserve(localFace.size());
          for (int localVertex : localFace)
          {
            assert(localVertex >= 0);
            assert(localVertex < static_cast<int>(cell.size()));
            vertices.push_back(cell[localVertex]);
          }

          std::array<Index, 3> key;
          if (localFace.size() == 3)
            key = getTriangleFaceKey(vertices[0], vertices[1], vertices[2]);
          else
            key = getQuadrilateralFaceKey(vertices[0], vertices[1], vertices[2], vertices[3]);

          const auto [seenIt, inserted] = seen.emplace(key, true);
          (void) seenIt;
          if (!inserted)
            return;

          const auto faceIt = faceIndex.find(key);
          assert(faceIt != faceIndex.end());
          out.push_back({ faceIt->second, std::move(vertices) });
        });
    }

    assert(out.size() == faceCount);
    return out;
  }

  template <size_t K>
  inline
  size_t getTriangleNodeIndex(size_t bary0, size_t bary1, size_t bary2)
  {
    assert(bary0 + bary1 + bary2 == K);
    const size_t i = bary1;
    const size_t j = bary2;
    return j * (K + 1) - (j * (j - 1)) / 2 + i;
  }

  template <size_t K>
  inline
  std::vector<size_t> getTriangleNodeMap(
      const Geometry::Polytope::Key& sourceVertices,
      const std::vector<Index>& targetVertices)
  {
    assert(sourceVertices.size() == 3);
    assert(targetVertices.size() == 3);

    constexpr size_t Count = (K + 1) * (K + 2) / 2;
    std::vector<size_t> map(Count);

    std::array<int, 3> targetToSource{};
    for (size_t t = 0; t < 3; ++t)
    {
      targetToSource[t] = -1;
      for (size_t s = 0; s < 3; ++s)
      {
        if (sourceVertices[s] == targetVertices[t])
        {
          targetToSource[t] = static_cast<int>(s);
          break;
        }
      }
      assert(targetToSource[t] >= 0);
    }

    for (size_t j = 0; j <= K; ++j)
    {
      for (size_t i = 0; i <= K - j; ++i)
      {
        const std::array<size_t, 3> targetBary{ K - i - j, i, j };
        std::array<size_t, 3> sourceBary{ 0, 0, 0 };
        for (size_t t = 0; t < 3; ++t)
          sourceBary[static_cast<size_t>(targetToSource[t])] = targetBary[t];

        const size_t target = getTriangleNodeIndex<K>(targetBary[0], targetBary[1], targetBary[2]);
        const size_t source = getTriangleNodeIndex<K>(sourceBary[0], sourceBary[1], sourceBary[2]);
        map[target] = source;
      }
    }

    return map;
  }

  template <size_t K>
  inline
  size_t getQuadrilateralNodeIndex(size_t i, size_t j)
  {
    assert(i <= K);
    assert(j <= K);
    return j * (K + 1) + i;
  }

  template <size_t K>
  inline
  std::pair<size_t, size_t> getQuadrilateralCornerCoordinates(int vertex)
  {
    switch (vertex)
    {
      case 0: return { 0, 0 };
      case 1: return { K, 0 };
      case 2: return { K, K };
      case 3: return { 0, K };
      default:
        assert(false && "Invalid quadrilateral vertex.");
        return { 0, 0 };
    }
  }

  template <size_t K>
  inline
  std::pair<size_t, size_t> applyQuadrilateralTransform(int transform, size_t i, size_t j)
  {
    switch (transform)
    {
      case 0: return { i, j };
      case 1: return { j, K - i };
      case 2: return { K - i, K - j };
      case 3: return { K - j, i };
      case 4: return { i, K - j };
      case 5: return { K - i, j };
      case 6: return { j, i };
      case 7: return { K - j, K - i };
      default:
        assert(false && "Invalid quadrilateral transform.");
        return { i, j };
    }
  }

  template <size_t K>
  inline
  std::vector<size_t> getQuadrilateralNodeMap(
      const Geometry::Polytope::Key& sourceVertices,
      const std::vector<Index>& targetVertices)
  {
    assert(sourceVertices.size() == 4);
    assert(targetVertices.size() == 4);

    constexpr size_t Count = (K + 1) * (K + 1);
    std::vector<size_t> map(Count);

    std::array<int, 4> targetToSource{};
    for (size_t t = 0; t < 4; ++t)
    {
      targetToSource[t] = -1;
      for (size_t s = 0; s < 4; ++s)
      {
        if (sourceVertices[s] == targetVertices[t])
        {
          targetToSource[t] = static_cast<int>(s);
          break;
        }
      }
      assert(targetToSource[t] >= 0);
    }

    std::array<std::pair<size_t, size_t>, 4> sourceCorners;
    for (size_t t = 0; t < 4; ++t)
      sourceCorners[t] = getQuadrilateralCornerCoordinates<K>(targetToSource[t]);

    int chosen = -1;
    for (int transform = 0; transform < 8; ++transform)
    {
      if (applyQuadrilateralTransform<K>(transform, 0, 0) == sourceCorners[0] &&
          applyQuadrilateralTransform<K>(transform, K, 0) == sourceCorners[1] &&
          applyQuadrilateralTransform<K>(transform, K, K) == sourceCorners[2] &&
          applyQuadrilateralTransform<K>(transform, 0, K) == sourceCorners[3])
      {
        chosen = transform;
        break;
      }
    }

    assert(chosen >= 0 && "Could not determine quadrilateral orientation transform.");

    for (size_t j = 0; j <= K; ++j)
    {
      for (size_t i = 0; i <= K; ++i)
      {
        const auto source = applyQuadrilateralTransform<K>(chosen, i, j);
        map[getQuadrilateralNodeIndex<K>(i, j)] =
          getQuadrilateralNodeIndex<K>(source.first, source.second);
      }
    }

    return map;
  }

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
        return Rodin::Geometry::Polytope::Type::Pyramid;
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
      case Geometry::Polytope::Type::Pyramid:
        return GeometryType::PYRAMID;
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
        const auto getDouble = [&](auto& ctx) {
          assert(i < m_sdim);
          res(i++) = _attr(ctx);
        };
        const auto p = repeat(m_sdim)[double_[getDouble]];
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
      const auto getAttribute = [&](auto& ctx) { res.attribute = _attr(ctx); };
      const auto assignGeometry = [&](auto& ctx) {
        geometry = static_cast<GeometryType>(_attr(ctx));
      };
      const auto p = uint_[getAttribute] >> uint_[assignGeometry];
      const bool r = boost::spirit::x3::phrase_parse(begin, end, p, space);

      if (!r)
        return {};
      auto g = getGeometry(geometry);
      if (!g)
        return {};
      res.geometry = *g;

     res.vertices.resize(Geometry::Polytope::Traits(res.geometry).getVertexCount());
     size_t i = 0;
     const auto getVertex = [&](auto& ctx) { res.vertices(i++) = _attr(ctx); };
     const auto pvs = repeat(res.vertices.size())[uint_[getVertex]];
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
      const auto getKeyword = [&](auto& ctx) { kw = _attr(ctx); };
      const auto p = (+alpha)[getKeyword];
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
        const auto getMajorVersion = [&](auto& ctx) { h.version.major = _attr(ctx); };
        const auto getMinorVersion = [&](auto& ctx) { h.version.minor = _attr(ctx); };
        const auto p = boost::spirit::x3::string("MFEM") >>
          boost::spirit::x3::string("mesh") >> boost::spirit::x3::char_('v') >>
          uint_[getMajorVersion] >> boost::spirit::x3::char_('.') >>
          uint_[getMinorVersion];
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
      /// @brief Number of MFEM H1 nodes on a triangle of order @p K.
      static constexpr size_t Count = (K + 1) * (K + 2) / 2;

      /// @brief Returns the MFEM reference nodes for the triangle.
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
      /// @brief Returns the Dubiner Vandermonde matrix evaluated at MFEM triangle nodes.
      static const Math::Matrix<Real>& getMatrix()
      {
        static thread_local Math::Matrix<Real> s_vandermonde;

        constexpr size_t N = TriangleNodes<K>::Count;

        if (s_vandermonde.size() == 0)
        {
          const auto& nodes = TriangleNodes<K>::getNodes();
          s_vandermonde.resize(N, N);

          size_t modeIdx = 0;
          Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
            constexpr size_t P = pIdx.value;
            Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
              constexpr size_t Q = qIdx.value;
              for (size_t nodeIdx = 0; nodeIdx < N; ++nodeIdx)
              {
                const auto& pt = nodes[nodeIdx];
                const Real x = pt.x();
                const Real y = pt.y();

                Real r, s;
                Variational::DubinerTriangle<K>::getCollapsed(r, s, x, y);

                Variational::DubinerTriangle<K>::template getBasis<P, Q>(
                  s_vandermonde(nodeIdx, modeIdx), r, s);
              }
              ++modeIdx;
            });
          });
        }

        return s_vandermonde;
      }

      /// @brief Returns the inverse of @ref getMatrix().
      static const Math::Matrix<Real>& getInverse()
      {
        static thread_local Math::Matrix<Real> s_inv;

        if (s_inv.size() == 0)
        {
          const auto& V = VandermondeTriangle<K>::getMatrix();
          Math::ThinSVD<Math::Matrix<Real>> svd(V);
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
      /// @brief Number of MFEM H1 nodes on a tetrahedron of order @p K.
      static constexpr size_t Count = (K + 1) * (K + 2) * (K + 3) / 6;

      /// @brief Returns the MFEM reference nodes for the tetrahedron.
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
      /// @brief Returns the Dubiner Vandermonde matrix evaluated at MFEM tetrahedron nodes.
      static const Math::Matrix<Real>& getMatrix()
      {
        static thread_local Math::Matrix<Real> s_vandermonde;

        constexpr size_t N = TetrahedronNodes<K>::Count;

        if (s_vandermonde.size() == 0)
        {
          const auto& nodes = TetrahedronNodes<K>::getNodes();
          s_vandermonde.resize(N, N);

          size_t modeIdx = 0;
          Rodin::Utility::ForIndex<K + 1>([&](auto pIdx) {
            constexpr size_t P = pIdx.value;
            Rodin::Utility::ForIndex<K + 1 - P>([&](auto qIdx) {
              constexpr size_t Q = qIdx.value;
              Rodin::Utility::ForIndex<K + 1 - P - Q>([&](auto rIdx) {
                constexpr size_t R = rIdx.value;
                for (size_t nodeIdx = 0; nodeIdx < N; ++nodeIdx)
                {
                  const auto& pt = nodes[nodeIdx];
                  const Real x = pt.x();
                  const Real y = pt.y();
                  const Real z = pt.z();

                  Real a, b, c;
                  Variational::DubinerTetrahedron<K>::getCollapsed(a, b, c, x, y, z);

                  Variational::DubinerTetrahedron<K>::template getBasis<P, Q, R>(
                    s_vandermonde(nodeIdx, modeIdx), a, b, c);
                }
                ++modeIdx;
              });
            });
          });
        }

        return s_vandermonde;
      }

      /// @brief Returns the inverse of @ref getMatrix().
      static const Math::Matrix<Real>& getInverse()
      {
        static thread_local Math::Matrix<Real> s_inv;

        if (s_inv.size() == 0)
        {
          const auto& V = VandermondeTetrahedron<K>::getMatrix();
          Math::ThinSVD<Math::Matrix<Real>> svd(V);
          const Math::Matrix<Real> I = Math::Matrix<Real>::Identity(V.rows(), V.cols());
          s_inv = svd.solve(I);
        }

        return s_inv;
      }
  };

  /**
   * @brief MFEM-style nodal set on the reference wedge / prism.
   *
   * This reproduces the local node ordering used by mfem::H1_WedgeElement:
   * vertices, edge interiors, triangular face interiors, quadrilateral face
   * interiors, then cell interiors. The triangular factor uses MFEM's closed
   * H1 triangle nodes and the segment factor uses MFEM's H1 segment ordering.
   */
  template <size_t K>
  class WedgeNodes
  {
    public:
      /// @brief Number of triangular nodes per wedge layer.
      static constexpr size_t TriangleCount = TriangleNodes<K>::Count;
      /// @brief Number of MFEM H1 nodes on a wedge of order @p K.
      static constexpr size_t Count = TriangleCount * (K + 1);

      /// @brief Returns the MFEM reference nodes for the wedge.
      static const std::vector<Math::SpatialPoint>& getNodes()
      {
        static thread_local std::vector<Math::SpatialPoint> s_nodes;

        if (s_nodes.empty())
        {
          constexpr int p = static_cast<int>(K);
          constexpr int ne = p - 1;
          constexpr int nt = (p - 1) * (p - 2) / 2;
          constexpr int nq = (p - 1) * (p - 1);

          const auto& triNodes = TriangleNodes<K>::getNodes();
          const auto& cp = Variational::GLL01<K>::getNodes();

          auto segmentNode = [&](int local) -> Real
          {
            if (local == 0)
              return cp[0];
            if (local == 1)
              return cp[p];
            return cp[local - 1];
          };

          std::vector<int> tDof(Count, -1);
          std::vector<int> sDof(Count, -1);

          tDof[0] = 0; sDof[0] = 0;
          tDof[1] = 1; sDof[1] = 0;
          tDof[2] = 2; sDof[2] = 0;
          tDof[3] = 0; sDof[3] = 1;
          tDof[4] = 1; sDof[4] = 1;
          tDof[5] = 2; sDof[5] = 1;

          int offset = 6;
          for (int i = 1; i < p; ++i)
          {
            const int k = i - 1;

            tDof[offset + 0 * ne + k] = 2 + 0 * ne + i; sDof[offset + 0 * ne + k] = 0;
            tDof[offset + 1 * ne + k] = 2 + 1 * ne + i; sDof[offset + 1 * ne + k] = 0;
            tDof[offset + 2 * ne + k] = 2 + 2 * ne + i; sDof[offset + 2 * ne + k] = 0;

            tDof[offset + 3 * ne + k] = 2 + 0 * ne + i; sDof[offset + 3 * ne + k] = 1;
            tDof[offset + 4 * ne + k] = 2 + 1 * ne + i; sDof[offset + 4 * ne + k] = 1;
            tDof[offset + 5 * ne + k] = 2 + 2 * ne + i; sDof[offset + 5 * ne + k] = 1;

            tDof[offset + 6 * ne + k] = 0; sDof[offset + 6 * ne + k] = i + 1;
            tDof[offset + 7 * ne + k] = 1; sDof[offset + 7 * ne + k] = i + 1;
            tDof[offset + 8 * ne + k] = 2; sDof[offset + 8 * ne + k] = i + 1;
          }

          offset += 9 * ne;

          int k = 0;
          for (int j = 1; j < p; ++j)
          {
            for (int i = 1; i < p - j; ++i)
            {
              const int l = j - p + (((2 * p - 1) - i) * i) / 2;
              tDof[offset + k]      = 3 * p + l; sDof[offset + k]      = 0;
              tDof[offset + nt + k] = 3 * p + k; sDof[offset + nt + k] = 1;
              ++k;
            }
          }

          offset += 2 * nt;

          k = 0;
          for (int j = 1; j < p; ++j)
          {
            for (int i = 1; i < p; ++i)
            {
              tDof[offset + 0 * nq + k] = 2 + 0 * ne + i;
              tDof[offset + 1 * nq + k] = 2 + 1 * ne + i;
              tDof[offset + 2 * nq + k] = 2 + 2 * ne + i;

              sDof[offset + 0 * nq + k] = 1 + j;
              sDof[offset + 1 * nq + k] = 1 + j;
              sDof[offset + 2 * nq + k] = 1 + j;
              ++k;
            }
          }

          offset += 3 * nq;

          int m = 0;
          for (int z = 1; z < p; ++z)
          {
            int l = 0;
            for (int j = 1; j < p; ++j)
            {
              for (int i = 1; i + j < p; ++i)
              {
                tDof[offset + m] = 3 * p + l;
                sDof[offset + m] = 1 + z;
                ++l;
                ++m;
              }
            }
          }

          s_nodes.reserve(Count);
          for (size_t i = 0; i < Count; ++i)
          {
            assert(tDof[i] >= 0 && sDof[i] >= 0);
            const auto& pt = triNodes[static_cast<size_t>(tDof[i])];
            s_nodes.emplace_back(Math::SpatialPoint{{ pt.x(), pt.y(), segmentNode(sDof[i]) }});
          }

          assert(s_nodes.size() == Count);
        }

        return s_nodes;
      }
  };

  /**
   * @brief Change of nodes from Rodin wedge nodes to MFEM wedge nodes.
   *
   * C(i,j) = phi_j(x_i), where phi_j is Rodin's wedge H1 nodal basis and x_i
   * is the i-th MFEM wedge node.
   */
  template <size_t K>
  class WedgeChange
  {
    public:
      /// @brief Returns the change-of-nodes matrix from Rodin to MFEM wedge nodes.
      static const Math::Matrix<Real>& getMatrix()
      {
        static thread_local Math::Matrix<Real> s_change;

        constexpr size_t N = WedgeNodes<K>::Count;

        if (s_change.size() == 0)
        {
          s_change.resize(N, N);
          const auto& nodes = WedgeNodes<K>::getNodes();
          const Variational::H1Element<K, Real> fe(Geometry::Polytope::Type::Wedge);

          for (size_t i = 0; i < N; ++i)
          {
            for (size_t j = 0; j < N; ++j)
              s_change(static_cast<Index>(i), static_cast<Index>(j)) = fe.getBasis(j)(nodes[i]);
          }
        }

        return s_change;
      }

      /// @brief Returns the inverse of @ref getMatrix().
      static const Math::Matrix<Real>& getInverse()
      {
        static thread_local Math::Matrix<Real> s_inv;

        if (s_inv.size() == 0)
        {
          const auto& C = WedgeChange<K>::getMatrix();
          Math::ThinSVD<Math::Matrix<Real>> svd(C);
          const Math::Matrix<Real> I = Math::Matrix<Real>::Identity(C.rows(), C.cols());
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
      /// @brief Execution context type.
      using ContextType = Context::Local;

      /// @brief Mesh type being loaded.
      using ObjectType = Geometry::Mesh<ContextType>;

      /// @brief Parent class type.
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
      /// @brief Execution context type.
      using ContextType = Context::Local;

      /// @brief Mesh type being printed.
      using ObjectType = Geometry::Mesh<ContextType>;

      /// @brief Parent class type.
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
        const auto getFec = [&](auto& ctx) { header.fec = _attr(ctx); };
        const auto getVdim = [&](auto& ctx) { header.vdim = _attr(ctx); };
        const auto getOrdering = [&](auto& ctx) {
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
        const auto pfec =
          boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[getFec];
        const bool rfec = boost::spirit::x3::phrase_parse(it, line.end(), pfec, space);
        (void) rfec;
        assert(it == line.end() && rfec);

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it = line.begin();
        const auto pvdim = boost::spirit::x3::string("VDim:") >> uint_[getVdim];
        const bool rvdim = boost::spirit::x3::phrase_parse(it, line.end(), pvdim, space);
        (void) rvdim;
        assert(it == line.end() && rvdim);

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it = line.begin();
        const auto pordering =
          boost::spirit::x3::string("Ordering:") >> uint_[getOrdering];
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
      /// @brief Finite element space type.
      using FESType    = Variational::H1<K, Range, Geometry::Mesh<Context::Local>>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
      /// @brief Coefficient data storage type.
      using DataType   = Math::Vector<ScalarType>;
      /// @brief Grid function type being loaded.
      /// @brief Grid function type being loaded.
      using ObjectType = Variational::GridFunction<FESType, DataType>;
      /// @brief Parent loader base type.
      using Parent     = GridFunctionLoaderBase<FESType, DataType>;

      /**
       * @brief Constructs an MFEM H1 grid-function loader.
       * @param[in,out] gf Grid function to populate.
       */
      GridFunctionLoader(ObjectType& gf)
        : Parent(gf)
      {}

      /**
       * @brief Loads H1 coefficient data from an MFEM grid-function stream.
       * @param[in,out] is Input stream containing MFEM grid-function data.
       */
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

        const auto getFec = [&](auto& ctx) { header.fec = _attr(ctx); };
        const auto getVdim = [&](auto& ctx) { header.vdim = _attr(ctx); };
        const auto getOrdering = [&](auto& ctx) {
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
          boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[getFec];
        const bool rfec = boost::spirit::x3::phrase_parse(it, line.end(), pfec, space);
        (void) rfec;
        assert(rfec && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pvdim = boost::spirit::x3::string("VDim:") >> uint_[getVdim];
        const bool rvdim = boost::spirit::x3::phrase_parse(it, line.end(), pvdim, space);
        (void) rvdim;
        assert(rvdim && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pord = boost::spirit::x3::string("Ordering:") >> uint_[getOrdering];
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
        std::vector<ScalarType> mfemValues;
        mfemValues.reserve(fes.getSize());

        ScalarType value;
        while (mfemValues.size() < fes.getSize() && (is >> value))
          mfemValues.push_back(value);

        assert(mfemValues.size() == fes.getSize() &&
          "Mismatch in number of coefficients read");

        // -------------------------------------------------------------
        // 2.a Normalize MFEM ordering into component-blocked storage:
        //      mfem_block[c][i] = value of component c at MFEM scalar position i
        // -------------------------------------------------------------
        std::vector<std::vector<ScalarType>> mfemBlock(
          vdim, std::vector<ScalarType>(scalarSize));

        if (header.ordering == MFEM::Ordering::Nodes) // 0: XXX..YYY..ZZZ..
        {
          for (size_t c = 0; c < vdim; ++c)
            for (size_t i = 0; i < scalarSize; ++i)
              mfemBlock[c][i] = mfemValues[c * scalarSize + i];
        }
        else // 1: XYZ,XYZ,XYZ...
        {
          for (size_t i = 0; i < scalarSize; ++i)
            for (size_t c = 0; c < vdim; ++c)
              mfemBlock[c][i] = mfemValues[i * vdim + c];
        }

        // Scalar-position cursor in the MFEM traversal (increments once per scalar DOF consumed)
        size_t pos = 0;

        // Helper: Rodin DOF -> scalar DOF index
        const auto toScalarDof = [&](Index dof) -> Index {
          if (vdim > 1 && dof >= static_cast<Index>(scalarSize))
            return dof % static_cast<Index>(scalarSize);
          return dof;
        };

        const auto scalarLocalSize = [&](const auto& dofs) -> size_t {
          assert(static_cast<size_t>(dofs.size()) % vdim == 0);
          return static_cast<size_t>(dofs.size()) / vdim;
        };

        const auto scalarLocalDof = [&](const auto& dofs, size_t local) -> Index {
          if (vdim == 1)
            return toScalarDof(dofs(static_cast<Index>(local)));

          assert(local * vdim < static_cast<size_t>(dofs.size()));
          return toScalarDof(dofs(static_cast<Index>(local * vdim)));
        };

        // Helper: assign (all components) for a Rodin scalar DOF from current MFEM scalar position
        auto setScalarDofFromPos = [&](Index scalarDof) {
          assert(pos < scalarSize);
          for (size_t c = 0; c < vdim; ++c)
            data.coeffRef(scalarDof + static_cast<Index>(c * scalarSize)) =
              mfemBlock[c][pos];
          ++pos;
        };

        // Track which Rodin scalar DOFs have been set
        std::vector<uint8_t> written(scalarSize, false);

        //--------------------------------------------------------------------
        // Precomputed change-of-nodes matrices (Rodin Fekete -> MFEM nodes)
        // and inverse matrices (MFEM nodes -> Rodin Fekete)
        //--------------------------------------------------------------------
        auto& sTriChangeScalar = []() -> const Math::Matrix<ScalarType>& {
          static thread_local Math::Matrix<ScalarType> C;
          if (C.size() == 0)
          {
            const auto& vMfem = MFEM::VandermondeTriangle<K>::getMatrix();
            const auto& vRodInv = Variational::VandermondeTriangle<K>::getInverse();
            const Math::Matrix<Real> cReal = vMfem * vRodInv;
            C = cReal.template cast<ScalarType>();
          }
          return C;
        }();

        auto& sTriInvChangeScalar = []() -> const Math::Matrix<ScalarType>& {
          static thread_local Math::Matrix<ScalarType> Cinv;
          if (Cinv.size() == 0)
          {
            const auto& vMfem = MFEM::VandermondeTriangle<K>::getMatrix();
            const auto& vRodInv = Variational::VandermondeTriangle<K>::getInverse();
            const Math::Matrix<Real> cReal = vMfem * vRodInv;
            const Math::Matrix<Real> cinvReal = cReal.inverse();
            Cinv = cinvReal.template cast<ScalarType>();
          }
          return Cinv;
        }();

        auto& sTetChangeScalar = []() -> const Math::Matrix<ScalarType>& {
          static thread_local Math::Matrix<ScalarType> C;
          if (C.size() == 0)
          {
            const auto& vMfem = MFEM::VandermondeTetrahedron<K>::getMatrix();
            const auto& vRodInv = Variational::VandermondeTetrahedron<K>::getInverse();
            const Math::Matrix<Real> cReal = vMfem * vRodInv;
            C = cReal.template cast<ScalarType>();
          }
          return C;
        }();

        auto& sTetInvChangeScalar = []() -> const Math::Matrix<ScalarType>& {
          static thread_local Math::Matrix<ScalarType> Cinv;
          if (Cinv.size() == 0)
          {
            const auto& vMfem = MFEM::VandermondeTetrahedron<K>::getMatrix();
            const auto& vRodInv = Variational::VandermondeTetrahedron<K>::getInverse();
            const Math::Matrix<Real> cReal = vMfem * vRodInv;
            const Math::Matrix<Real> cinvReal = cReal.inverse();
            Cinv = cinvReal.template cast<ScalarType>();
          }
          return Cinv;
        }();

        auto& sWedgeChangeScalar = []() -> const Math::Matrix<ScalarType>& {
          static thread_local Math::Matrix<ScalarType> C;
          if (C.size() == 0)
            C = MFEM::WedgeChange<K>::getMatrix().template cast<ScalarType>();
          return C;
        }();

        auto& sWedgeInvChangeScalar = []() -> const Math::Matrix<ScalarType>& {
          static thread_local Math::Matrix<ScalarType> Cinv;
          if (Cinv.size() == 0)
            Cinv = MFEM::WedgeChange<K>::getInverse().template cast<ScalarType>();
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
          vertexScalarDof[static_cast<size_t>(v)] = toScalarDof(vdofs(0));
        }

        for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
        {
          const Index sdof = vertexScalarDof[static_cast<size_t>(v)];
          const size_t s = static_cast<size_t>(sdof);
          if (s < scalarSize && !written[s])
          {
            setScalarDofFromPos(sdof);
            written[s] = true;
          }
        }

        //--------------------------------------------------------------------
        // 2. Read edges: interior DOFs, oriented vmin -> vmax
        //--------------------------------------------------------------------
        if (D >= 1)
        {
          const auto& conn10  = mesh.getConnectivity().getIncidence(1, 0);
          const auto mfemEdges = MFEM::getMFEMEdgeOrder(mesh);

          std::vector<Index> interior;
          for (Index e : mfemEdges)
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
            const size_t edgeScalarDofs = scalarLocalSize(edofs);

            interior.clear();
            for (size_t k = 0; k < edgeScalarDofs; ++k)
            {
              const Index sd = scalarLocalDof(edofs, k);
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

              setScalarDofFromPos(sdof);
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
          [[maybe_unused]] const size_t faceCount = mesh.getConnectivity().getCount(faceDim);

          if (D == 3)
          {
            // Triangle face parameters
            constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
            [[maybe_unused]] constexpr size_t QuadN = (K + 1) * (K + 1);
            const int p  = static_cast<int>(K);
            const int nV = 3;
            const int nE = 3 * (p - 1);
            const int triInteriorOffset = nV + nE;

            std::vector<Math::Vector<ScalarType>> uRFace(
              vdim, Math::Vector<ScalarType>(TriN));
            const auto mfemFaces = MFEM::getMFEMFaceOrder(mesh);

            for (const auto& mfemFace : mfemFaces)
            {
              const Index f = mfemFace.index;
              const auto faceGeom = mesh.getGeometry(2, f);
              const auto& fdofs   = fes.getDOFs(2, f);

              switch (faceGeom)
              {
                case Geometry::Polytope::Type::Triangle:
                {
                  assert(scalarLocalSize(fdofs) == TriN);
                  const auto& faceVertices = mesh.getConnectivity().getPolytope(2, f);
                  const auto faceNodeMap = MFEM::getTriangleNodeMap<K>(faceVertices, mfemFace.vertices);
                  const int numInterior = static_cast<int>(TriN) - triInteriorOffset;

                  if (numInterior > 0)
                  {
                    std::vector<Math::Vector<ScalarType>> uMByComp(
                      vdim, Math::Vector<ScalarType>(TriN));

                    for (size_t comp = 0; comp < vdim; ++comp)
                    {
                      uMByComp[comp].setZero();

                      Math::Vector<ScalarType> tempUR(TriN);
                      for (size_t k = 0; k < TriN; ++k)
                      {
                        const Index sd = scalarLocalDof(fdofs, faceNodeMap[k]);
                        tempUR(static_cast<Index>(k)) =
                          (static_cast<size_t>(sd) < scalarSize)
                          ? data.coeffRef(sd + static_cast<Index>(comp * scalarSize))
                          : ScalarType(0);
                      }

                      Math::Vector<ScalarType> tempUM = sTriChangeScalar * tempUR;

                      for (int k = 0; k < triInteriorOffset; ++k)
                        uMByComp[comp](k) = tempUM(k);
                    }

                    // IMPORTANT: consume from the MFEM scalar stream once per
                    // scalar node, then assign all components at that position.
                    for (int k = 0; k < numInterior; ++k)
                    {
                      assert(pos < scalarSize);
                      for (size_t comp = 0; comp < vdim; ++comp)
                        uMByComp[comp](triInteriorOffset + k) = mfemBlock[comp][pos];
                      ++pos;
                    }

                    for (size_t comp = 0; comp < vdim; ++comp)
                      uRFace[comp] = sTriInvChangeScalar * uMByComp[comp];

                    for (size_t k = 0; k < TriN; ++k)
                    {
                      const Index sd = scalarLocalDof(fdofs, faceNodeMap[k]);
                      const size_t s = static_cast<size_t>(sd);
                      if (s >= scalarSize)
                        continue;

                      for (size_t comp = 0; comp < vdim; ++comp)
                        data.coeffRef(sd + static_cast<Index>(comp * scalarSize)) =
                          uRFace[comp](static_cast<Index>(k));
                      written[s] = true;
                    }
                  }
                  else
                  {
                    for (size_t k = 0; k < TriN; ++k)
                    {
                      const Index sd = scalarLocalDof(fdofs, faceNodeMap[k]);
                      const size_t s = static_cast<size_t>(sd);
                      if (s < scalarSize)
                        written[s] = true;
                    }
                  }
                  break;
                }

                case Geometry::Polytope::Type::Quadrilateral:
                {
                  assert(scalarLocalSize(fdofs) == QuadN);
                  const auto& faceVertices = mesh.getConnectivity().getPolytope(2, f);
                  const auto faceNodeMap = MFEM::getQuadrilateralNodeMap<K>(faceVertices, mfemFace.vertices);

                  for (size_t j = 1; j < K; ++j)
                  {
                    for (size_t i = 1; i < K; ++i)
                    {
                      const size_t targetLocal = MFEM::getQuadrilateralNodeIndex<K>(i, j);
                      const Index sdof = scalarLocalDof(fdofs, faceNodeMap[targetLocal]);
                      const size_t s = static_cast<size_t>(sdof);
                      assert(pos < scalarSize);
                      if (s < scalarSize)
                      {
                        for (size_t comp = 0; comp < vdim; ++comp)
                          data.coeffRef(sdof + static_cast<Index>(comp * scalarSize)) =
                            mfemBlock[comp][pos];
                        written[s] = true;
                      }
                      ++pos;
                    }
                  }
                  break;
                }

                default:
                {
                  const size_t nScalarLocal = scalarLocalSize(fdofs);
                  for (size_t k = 0; k < nScalarLocal; ++k)
                  {
                    const Index sdof = scalarLocalDof(fdofs, k);
                    const size_t s = static_cast<size_t>(sdof);
                    if (s >= scalarSize || written[s])
                      continue;

                    setScalarDofFromPos(sdof);
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

          std::vector<Math::Vector<ScalarType>> uRElem(
            vdim, Math::Vector<ScalarType>(TriN));

          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto geom = mesh.getGeometry(2, c);
            const auto& cdofs = fes.getDOFs(2, c);

            switch (geom)
            {
              case Geometry::Polytope::Type::Triangle:
              {
                assert(scalarLocalSize(cdofs) == TriN);
                const int numInterior = static_cast<int>(TriN) - triInteriorOffset;

                if (numInterior > 0)
                {
                  std::vector<Math::Vector<ScalarType>> uMByComp(
                    vdim, Math::Vector<ScalarType>(TriN));

                  for (size_t comp = 0; comp < vdim; ++comp)
                  {
                    uMByComp[comp].setZero();

                    Math::Vector<ScalarType> tempUR(TriN);
                    for (size_t k = 0; k < TriN; ++k)
                    {
                      const Index sd = scalarLocalDof(cdofs, k);
                      tempUR(static_cast<Index>(k)) =
                        (static_cast<size_t>(sd) < scalarSize)
                        ? data.coeffRef(sd + static_cast<Index>(comp * scalarSize))
                        : ScalarType(0);
                    }

                    Math::Vector<ScalarType> tempUM = sTriChangeScalar * tempUR;

                    for (int k = 0; k < triInteriorOffset; ++k)
                      uMByComp[comp](k) = tempUM(k);
                  }

                  for (int k = 0; k < numInterior; ++k)
                  {
                    assert(pos < scalarSize);
                    for (size_t comp = 0; comp < vdim; ++comp)
                      uMByComp[comp](triInteriorOffset + k) = mfemBlock[comp][pos];
                    ++pos;
                  }

                  for (size_t comp = 0; comp < vdim; ++comp)
                    uRElem[comp] = sTriInvChangeScalar * uMByComp[comp];

                  for (size_t k = 0; k < TriN; ++k)
                  {
                    const Index sd = scalarLocalDof(cdofs, k);
                    const size_t s = static_cast<size_t>(sd);
                    if (s >= scalarSize)
                      continue;

                    for (size_t comp = 0; comp < vdim; ++comp)
                      data.coeffRef(sd + static_cast<Index>(comp * scalarSize)) =
                        uRElem[comp](static_cast<Index>(k));
                    written[s] = true;
                  }
                }
                else
                {
                  for (size_t k = 0; k < TriN; ++k)
                  {
                    const Index sd = scalarLocalDof(cdofs, k);
                    const size_t s = static_cast<size_t>(sd);
                    if (s < scalarSize)
                      written[s] = true;
                  }
                }
                break;
              }

              default:
              {
                const size_t nScalarLocal = scalarLocalSize(cdofs);
                for (size_t k = 0; k < nScalarLocal; ++k)
                {
                  const Index sdof = scalarLocalDof(cdofs, k);
                  const size_t s = static_cast<size_t>(sdof);
                  if (s >= scalarSize || written[s])
                    continue;

                  setScalarDofFromPos(sdof);
                  written[s] = true;
                }
                break;
              }
            }
          }
        }
        else if (D == 3)
        {
          constexpr size_t TetN = MFEM::TetrahedronNodes<K>::Count;
          const int p = static_cast<int>(K);

          const int nV    = 4;
          const int nE    = 6 * (p - 1);
          const int nF    = 2 * (p - 1) * (p - 2);
          const int tetInteriorOffset = nV + nE + nF;

          std::vector<Math::Vector<ScalarType>> uRElem(
            vdim, Math::Vector<ScalarType>(TetN));

          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto  cellGeom = mesh.getGeometry(3, c);
            const auto& cdofs    = fes.getDOFs(3, c);

            switch (cellGeom)
            {
              case Geometry::Polytope::Type::Tetrahedron:
              {
                assert(scalarLocalSize(cdofs) == TetN);

                const int numInterior = static_cast<int>(TetN) - tetInteriorOffset;

                if (numInterior > 0)
                {
                  std::vector<Math::Vector<ScalarType>> uMByComp(
                    vdim, Math::Vector<ScalarType>(TetN));

                  for (size_t comp = 0; comp < vdim; ++comp)
                  {
                    uMByComp[comp].setZero();

                    Math::Vector<ScalarType> tempUR(TetN);
                    for (size_t k = 0; k < TetN; ++k)
                    {
                      const Index sd = scalarLocalDof(cdofs, k);
                      tempUR(static_cast<Index>(k)) =
                        (static_cast<size_t>(sd) < scalarSize)
                        ? data.coeffRef(sd + static_cast<Index>(comp * scalarSize))
                        : ScalarType(0);
                    }

                    Math::Vector<ScalarType> tempUM = sTetChangeScalar * tempUR;

                    for (int k = 0; k < tetInteriorOffset; ++k)
                      uMByComp[comp](k) = tempUM(k);
                  }

                  for (int k = 0; k < numInterior; ++k)
                  {
                    assert(pos < scalarSize);
                    for (size_t comp = 0; comp < vdim; ++comp)
                      uMByComp[comp](tetInteriorOffset + k) = mfemBlock[comp][pos];
                    ++pos;
                  }

                  for (size_t comp = 0; comp < vdim; ++comp)
                    uRElem[comp] = sTetInvChangeScalar * uMByComp[comp];

                  for (size_t k = 0; k < TetN; ++k)
                  {
                    const Index sd = scalarLocalDof(cdofs, k);
                    const size_t s = static_cast<size_t>(sd);
                    if (s >= scalarSize)
                      continue;

                    for (size_t comp = 0; comp < vdim; ++comp)
                      data.coeffRef(sd + static_cast<Index>(comp * scalarSize)) =
                        uRElem[comp](static_cast<Index>(k));
                    written[s] = true;
                  }
                }
                else
                {
                  for (size_t k = 0; k < TetN; ++k)
                  {
                    const Index sd = scalarLocalDof(cdofs, k);
                    const size_t s = static_cast<size_t>(sd);
                    if (s < scalarSize)
                      written[s] = true;
                  }
                }
                break;
              }

              case Geometry::Polytope::Type::Wedge:
              {
                constexpr size_t WedgeN = MFEM::WedgeNodes<K>::Count;
                const int ne = p - 1;
                const int nt = (p - 1) * (p - 2) / 2;
                const int nq = (p - 1) * (p - 1);
                const int wedgeInteriorOffset = 6 + 9 * ne + 2 * nt + 3 * nq;

                assert(scalarLocalSize(cdofs) == WedgeN);

                const int numInterior = static_cast<int>(WedgeN) - wedgeInteriorOffset;
                if (numInterior > 0)
                {
                  std::vector<Math::Vector<ScalarType>> uMByComp(
                    vdim, Math::Vector<ScalarType>(WedgeN));
                  std::vector<Math::Vector<ScalarType>> uRByComp(
                    vdim, Math::Vector<ScalarType>(WedgeN));

                  for (size_t comp = 0; comp < vdim; ++comp)
                  {
                    Math::Vector<ScalarType> tempUR(WedgeN);
                    for (size_t k = 0; k < WedgeN; ++k)
                    {
                      const Index sd = scalarLocalDof(cdofs, k);
                      tempUR(static_cast<Index>(k)) =
                        (static_cast<size_t>(sd) < scalarSize)
                        ? data.coeffRef(sd + static_cast<Index>(comp * scalarSize))
                        : ScalarType(0);
                    }

                    const Math::Vector<ScalarType> tempUM = sWedgeChangeScalar * tempUR;
                    uMByComp[comp].setZero();
                    for (int k = 0; k < wedgeInteriorOffset; ++k)
                      uMByComp[comp](k) = tempUM(k);
                  }

                  for (int k = 0; k < numInterior; ++k)
                  {
                    assert(pos < scalarSize);
                    for (size_t comp = 0; comp < vdim; ++comp)
                      uMByComp[comp](wedgeInteriorOffset + k) = mfemBlock[comp][pos];
                    ++pos;
                  }

                  for (size_t comp = 0; comp < vdim; ++comp)
                    uRByComp[comp] = sWedgeInvChangeScalar * uMByComp[comp];

                  for (size_t k = 0; k < WedgeN; ++k)
                  {
                    const Index sd = scalarLocalDof(cdofs, k);
                    const size_t s = static_cast<size_t>(sd);
                    if (s >= scalarSize)
                      continue;

                    for (size_t comp = 0; comp < vdim; ++comp)
                      data.coeffRef(sd + static_cast<Index>(comp * scalarSize)) =
                        uRByComp[comp](static_cast<Index>(k));
                    written[s] = true;
                  }
                }
                else
                {
                  for (size_t k = 0; k < WedgeN; ++k)
                  {
                    const Index sd = scalarLocalDof(cdofs, k);
                    const size_t s = static_cast<size_t>(sd);
                    if (s < scalarSize)
                      written[s] = true;
                  }
                }
                break;
              }

              default:
              {
                const size_t nScalarLocal = scalarLocalSize(cdofs);
                for (size_t k = 0; k < nScalarLocal; ++k)
                {
                  const Index sdof = scalarLocalDof(cdofs, k);
                  const size_t s = static_cast<size_t>(sdof);
                  if (s >= scalarSize || written[s])
                    continue;

                  setScalarDofFromPos(sdof);
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
            const size_t nScalarLocal = scalarLocalSize(cdofs);
            for (size_t k = 0; k < nScalarLocal; ++k)
            {
              const Index sdof = scalarLocalDof(cdofs, k);
              const size_t s = static_cast<size_t>(sdof);
              if (s >= scalarSize || written[s])
                continue;

              setScalarDofFromPos(sdof);
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
      /// @brief Finite element space type.
      using FESType = Variational::P0<Range, Geometry::Mesh<Context::Local>>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;

      /// @brief Coefficient data storage type.
      using DataType = Math::Vector<ScalarType>;

      using ObjectType = Variational::GridFunction<FESType,
        DataType>; ///< Grid function type being loaded.

      /// @brief Parent class type.
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

        const auto getFec = [&](auto& ctx) { header.fec = _attr(ctx); };
        const auto getVdim = [&](auto& ctx) { header.vdim = _attr(ctx); };
        const auto getOrdering = [&](auto& ctx) {
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
          boost::spirit::x3::string("FiniteElementCollection: ") >> (+char_)[getFec];
        const bool rfec  = boost::spirit::x3::phrase_parse(it, line.end(), pfec, space);
        (void) rfec;
        assert(rfec && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pvdim = boost::spirit::x3::string("VDim:") >> uint_[getVdim];
        const bool rvdim = boost::spirit::x3::phrase_parse(it, line.end(), pvdim, space);
        (void) rvdim;
        assert(rvdim && it == line.end());

        line = MFEM::skipEmptyLinesAndComments(is, m_currentLineNumber);
        it   = line.begin();
        const auto pord = boost::spirit::x3::string("Ordering:") >> uint_[getOrdering];
        const bool rord = boost::spirit::x3::phrase_parse(it, line.end(), pord, space);
        (void) rord;
        assert(rord && it == line.end());

        // -------------------------------------------------------------
        // 2. Read coefficient data
        // -------------------------------------------------------------
        auto& gf  = this->getObject();
        auto& fes = gf.getFiniteElementSpace();
        auto& data = gf.getData();

        const size_t vdim = fes.getVectorDimension();
        assert(header.vdim == vdim);

        const size_t dofCount = fes.getSize();
        assert(dofCount % vdim == 0);
        const size_t scalarDofCount = dofCount / vdim;

        // Resize data vector
        data.resize(dofCount);

        std::vector<ScalarType> tmp;
        tmp.reserve(dofCount);

        ScalarType value;
        while (is >> value)
          tmp.push_back(value);

        if (tmp.size() != dofCount)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Expected " << dofCount << " coefficient values, but found "
            << tmp.size() << "." << Alert::Raise;
        }

        if (header.ordering == MFEM::Ordering::Nodes)
        {
          for (size_t c = 0; c < vdim; ++c)
          {
            for (size_t i = 0; i < scalarDofCount; ++i)
              data(static_cast<Index>(i + c * scalarDofCount)) =
                tmp[i + c * scalarDofCount];
          }
        }
        else
        {
          for (size_t i = 0; i < scalarDofCount; ++i)
          {
            for (size_t c = 0; c < vdim; ++c)
              data(static_cast<Index>(i + c * scalarDofCount)) =
                tmp[c + i * vdim];
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
      /// @brief Range (evaluation value) type.
      using RangeType = Range;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      /// @brief Coefficient data storage type.
      using DataType = Data;

      /// @brief Mesh type.
      using MeshType = Geometry::Mesh<Context>;

      /// @brief Finite element space type.
      using FESType = Variational::P0<Range, MeshType>;

      /// @brief File format handled by this printer base.
      static constexpr FileFormat Format = FileFormat::MFEM;

      /// @brief Grid function type being printed.
      using ObjectType = Variational::GridFunction<FESType, DataType>;

      /// @brief Parent class type.
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
      /// @brief Range (evaluation value) type.
      using RangeType = Range;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      /// @brief Coefficient data storage type.
      using DataType = Data;

      /// @brief Mesh type.
      using MeshType = Geometry::Mesh<Context>;

      /// @brief Finite element space type.
      using FESType = Variational::P1<Range, MeshType>;

      /// @brief File format handled by this printer base.
      static constexpr FileFormat Format = FileFormat::MFEM;

      /// @brief Grid function type being printed.
      using ObjectType = Variational::GridFunction<FESType, DataType>;

      /// @brief Parent class type.
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
      /// @brief Range (evaluation value) type.
      using RangeType = Range;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<RangeType>::ScalarType;

      /// @brief Coefficient data storage type.
      using DataType = Data;

      /// @brief Mesh type.
      using MeshType = Geometry::Mesh<Context>;

      /// @brief Finite element space type.
      using FESType = Variational::H1<K, Range, MeshType>;

      /// @brief File format handled by this printer base.
      static constexpr FileFormat Format = FileFormat::MFEM;

      /// @brief Grid function type being printed.
      using ObjectType = Variational::GridFunction<FESType, DataType>;

      /// @brief Parent class type.
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
      /// @brief Coefficient data storage type.
      using DataType = Math::Vector<Scalar>;

      /// @brief Grid function type being printed.
      using ObjectType = Variational::GridFunction<FES, DataType>;

      /// @brief Parent class type.
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
      /// @brief Finite element space type.
      using FESType    = Variational::H1<K, Range, Geometry::Mesh<Context::Local>>;
      /// @brief Coefficient data storage type.
      using DataType   = Math::Vector<Scalar>;
      /// @brief Grid function type being printed.
      using ObjectType = Variational::GridFunction<FESType, DataType>;
      /// @brief Parent printer base type.
      using Parent     = GridFunctionPrinterBase<FileFormat::MFEM, FESType, DataType>;

      /**
       * @brief Constructs an MFEM H1 grid-function printer.
       * @param[in] gf Grid function to print.
       */
      GridFunctionPrinter(const ObjectType& gf)
        : Parent(gf)
      {}

      /**
       * @brief Prints H1 coefficient data in MFEM node ordering.
       * @param[in,out] os Output stream receiving coefficient data.
       */
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

        const auto toScalarDof = [&](Index dof) -> Index {
          if (vdim > 1 && dof >= scalarSize)
            return dof % scalarSize;
          return dof;
        };

        const auto scalarLocalSize = [&](const auto& dofs) -> size_t {
          assert(static_cast<size_t>(dofs.size()) % vdim == 0);
          return static_cast<size_t>(dofs.size()) / vdim;
        };

        const auto scalarLocalDof = [&](const auto& dofs, size_t local) -> Index {
          if (vdim == 1)
            return toScalarDof(dofs(static_cast<Index>(local)));

          assert(local * vdim < static_cast<size_t>(dofs.size()));
          return toScalarDof(dofs(static_cast<Index>(local * vdim)));
        };

        size_t activeComponent = 0;
        const auto emitScalarDof = [&](Index rodinDof) {
          const Index scalarDof = toScalarDof(rodinDof);
          const size_t s = static_cast<size_t>(scalarDof);
          if (s >= scalarSize)
            return;
          if (written[s])
            return;
          os << data.coeffRef(scalarDof + activeComponent * scalarSize) << '\n';
          written[s] = true;
        };

        //--------------------------------------------------------------------
        // Precomputed change-of-nodes matrices (Rodin Fekete -> MFEM nodes)
        //--------------------------------------------------------------------

        auto& sTriChangeScalar = []() -> const Math::Matrix<Scalar>& {
          static thread_local Math::Matrix<Scalar> C;
          if (C.size() == 0)
          {
            const auto& vMfem = MFEM::VandermondeTriangle<K>::getMatrix();
            const auto& vRodInv = Variational::VandermondeTriangle<K>::getInverse();
            const Math::Matrix<Real> cReal = vMfem * vRodInv;
            C = cReal.template cast<Scalar>();
          }
          return C;
        }();

        auto& sTetChangeScalar = []() -> const Math::Matrix<Scalar>& {
          static thread_local Math::Matrix<Scalar> C;
          if (C.size() == 0)
          {
            const auto& vMfem = MFEM::VandermondeTetrahedron<K>::getMatrix();
            const auto& vRodInv = Variational::VandermondeTetrahedron<K>::getInverse();
            const Math::Matrix<Real> cReal = vMfem * vRodInv;
            C = cReal.template cast<Scalar>();
          }
          return C;
        }();

        auto& sWedgeChangeScalar = []() -> const Math::Matrix<Scalar>& {
          static thread_local Math::Matrix<Scalar> C;
          if (C.size() == 0)
            C = MFEM::WedgeChange<K>::getMatrix().template cast<Scalar>();
          return C;
        }();

        for (activeComponent = 0; activeComponent < vdim; ++activeComponent)
        {
          std::fill(written.begin(), written.end(), false);

          //--------------------------------------------------------------------
          // 1. Vertices
          //--------------------------------------------------------------------

          const size_t nVertices = mesh.getConnectivity().getCount(0);
          std::vector<Index> vertexScalarDof(nVertices);
          for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
          {
            const auto& vdofs = fes.getDOFs(0, v);
            assert(vdofs.size() >= 1 && "H1 vertex should have at least one DOF.");
            vertexScalarDof[v] = toScalarDof(vdofs(0));
          }

          for (Index v = 0; v < static_cast<Index>(nVertices); ++v)
            emitScalarDof(vertexScalarDof[v]);

        //--------------------------------------------------------------------
        // 2. Edges: interior DOFs, oriented vmin -> vmax
        //--------------------------------------------------------------------

        if (D >= 1)
        {
          const auto& conn10  = mesh.getConnectivity().getIncidence(1, 0);
          const auto mfemEdges = MFEM::getMFEMEdgeOrder(mesh);

          std::vector<Index> interior;
          for (Index e : mfemEdges)
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
            const size_t edgeScalarDofs = scalarLocalSize(edofs);

            interior.clear();
            for (size_t k = 0; k < edgeScalarDofs; ++k)
            {
              Index d = scalarLocalDof(edofs, k);
              if (d != vminDof && d != vmaxDof)
                interior.push_back(d);
            }

            if (v0 > v1)
              std::reverse(interior.begin(), interior.end());

            for (Index d : interior)
              emitScalarDof(d);
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
            // Triangle face size / offsets in MFEM ordering
            constexpr size_t TriN = MFEM::TriangleNodes<K>::Count;
            [[maybe_unused]] constexpr size_t QuadN = (K + 1) * (K + 1);
            const int p  = static_cast<int>(K);
            const int nV = 3;
            const int nE = 3 * (p - 1);
            const int triInteriorOffset = nV + nE;

            Math::Vector<Scalar> uRFace(TriN);
            Math::Vector<Scalar> uMFace(TriN);
            const auto mfemFaces = MFEM::getMFEMFaceOrder(mesh);

            for (const auto& mfemFace : mfemFaces)
            {
              const Index f = mfemFace.index;
              const auto faceGeom = mesh.getGeometry(2, f);
              const auto& fdofs   = fes.getDOFs(2, f);

              switch (faceGeom)
              {
                case Geometry::Polytope::Type::Triangle:
                {
                  assert(scalarLocalSize(fdofs) == TriN);
                  const auto& faceVertices = mesh.getConnectivity().getPolytope(2, f);
                  const auto faceNodeMap = MFEM::getTriangleNodeMap<K>(faceVertices, mfemFace.vertices);

                  for (size_t k = 0; k < TriN; ++k)
                  {
                    const Index d = scalarLocalDof(fdofs, faceNodeMap[k]);
                    uRFace(static_cast<Index>(k)) =
                      data.coeffRef(d + static_cast<Index>(activeComponent * scalarSize));
                  }
                  uMFace = sTriChangeScalar * uRFace;

                  int loc = 0;
                  for (int j = 1; j < p; ++j)
                  {
                    for (int i = 1; i + j < p; ++i)
                    {
                      const int idx = triInteriorOffset + loc++;
                      os << uMFace(static_cast<Index>(idx)) << '\n';
                    }
                  }
                  break;
                }

                case Geometry::Polytope::Type::Quadrilateral:
                {
                  assert(scalarLocalSize(fdofs) == QuadN);
                  const auto& faceVertices = mesh.getConnectivity().getPolytope(2, f);
                  const auto faceNodeMap = MFEM::getQuadrilateralNodeMap<K>(faceVertices, mfemFace.vertices);

                  for (size_t j = 1; j < K; ++j)
                  {
                    for (size_t i = 1; i < K; ++i)
                    {
                      const size_t targetLocal = MFEM::getQuadrilateralNodeIndex<K>(i, j);
                      const Index sdof = scalarLocalDof(fdofs, faceNodeMap[targetLocal]);
                      const size_t s = static_cast<size_t>(sdof);
                      if (s < scalarSize && !written[s])
                      {
                        os << data.coeffRef(sdof + static_cast<Index>(activeComponent * scalarSize)) << '\n';
                        written[s] = true;
                      }
                    }
                  }
                  break;
                }

                default:
                {
                  const size_t nScalarLocal = scalarLocalSize(fdofs);
                  for (size_t k = 0; k < nScalarLocal; ++k)
                    emitScalarDof(scalarLocalDof(fdofs, k));
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
              const size_t nScalarLocal = scalarLocalSize(fdofs);
              for (size_t k = 0; k < nScalarLocal; ++k)
                emitScalarDof(scalarLocalDof(fdofs, k));
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

          Math::Vector<Scalar> uRElem(TriN);
          Math::Vector<Scalar> uMElem(TriN);

          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto geom = mesh.getGeometry(2, c);
            const auto& cdofs = fes.getDOFs(2, c);

            switch (geom)
            {
              case Geometry::Polytope::Type::Triangle:
              {
                assert(scalarLocalSize(cdofs) == TriN);

                for (size_t k = 0; k < TriN; ++k)
                {
                  const Index d = scalarLocalDof(cdofs, k);
                  uRElem(static_cast<Index>(k)) =
                    data.coeffRef(d + static_cast<Index>(activeComponent * scalarSize));
                }
                uMElem = sTriChangeScalar * uRElem;

                int loc = 0;
                for (int j = 1; j < p; ++j)
                {
                  for (int i = 1; i + j < p; ++i)
                  {
                    const int idx = triInteriorOffset + loc++;
                    os << uMElem(static_cast<Index>(idx)) << '\n';
                  }
                }
                break;
              }

              default:
              {
                const size_t nScalarLocal = scalarLocalSize(cdofs);
                for (size_t k = 0; k < nScalarLocal; ++k)
                  emitScalarDof(scalarLocalDof(cdofs, k));
                break;
              }
            }
          }
        }
        else if (D == 3)
        {
          constexpr size_t TetN = MFEM::TetrahedronNodes<K>::Count;
          const int p = static_cast<int>(K);

          const int nV    = 4;
          const int nE    = 6 * (p - 1);
          const int nF    = 2 * (p - 1) * (p - 2); // 4 faces × (p-1)(p-2)/2
          const int tetInteriorOffset = nV + nE + nF;

          Math::Vector<Scalar> uRElem(TetN);
          Math::Vector<Scalar> uMElem(TetN);

          for (Index c = 0; c < static_cast<Index>(nCells); ++c)
          {
            const auto  cellGeom = mesh.getGeometry(3, c);
            const auto& cdofs    = fes.getDOFs(3, c);

            switch (cellGeom)
            {
              //---------------------------------------------------------------
              // Tetrahedron: change-of-nodes for interior DOFs
              //---------------------------------------------------------------
              case Geometry::Polytope::Type::Tetrahedron:
              {
                assert(scalarLocalSize(cdofs) == TetN);

                for (size_t k = 0; k < TetN; ++k)
                {
                  const Index d = scalarLocalDof(cdofs, k);
                  uRElem(static_cast<Index>(k)) =
                    data.coeffRef(d + static_cast<Index>(activeComponent * scalarSize));
                }
                uMElem = sTetChangeScalar * uRElem;

                for (int idx = tetInteriorOffset; idx < static_cast<int>(TetN); ++idx)
                  os << uMElem(static_cast<Index>(idx)) << '\n';
                break;
              }

              //---------------------------------------------------------------
              // Wedge / Prism: change-of-nodes for interior DOFs
              //---------------------------------------------------------------
              case Geometry::Polytope::Type::Wedge:
              {
                constexpr size_t WedgeN = MFEM::WedgeNodes<K>::Count;
                const int ne = p - 1;
                const int nt = (p - 1) * (p - 2) / 2;
                const int nq = (p - 1) * (p - 1);
                const int wedgeInteriorOffset = 6 + 9 * ne + 2 * nt + 3 * nq;

                assert(scalarLocalSize(cdofs) == WedgeN);

                Math::Vector<Scalar> uRWedge(WedgeN);
                Math::Vector<Scalar> uMWedge(WedgeN);

                for (size_t k = 0; k < WedgeN; ++k)
                {
                  const Index d = scalarLocalDof(cdofs, k);
                  uRWedge(static_cast<Index>(k)) =
                    data.coeffRef(d + static_cast<Index>(activeComponent * scalarSize));
                }
                uMWedge = sWedgeChangeScalar * uRWedge;

                for (int idx = wedgeInteriorOffset; idx < static_cast<int>(WedgeN); ++idx)
                  os << uMWedge(static_cast<Index>(idx)) << '\n';
                break;
              }

              //---------------------------------------------------------------
              // Other 3D cell types (including Hexahedron): DOF-based interiors
              //---------------------------------------------------------------
              default:
              {
                const size_t nScalarLocal = scalarLocalSize(cdofs);
                for (size_t k = 0; k < nScalarLocal; ++k)
                  emitScalarDof(scalarLocalDof(cdofs, k));
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
            const size_t nScalarLocal = scalarLocalSize(cdofs);
            for (size_t k = 0; k < nScalarLocal; ++k)
              emitScalarDof(scalarLocalDof(cdofs, k));
          }
        }
        }
      }
  };
}

#endif
