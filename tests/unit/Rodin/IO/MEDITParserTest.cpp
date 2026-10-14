/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <array>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include "Rodin/IO/MEDIT.h"

using namespace Rodin;
using namespace Rodin::IO;

/// @brief Every declared MEDIT keyword has an exact, reversible spelling.
TEST(MEDITParserContracts, KeywordSpellingsAndUnknownNames)
{
  using K = MEDIT::Keyword;
  const std::array<std::pair<K, const char*>, 28> keywords{
    {{K::MeshVersionFormatted, "MeshVersionFormatted"}, {K::Dimension, "Dimension"},
      {K::Vertices, "Vertices"}, {K::Triangles, "Triangles"},
      {K::Quadrilaterals, "Quadrilaterals"}, {K::Tetrahedra, "Tetrahedra"},
      {K::Hexahedra, "Hexahedra"}, {K::Wedges, "Wedges"}, {K::Pyramids, "Pyramids"},
      {K::Corners, "Corners"}, {K::Ridges, "Ridges"}, {K::Edges, "Edges"},
      {K::SolAtVertices, "SolAtVertices"}, {K::SolAtEdges, "SolAtEdges"},
      {K::SolAtTriangles, "SolAtTriangles"},
      {K::SolAtQuadrilaterals, "SolAtQuadrilaterals"},
      {K::SolAtTetrahedra, "SolAtTetrahedra"}, {K::SolAtPentahedra, "SolAtPentahedra"},
      {K::SolAtHexahedra, "SolAtHexahedra"}, {K::RequiredVertices, "RequiredVertices"},
      {K::RequiredEdges, "RequiredEdges"}, {K::RequiredTriangles, "RequiredTriangles"},
      {K::RequiredTetrahedra, "RequiredTetrahedra"}, {K::Normals, "Normals"},
      {K::NormalAtVertices, "NormalAtVertices"}, {K::Tangents, "Tangents"},
      {K::TangentAtVertices, "TangentAtVertices"}, {K::End, "End"}}};
  for (const auto& [keyword, spelling] : keywords)
  {
    SCOPED_TRACE(spelling);
    EXPECT_STREQ(MEDIT::toCharString(keyword), spelling);
    ASSERT_TRUE(MEDIT::toKeyword(spelling));
    EXPECT_EQ(*MEDIT::toKeyword(spelling), keyword);
    EXPECT_TRUE(std::string(spelling) == keyword);
    EXPECT_FALSE(std::string(spelling) != keyword);
    EXPECT_TRUE((std::string(spelling) + "!") != keyword);
    std::ostringstream stream;
    stream << keyword;
    EXPECT_TRUE(stream.good());
    EXPECT_EQ(stream.str(), spelling);
  }
  for (const char* spelling : {"", "triangles", "Triangle", "Dimension2", "unknown"})
    EXPECT_FALSE(MEDIT::toKeyword(spelling));
  EXPECT_EQ(MEDIT::toCharString(static_cast<K>(-1)), nullptr);
}

template <class PointerRange>
class MEDITIteratorContracts : public ::testing::Test
{
  protected:
    /// @brief Evaluates the parser over either raw characters or string iterators.
    template <class Parser>
    auto parse(const Parser& parser, const std::string& text) const
    {
      if constexpr (PointerRange::value)
        return parser(text.data(), text.data() + text.size());
      else
        return parser(text.cbegin(), text.cend());
    }
};

using MEDITIteratorTypes = ::testing::Types<std::true_type, std::false_type>;
TYPED_TEST_SUITE(MEDITIteratorContracts, MEDITIteratorTypes);

/// @brief Whitespace and alphabetic keywords require complete input consumption.
TYPED_TEST(MEDITIteratorContracts, EmptyLinesAndKeywords)
{
  for (const char* text : {"", " ", "\t  "})
    EXPECT_TRUE(this->parse(MEDIT::ParseEmptyLine{}, text));
  for (const char* text : {"x", " 42", " Vertices"})
    EXPECT_FALSE(this->parse(MEDIT::ParseEmptyLine{}, text));
  const auto keyword = this->parse(MEDIT::ParseKeyword{}, " \tTriangles  ");
  ASSERT_TRUE(keyword);
  EXPECT_EQ(*keyword, "Triangles");
  for (const char* text : {"", " ", "42", "Triangles7", "Triangles tail", "End#"})
  {
    SCOPED_TRACE(text);
    EXPECT_FALSE(this->parse(MEDIT::ParseKeyword{}, text));
  }
}

/// @brief Integer limits are accepted exactly; overflow and trailing tokens are rejected.
TYPED_TEST(MEDITIteratorContracts, SignedAndUnsignedIntegerLimits)
{
  for (int value :
    {std::numeric_limits<int>::min(), -7, 0, 7, std::numeric_limits<int>::max()})
  {
    const auto parsed = this->parse(MEDIT::ParseInteger{}, " " + std::to_string(value));
    ASSERT_TRUE(parsed);
    EXPECT_EQ(*parsed, value);
  }
  for (unsigned int value : {0u, 7u, std::numeric_limits<unsigned int>::max()})
  {
    const auto parsed =
      this->parse(MEDIT::ParseUnsignedInteger{}, std::to_string(value) + "\t");
    ASSERT_TRUE(parsed);
    EXPECT_EQ(*parsed, value);
  }
  for (const char* text : {"", " ", "x", "2 3", "12x", "1.0", "9999999999999999999"})
  {
    EXPECT_FALSE(this->parse(MEDIT::ParseInteger{}, text));
    EXPECT_FALSE(this->parse(MEDIT::ParseUnsignedInteger{}, text));
  }
  EXPECT_FALSE(this->parse(MEDIT::ParseUnsignedInteger{}, "-1"));
}

/// @brief Structured declarations distinguish missing values, unexpected keywords and malformed keyword tokens.
TYPED_TEST(MEDITIteratorContracts, StructuredHeaderDeclarations)
{
  const auto version =
    this->parse(MEDIT::ParseMeshVersionFormatted{}, " MeshVersionFormatted 2 ");
  ASSERT_TRUE(version);
  EXPECT_EQ(*version, 2u);
  for (unsigned int dimension : {0u, 1u, 2u, 3u})
  {
    const auto parsed =
      this->parse(MEDIT::ParseDimension{}, "Dimension " + std::to_string(dimension));
    ASSERT_TRUE(parsed);
    EXPECT_EQ(*parsed, dimension);
  }
  for (const char* suffix : {"", " -1", " 1 2", " 9999999999999999999"})
  {
    EXPECT_FALSE(this->parse(MEDIT::ParseDimension{}, std::string("Dimension") + suffix));
    EXPECT_FALSE(this->parse(
      MEDIT::ParseMeshVersionFormatted{}, std::string("MeshVersionFormatted") + suffix));
  }
  EXPECT_THROW(this->parse(MEDIT::ParseDimension{}, "MeshVersionFormatted 2"),
    MEDIT::UnexpectedKeywordException);
  EXPECT_THROW(this->parse(MEDIT::ParseMeshVersionFormatted{}, "Dimension 2"),
    MEDIT::UnexpectedKeywordException);
  EXPECT_THROW(this->parse(MEDIT::ParseDimension{}, "Dim ension 3"),
    MEDIT::UnexpectedKeywordException);
  EXPECT_THROW(this->parse(MEDIT::ParseMeshVersionFormatted{}, "MeshVersion Formatted 2"),
    MEDIT::UnexpectedKeywordException);
  for (const char* text : {"", " ", "123 2"})
  {
    EXPECT_THROW(this->parse(MEDIT::ParseDimension{}, text), Alert::Exception);
    EXPECT_THROW(this->parse(MEDIT::ParseMeshVersionFormatted{}, text), Alert::Exception);
  }
}

/// @brief Entity and zero-through-three-dimensional vertex parsers preserve exact coordinates, indices and attributes.
TYPED_TEST(MEDITIteratorContracts, VertexAndConnectivityRecords)
{
  const std::array<std::string, 4> coordinates{
    "", "1.25 ", "1.25 -2.5 ", "1.25 -2.5 3.75 "};
  const std::array<Real, 3> expectedCoordinates{1.25, -2.5, 3.75};
  for (size_t dimension = 0; dimension <= 3; ++dimension)
  {
    SCOPED_TRACE(dimension);
    const MEDIT::ParseVertex parser(dimension);
    const auto vertex = this->parse(parser, coordinates[dimension] + "19");
    ASSERT_TRUE(vertex);
    EXPECT_EQ(vertex->vertex.size(), dimension);
    EXPECT_EQ(vertex->attribute, Geometry::Attribute(19));
    for (size_t axis = 0; axis < dimension; ++axis)
      EXPECT_EQ(vertex->vertex[axis], expectedCoordinates[axis]);
    EXPECT_FALSE(this->parse(parser, coordinates[dimension]));
    EXPECT_FALSE(this->parse(parser, coordinates[dimension] + "19 extra"));
  }
  for (size_t count = 1; count <= 8; ++count)
  {
    SCOPED_TRACE(count);
    std::string indices;
    for (size_t vertex = 1; vertex <= count; ++vertex)
      indices += std::to_string(vertex) + ' ';
    const MEDIT::ParseEntity parser(count);
    const auto entity = this->parse(parser, indices + "23");
    ASSERT_TRUE(entity);
    EXPECT_EQ(entity->vertices.size(), count);
    EXPECT_EQ(entity->attribute, Geometry::Attribute(23));
    for (size_t vertex = 0; vertex < count; ++vertex)
      EXPECT_EQ(entity->vertices(vertex), vertex + 1);
    EXPECT_FALSE(this->parse(parser, ""));
    EXPECT_FALSE(this->parse(parser, indices));
    EXPECT_FALSE(this->parse(parser, indices + "23 extra"));
  }
}
