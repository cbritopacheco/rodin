/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <sstream>
#include <fstream>

#include "Rodin/Geometry.h"
#include "Rodin/IO/EnSight6.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::IO;

class EnSight6Test : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(EnSight6Test, MeshPrinterBasic)
{
  // Create a simple triangle mesh
  Mesh<Context::Local> mesh;
  auto& builder = mesh.getBuilder();
  
  // Add vertices
  Math::Vector<Real> v1(2); v1 << 0.0, 0.0;
  Math::Vector<Real> v2(2); v2 << 1.0, 0.0;
  Math::Vector<Real> v3(2); v3 << 0.5, 1.0;
  
  builder.reserve(0, 3); // Reserve space for 3 vertices
  builder.vertex(v1);
  builder.vertex(v2);
  builder.vertex(v3);
  
  // Add triangle
  Array<Index> vertices{0, 1, 2};
  builder.reserve(2, 1); // Reserve space for 1 triangle (2D element)
  builder.polytope(Polytope::Type::Triangle, std::move(vertices));
  builder.attribute({2, 0}, 1); // Set attribute for the triangle
  
  mesh.finalize();

  // Test mesh printing
  std::ostringstream oss;
  MeshPrinter<FileFormat::ENSIGHT6, Context::Local> printer(mesh);
  printer.print(oss);
  
  std::string output = oss.str();
  
  // Basic checks for EnSight6 format
  EXPECT_TRUE(output.find("EnSight6 Geometry File Format") != std::string::npos);
  EXPECT_TRUE(output.find("coordinates") != std::string::npos);
  EXPECT_TRUE(output.find("part") != std::string::npos);
  EXPECT_TRUE(output.find("tria3") != std::string::npos);
}

TEST_F(EnSight6Test, MeshLoaderBasic)
{
  // Create EnSight6 geometry data
  std::string ensightData = R"(EnSight6 Geometry File Format
Rodin v0.0.1
node id given
element id given
coordinates
3
0  0.0000e+00  0.0000e+00  0.0000e+00
1  1.0000e+00  0.0000e+00  0.0000e+00
2  5.0000e-01  1.0000e+00  0.0000e+00
part 1
Triangle_Part
tria3
1
0 0 1 2
)";

  std::istringstream iss(ensightData);
  
  Mesh<Context::Local> mesh;
  MeshLoader<FileFormat::ENSIGHT6, Context::Local> loader(mesh);
  
  // Test that loading doesn't crash
  EXPECT_NO_THROW(loader.load(iss));
  
  // Basic checks
  EXPECT_EQ(mesh.getVertexCount(), 3);
  EXPECT_EQ(mesh.getCellCount(), 1);
}

TEST_F(EnSight6Test, ParserTests)
{
  // Test keyword parsing
  {
    std::string text = "coordinates";
    auto result = EnSight6::ParseKeyword()(text.begin(), text.end());
    EXPECT_TRUE(result);
    EXPECT_EQ(*result, "coordinates");
  }
  
  // Test unsigned integer parsing
  {
    std::string text = "123";
    auto result = EnSight6::ParseUnsignedInteger()(text.begin(), text.end());
    EXPECT_TRUE(result);
    EXPECT_EQ(*result, 123u);
  }
  
  // Test real number parsing
  {
    std::string text = "1.23e+02";
    auto result = EnSight6::ParseReal()(text.begin(), text.end());
    EXPECT_TRUE(result);
    EXPECT_DOUBLE_EQ(*result, 123.0);
  }
  
  // Test empty line parsing
  {
    std::string emptyLine = "   ";
    EXPECT_TRUE(EnSight6::ParseEmptyLine()(emptyLine.begin(), emptyLine.end()));
    
    std::string commentLine = "# This is a comment";
    EXPECT_TRUE(EnSight6::ParseEmptyLine()(commentLine.begin(), commentLine.end()));
    
    std::string contentLine = "coordinates";
    EXPECT_FALSE(EnSight6::ParseEmptyLine()(contentLine.begin(), contentLine.end()));
  }
}

TEST_F(EnSight6Test, ElementTypeConversion)
{
  // Test ElementType to PolytopeType conversion
  EXPECT_EQ(EnSight6::toPolytopeType(EnSight6::ElementType::point), Polytope::Type::Point);
  EXPECT_EQ(EnSight6::toPolytopeType(EnSight6::ElementType::bar2), Polytope::Type::Segment);
  EXPECT_EQ(EnSight6::toPolytopeType(EnSight6::ElementType::tria3), Polytope::Type::Triangle);
  EXPECT_EQ(EnSight6::toPolytopeType(EnSight6::ElementType::quad4), Polytope::Type::Quadrilateral);
  EXPECT_EQ(EnSight6::toPolytopeType(EnSight6::ElementType::tetra4), Polytope::Type::Tetrahedron);
  EXPECT_EQ(EnSight6::toPolytopeType(EnSight6::ElementType::penta6), Polytope::Type::Wedge);
  
  // Test string to ElementType conversion
  EXPECT_EQ(EnSight6::toElementType("point"), EnSight6::ElementType::point);
  EXPECT_EQ(EnSight6::toElementType("bar2"), EnSight6::ElementType::bar2);
  EXPECT_EQ(EnSight6::toElementType("tria3"), EnSight6::ElementType::tria3);
  EXPECT_EQ(EnSight6::toElementType("quad4"), EnSight6::ElementType::quad4);
  EXPECT_EQ(EnSight6::toElementType("tetra4"), EnSight6::ElementType::tetra4);
  EXPECT_EQ(EnSight6::toElementType("penta6"), EnSight6::ElementType::penta6);
  
  // Test invalid conversion
  EXPECT_FALSE(EnSight6::toElementType("invalid"));
}

TEST_F(EnSight6Test, RoundTripTest)
{
  // Create a simple mesh
  Mesh<Context::Local> originalMesh;
  auto& builder = originalMesh.getBuilder();
  
  // Add vertices for a tetrahedron
  Math::Vector<Real> v1(3); v1 << 0.0, 0.0, 0.0;
  Math::Vector<Real> v2(3); v2 << 1.0, 0.0, 0.0;
  Math::Vector<Real> v3(3); v3 << 0.5, 1.0, 0.0;
  Math::Vector<Real> v4(3); v4 << 0.5, 0.5, 1.0;
  
  builder.reserve(0, 4); // Reserve space for 4 vertices
  builder.vertex(v1);
  builder.vertex(v2);
  builder.vertex(v3);
  builder.vertex(v4);
  
  // Add tetrahedron
  Array<Index> vertices{0, 1, 2, 3};
  builder.reserve(3, 1); // Reserve space for 1 tetrahedron (3D element)
  builder.polytope(Polytope::Type::Tetrahedron, std::move(vertices));
  builder.attribute({3, 0}, 1); // Set attribute for the tetrahedron
  
  originalMesh.finalize();

  // Export to EnSight6
  std::ostringstream oss;
  MeshPrinter<FileFormat::ENSIGHT6, Context::Local> printer(originalMesh);
  printer.print(oss);
  
  // Import back from EnSight6
  std::istringstream iss(oss.str());
  Mesh<Context::Local> loadedMesh;
  MeshLoader<FileFormat::ENSIGHT6, Context::Local> loader(loadedMesh);
  
  EXPECT_NO_THROW(loader.load(iss));
  
  // Verify basic properties match
  EXPECT_EQ(loadedMesh.getVertexCount(), originalMesh.getVertexCount());
  EXPECT_EQ(loadedMesh.getCellCount(), originalMesh.getCellCount());
}