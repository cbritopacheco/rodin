/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <fstream>
#include <boost/filesystem/fstream.hpp>
#include <gtest/gtest.h>

#include <Rodin/IO.h>
#include <Rodin/IO/MFEM.h>
#include <Rodin/IO/MEDIT.h>
#include <Rodin/IO/GMSH.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin::IO;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_IO_MeshLoader, SanityTest_MEDIT_2D_Square)
  {
    static constexpr const char* filename = "mmg/Square.medit.mesh";
    static constexpr size_t sdim = 2;
    boost::filesystem::path meshfile;
    meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
    meshfile.append(filename);

    boost::filesystem::ifstream in(meshfile);

    Mesh mesh;
    MeshLoader<FileFormat::MEDIT, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    EXPECT_EQ(mesh.getSpaceDimension(), sdim);

    size_t d = 0;
    EXPECT_EQ(mesh.getPolytopeCount(d), 4);
    EXPECT_EQ(mesh.getPolytopeCount(d), mesh.getVertexCount());
    EXPECT_EQ(mesh.getAttribute(d, 0), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
    EXPECT_EQ(mesh.getAttribute(d, 1), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
    EXPECT_EQ(mesh.getAttribute(d, 2), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
    EXPECT_EQ(mesh.getAttribute(d, 3), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);

    d = 1;
    EXPECT_EQ(mesh.getPolytopeCount(d), 5);
    EXPECT_EQ(mesh.getAttribute(d, 0), 1);
    EXPECT_EQ(mesh.getAttribute(d, 1), 2);

    d = 2;
    EXPECT_EQ(mesh.getPolytopeCount(d), 2);
    EXPECT_EQ(mesh.getPolytopeCount(d), mesh.getCellCount());
    EXPECT_EQ(mesh.getAttribute(d, 0), 1);
    EXPECT_EQ(mesh.getAttribute(d, 1), 2);
  }

  TEST(Rodin_IO_MeshLoader, SanityTest_MFEM_2D_Square)
  {
    static constexpr const char* filename = "mfem/Square.mfem.mesh";
    boost::filesystem::path meshfile;
    meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
    meshfile.append(filename);

    boost::filesystem::ifstream in(meshfile);

    Mesh mesh;
    MeshLoader<FileFormat::MFEM, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    size_t d = 0;
    EXPECT_EQ(mesh.getPolytopeCount(d), 4);
    EXPECT_EQ(mesh.getPolytopeCount(d), mesh.getVertexCount());

    // MFEM format does not support vertex attributes
    EXPECT_EQ(mesh.getAttribute(d, 0), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
    EXPECT_EQ(mesh.getAttribute(d, 1), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
    EXPECT_EQ(mesh.getAttribute(d, 2), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
    EXPECT_EQ(mesh.getAttribute(d, 3), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);

    d = 1;
    EXPECT_EQ(mesh.getPolytopeCount(d), 5);
    EXPECT_EQ(mesh.getAttribute(d, 0), 1);
    EXPECT_EQ(mesh.getAttribute(d, 1), 2);

    d = 2;
    EXPECT_EQ(mesh.getPolytopeCount(d), 2);
    EXPECT_EQ(mesh.getPolytopeCount(d), mesh.getCellCount());
    EXPECT_EQ(mesh.getAttribute(d, 0), 1);
    EXPECT_EQ(mesh.getAttribute(d, 1), 2);
  }

  TEST(Rodin_IO_MeshLoader, SanityTest_MFEM_2D_StarSquare)
  {
    static constexpr const char* filename = "mfem/StarSquare.mfem.mesh";
    static constexpr size_t vcount = 101;
    static constexpr size_t ecount = 80;
    static constexpr size_t sdim = 2;

    boost::filesystem::path meshfile;
    meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
    meshfile.append(filename);

    boost::filesystem::ifstream in(meshfile);

    Mesh mesh;
    MeshLoader<FileFormat::MFEM, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    EXPECT_EQ(mesh.getSpaceDimension(), sdim);

    size_t d = 0;
    EXPECT_EQ(mesh.getPolytopeCount(d), vcount);
    EXPECT_EQ(mesh.getPolytopeCount(d), mesh.getVertexCount());

    for (size_t i = 0; i < vcount; i++)
      EXPECT_EQ(mesh.getAttribute(d, i), RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);

    d = 2;
    for (size_t i = 0; i < ecount; i++)
      EXPECT_EQ(mesh.getAttribute(d, i), 1);

    for (size_t i = 0; i < ecount; i++)
      EXPECT_EQ(mesh.getGeometry(d, i), Polytope::Type::Quadrilateral);
  }

  TEST(Rodin_IO_MeshLoader, SanityTest_GMSH_2D_Triangle)
  {
    // Create a simple GMSH format string in memory
    std::stringstream gmshContent;
    gmshContent << "$MeshFormat\n";
    gmshContent << "2.2 0 8\n";
    gmshContent << "$EndMeshFormat\n";
    gmshContent << "$Nodes\n";
    gmshContent << "3\n";
    gmshContent << "1 0.0 0.0 0.0\n";
    gmshContent << "2 1.0 0.0 0.0\n";
    gmshContent << "3 0.0 1.0 0.0\n";
    gmshContent << "$EndNodes\n";
    gmshContent << "$Elements\n";
    gmshContent << "1\n";
    gmshContent << "1 2 2 1 1 1 2 3\n";  // Triangle element
    gmshContent << "$EndElements\n";

    Mesh mesh;
    MeshLoader<FileFormat::GMSH, Rodin::Context::Local> loader(mesh);
    loader.load(gmshContent);

    // Verify mesh structure
    size_t d = 0;  // Vertices
    EXPECT_EQ(mesh.getVertexCount(), 3);
    EXPECT_EQ(mesh.getPolytopeCount(d), 3);

    d = 2;  // Cells (triangles)
    EXPECT_EQ(mesh.getCellCount(), 1);
    EXPECT_EQ(mesh.getPolytopeCount(d), 1);
    EXPECT_EQ(mesh.getGeometry(d, 0), Polytope::Type::Triangle);
    EXPECT_EQ(mesh.getAttribute(d, 0), 1);
  }

  TEST(Rodin_IO_MeshLoader, SanityTest_GMSH_2D_Quadrilateral)
  {
    // Create a simple GMSH format string for a quad
    std::stringstream gmshContent;
    gmshContent << "$MeshFormat\n";
    gmshContent << "2.2 0 8\n";
    gmshContent << "$EndMeshFormat\n";
    gmshContent << "$Nodes\n";
    gmshContent << "4\n";
    gmshContent << "1 0.0 0.0 0.0\n";
    gmshContent << "2 1.0 0.0 0.0\n";
    gmshContent << "3 1.0 1.0 0.0\n";
    gmshContent << "4 0.0 1.0 0.0\n";
    gmshContent << "$EndNodes\n";
    gmshContent << "$Elements\n";
    gmshContent << "1\n";
    gmshContent << "1 3 2 2 2 1 2 3 4\n";  // Quadrilateral element with attribute 2
    gmshContent << "$EndElements\n";

    Mesh mesh;
    MeshLoader<FileFormat::GMSH, Rodin::Context::Local> loader(mesh);
    loader.load(gmshContent);

    // Verify mesh structure
    size_t d = 0;  // Vertices
    EXPECT_EQ(mesh.getVertexCount(), 4);
    EXPECT_EQ(mesh.getPolytopeCount(d), 4);

    d = 2;  // Cells (quads)
    EXPECT_EQ(mesh.getCellCount(), 1);
    EXPECT_EQ(mesh.getPolytopeCount(d), 1);
    EXPECT_EQ(mesh.getGeometry(d, 0), Polytope::Type::Quadrilateral);
    EXPECT_EQ(mesh.getAttribute(d, 0), 2);
  }

  TEST(Rodin_IO_MeshPrinter, SanityTest_GMSH_2D_Triangle)
  {
    // Create a mesh programmatically
    Mesh mesh;
    auto builder = mesh.newBuilder();
    builder.initialize(2, 3);  // 2D space, 3 vertices
    
    // Add vertices
    builder.vertex(0, 0.0, 0.0, 0.0);
    builder.vertex(1, 1.0, 0.0, 0.0);
    builder.vertex(2, 0.0, 1.0, 0.0);
    
    // Add triangle
    builder.polytope(Polytope::Type::Triangle, {0, 1, 2}, 1);
    
    mesh = builder.finalize();

    // Print to GMSH format
    std::stringstream output;
    MeshPrinter<FileFormat::GMSH, Rodin::Context::Local> printer(mesh);
    printer.print(output);

    // Verify the output contains expected sections
    std::string result = output.str();
    EXPECT_TRUE(result.find("$MeshFormat") != std::string::npos);
    EXPECT_TRUE(result.find("$Nodes") != std::string::npos);
    EXPECT_TRUE(result.find("$Elements") != std::string::npos);
    EXPECT_TRUE(result.find("3") != std::string::npos);  // Number of nodes
    EXPECT_TRUE(result.find("1") != std::string::npos);  // Number of elements
  }
}
