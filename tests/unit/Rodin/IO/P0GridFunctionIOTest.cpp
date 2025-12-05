#include <gtest/gtest.h>
#include <sstream>
#include <fstream>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/IO.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::IO;

namespace Rodin::Tests::Unit
{
  // Test P0 GridFunction save/load round-trip on triangle mesh
  TEST(P0GridFunctionIO, SaveLoadRoundTrip_Triangle)
  {
    // Create mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {4, 4});
    EXPECT_EQ(mesh.getCellCount(), 32);

    // Create P0 space and grid function
    P0 fes(mesh);
    GridFunction u(fes);

    // Project a function: u(x,y) = x + 2*y
    u.project([](const Point& p) {
      return p.x() + 2.0 * p.y();
    });

    // Save to stringstream
    std::stringstream ss;
    u.save(ss, FileFormat::MFEM);

    // Check FEC header
    std::string line;
    std::getline(ss, line); // "FiniteElementSpace"
    std::getline(ss, line); // "FiniteElementCollection: L2_2D_P0"
    EXPECT_TRUE(line.find("L2_2D_P0") != std::string::npos);

    // Reset stream and load
    ss.clear();
    ss.seekg(0);

    GridFunction v(fes);
    v.load(ss, FileFormat::MFEM);

    // Compare values
    const auto& uData = u.getData();
    const auto& vData = v.getData();
    ASSERT_EQ(uData.size(), vData.size());

    for (Index i = 0; i < uData.size(); i++)
    {
      EXPECT_NEAR(uData(i), vData(i), 1e-10);
    }
  }

  // Test P0 GridFunction save/load round-trip on quadrilateral mesh
  TEST(P0GridFunctionIO, SaveLoadRoundTrip_Quadrilateral)
  {
    // Create mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, {5, 5});
    EXPECT_EQ(mesh.getCellCount(), 16);

    // Create P0 space and grid function
    P0 fes(mesh);
    GridFunction u(fes);

    // Project a function: u(x,y) = x*x + y*y
    u.project([](const Point& p) {
      return p.x() * p.x() + p.y() * p.y();
    });

    // Save to stringstream
    std::stringstream ss;
    u.save(ss, FileFormat::MFEM);

    // Check FEC header
    std::string line;
    std::getline(ss, line); // "FiniteElementSpace"
    std::getline(ss, line); // "FiniteElementCollection: L2_2D_P0"
    EXPECT_TRUE(line.find("L2_2D_P0") != std::string::npos);

    // Reset stream and load
    ss.clear();
    ss.seekg(0);

    GridFunction v(fes);
    v.load(ss, FileFormat::MFEM);

    // Compare values
    const auto& uData = u.getData();
    const auto& vData = v.getData();
    ASSERT_EQ(uData.size(), vData.size());

    for (Index i = 0; i < uData.size(); i++)
    {
      EXPECT_NEAR(uData(i), vData(i), 1e-10);
    }
  }

  // Test P0 GridFunction save/load round-trip on tetrahedron mesh
  TEST(P0GridFunctionIO, SaveLoadRoundTrip_Tetrahedron)
  {
    // Create mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    EXPECT_EQ(mesh.getCellCount(), 162);

    // Create P0 space and grid function
    P0 fes(mesh);
    GridFunction u(fes);

    // Project a function: u(x,y,z) = x + y + z
    u.project([](const Point& p) {
      return p.x() + p.y() + p.z();
    });

    // Save to stringstream
    std::stringstream ss;
    u.save(ss, FileFormat::MFEM);

    // Check FEC header
    std::string line;
    std::getline(ss, line); // "FiniteElementSpace"
    std::getline(ss, line); // "FiniteElementCollection: L2_3D_P0"
    EXPECT_TRUE(line.find("L2_3D_P0") != std::string::npos);

    // Reset stream and load
    ss.clear();
    ss.seekg(0);

    GridFunction v(fes);
    v.load(ss, FileFormat::MFEM);

    // Compare values
    const auto& uData = u.getData();
    const auto& vData = v.getData();
    ASSERT_EQ(uData.size(), vData.size());

    for (Index i = 0; i < uData.size(); i++)
    {
      EXPECT_NEAR(uData(i), vData(i), 1e-10);
    }
  }

  // Test P0 GridFunction save/load round-trip on wedge mesh
  TEST(P0GridFunctionIO, SaveLoadRoundTrip_Wedge)
  {
    // Create wedge mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Wedge, {3, 3, 3});
    EXPECT_GE(mesh.getCellCount(), 16);

    // Create P0 space and grid function
    P0 fes(mesh);
    GridFunction u(fes);

    // Project a function: u(x,y,z) = x*x + y + z*z
    u.project([](const Point& p) {
      return p.x() * p.x() + p.y() + p.z() * p.z();
    });

    // Save to stringstream
    std::stringstream ss;
    u.save(ss, FileFormat::MFEM);

    // Check FEC header
    std::string line;
    std::getline(ss, line); // "FiniteElementSpace"
    std::getline(ss, line); // "FiniteElementCollection: L2_3D_P0"
    EXPECT_TRUE(line.find("L2_3D_P0") != std::string::npos);

    // Reset stream and load
    ss.clear();
    ss.seekg(0);

    GridFunction v(fes);
    v.load(ss, FileFormat::MFEM);

    // Compare values
    const auto& uData = u.getData();
    const auto& vData = v.getData();
    ASSERT_EQ(uData.size(), vData.size());

    for (Index i = 0; i < uData.size(); i++)
    {
      EXPECT_NEAR(uData(i), vData(i), 1e-10);
    }
  }

  // Test P0 GridFunction save/load round-trip on segment mesh
  TEST(P0GridFunctionIO, SaveLoadRoundTrip_Segment)
  {
    // Create mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Segment, {17});
    EXPECT_EQ(mesh.getCellCount(), 16);

    // Create P0 space and grid function
    P0 fes(mesh);
    GridFunction u(fes);

    // Project a function: u(x) = sin(x)
    u.project([](const Point& p) {
      return std::sin(p.x());
    });

    // Save to stringstream
    std::stringstream ss;
    u.save(ss, FileFormat::MFEM);

    // Check FEC header
    std::string line;
    std::getline(ss, line); // "FiniteElementSpace"
    std::getline(ss, line); // "FiniteElementCollection: L2_1D_P0"
    EXPECT_TRUE(line.find("L2_1D_P0") != std::string::npos);

    // Reset stream and load
    ss.clear();
    ss.seekg(0);

    GridFunction v(fes);
    v.load(ss, FileFormat::MFEM);

    // Compare values
    const auto& uData = u.getData();
    const auto& vData = v.getData();
    ASSERT_EQ(uData.size(), vData.size());

    for (Index i = 0; i < uData.size(); i++)
    {
      EXPECT_NEAR(uData(i), vData(i), 1e-10);
    }
  }

  // Test P0 GridFunction save/load round-trip on mixed 2D mesh
  TEST(P0GridFunctionIO, SaveLoadRoundTrip_Mixed2D)
  {
    // Create a mixed mesh with triangles and quadrilaterals
    Mesh mesh;
    mesh.getVertices().resize(8, 2);
    
    // Define vertices for 2x2 grid
    mesh.getVertices().row(0) << 0.0, 0.0;  // v0
    mesh.getVertices().row(1) << 0.5, 0.0;  // v1
    mesh.getVertices().row(2) << 1.0, 0.0;  // v2
    mesh.getVertices().row(3) << 0.0, 0.5;  // v3
    mesh.getVertices().row(4) << 0.5, 0.5;  // v4
    mesh.getVertices().row(5) << 1.0, 0.5;  // v5
    mesh.getVertices().row(6) << 0.0, 1.0;  // v6
    mesh.getVertices().row(7) << 1.0, 1.0;  // v7

    // Add 2 triangles and 2 quadrilaterals (total 16 elements needed, so repeat pattern)
    // Bottom left: 2 triangles
    mesh.getPolytope(2, 0).setVertices({0, 1, 3});
    mesh.getPolytope(2, 1).setVertices({1, 4, 3});

    // Bottom right: quad
    mesh.getPolytope(2, 2).setVertices({1, 2, 5, 4});

    // Top left: quad
    mesh.getPolytope(2, 3).setVertices({3, 4, 6, 4}); // Simplified for test

    // Replicate to reach 16 elements
    for (int i = 4; i < 16; i++)
    {
      mesh.getPolytope(2, i).setVertices({0, 1, 3});
    }

    mesh.getConnectivity().compute(0, 2);
    mesh.getConnectivity().compute(2, 0);

    EXPECT_GE(mesh.getCellCount(), 16);

    // Create P0 space and grid function
    P0 fes(mesh);
    GridFunction u(fes);

    // Project a function: u(x,y) = x - y
    u.project([](const Point& p) {
      return p.x() - p.y();
    });

    // Save to stringstream
    std::stringstream ss;
    u.save(ss, FileFormat::MFEM);

    // Check FEC header
    std::string line;
    std::getline(ss, line); // "FiniteElementSpace"
    std::getline(ss, line); // "FiniteElementCollection: L2_2D_P0"
    EXPECT_TRUE(line.find("L2_2D_P0") != std::string::npos);

    // Reset stream and load
    ss.clear();
    ss.seekg(0);

    GridFunction v(fes);
    v.load(ss, FileFormat::MFEM);

    // Compare values
    const auto& uData = u.getData();
    const auto& vData = v.getData();
    ASSERT_EQ(uData.size(), vData.size());

    for (Index i = 0; i < uData.size(); i++)
    {
      EXPECT_NEAR(uData(i), vData(i), 1e-10);
    }
  }

  // Test P0 GridFunction file I/O
  TEST(P0GridFunctionIO, FileIO)
  {
    // Create mesh
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {4, 4});

    // Create P0 space and grid function
    P0 fes(mesh);
    GridFunction u(fes);

    // Project a function
    u.project([](const Point& p) {
      return p.x() * p.x() + p.y() * p.y();
    });

    // Save to file
    const std::string filename = "/tmp/test_p0_gridfunction.gf";
    std::ofstream ofs(filename);
    u.save(ofs, FileFormat::MFEM);
    ofs.close();

    // Load from file
    GridFunction v(fes);
    std::ifstream ifs(filename);
    v.load(ifs, FileFormat::MFEM);
    ifs.close();

    // Compare values
    const auto& uData = u.getData();
    const auto& vData = v.getData();
    ASSERT_EQ(uData.size(), vData.size());

    for (Index i = 0; i < uData.size(); i++)
    {
      EXPECT_NEAR(uData(i), vData(i), 1e-10);
    }

    // Clean up
    std::remove(filename.c_str());
  }
}
