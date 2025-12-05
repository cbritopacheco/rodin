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
#include "Rodin/Variational.h"
#include "Rodin/IO.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /**
   * @brief Test saving and loading P1 GridFunction on triangle mesh (2D)
   */
  TEST(Rodin_IO_MFEM_P1_GridFunction, SaveLoadRoundTrip_Triangle)
  {
    // Create 2D triangle mesh with >= 16 elements (4x4 grid = 32 triangles)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    ASSERT_GE(mesh.getCellCount(), 16u);
    
    P1 fes(mesh);
    GridFunction gf(fes);
    
    // P1 is linear, so use a linear function
    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y(); 
    });
    gf.project(func);
    
    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);
    
    // Verify header
    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P1");
    
    ss.clear();
    ss.seekg(0);
    
    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);
    
    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test saving and loading P1 GridFunction on quadrilateral mesh (2D)
   */
  TEST(Rodin_IO_MFEM_P1_GridFunction, SaveLoadRoundTrip_Quadrilateral)
  {
    // Create 2D quad mesh with >= 16 elements (5x5 grid = 16 quads)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 5, 5 });
    
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    ASSERT_GE(mesh.getCellCount(), 16u);
    
    P1 fes(mesh);
    GridFunction gf(fes);
    
    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y(); 
    });
    gf.project(func);
    
    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);
    
    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P1");
    
    ss.clear();
    ss.seekg(0);
    
    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);
    
    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test saving and loading P1 GridFunction on tetrahedron mesh (3D)
   */
  TEST(Rodin_IO_MFEM_P1_GridFunction, SaveLoadRoundTrip_Tetrahedron)
  {
    // Create 3D tet mesh with >= 16 elements (3x3x3 grid = 162 tets)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });
    
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    ASSERT_GE(mesh.getCellCount(), 16u);
    
    P1 fes(mesh);
    GridFunction gf(fes);
    
    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y() + 3.0 * p.z(); 
    });
    gf.project(func);
    
    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);
    
    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P1");
    
    ss.clear();
    ss.seekg(0);
    
    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);
    
    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test saving and loading P1 GridFunction on wedge mesh (3D)
   */
  TEST(Rodin_IO_MFEM_P1_GridFunction, SaveLoadRoundTrip_Wedge)
  {
    // Create 3D wedge mesh with >= 16 elements (3x3x3 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 3, 3, 3 });
    
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    ASSERT_GE(mesh.getCellCount(), 16u);
    
    P1 fes(mesh);
    GridFunction gf(fes);
    
    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y() + 3.0 * p.z(); 
    });
    gf.project(func);
    
    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);
    
    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P1");
    
    ss.clear();
    ss.seekg(0);
    
    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);
    
    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test saving and loading P1 GridFunction on segment mesh (1D)
   */
  TEST(Rodin_IO_MFEM_P1_GridFunction, SaveLoadRoundTrip_Segment)
  {
    // Create 1D segment mesh with >= 16 elements (17 divisions = 16 segments)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, { 17 });
    
    mesh.getConnectivity().compute(1, 0);
    
    ASSERT_GE(mesh.getCellCount(), 16u);
    
    P1 fes(mesh);
    GridFunction gf(fes);
    
    RealFunction func([](const Geometry::Point& p) { 
      return 2.0 * p.x() + 1.0; 
    });
    gf.project(func);
    
    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);
    
    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_1D_P1");
    
    ss.clear();
    ss.seekg(0);
    
    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);
    
    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test saving and loading P1 GridFunction on mixed 2D mesh
   */
  TEST(Rodin_IO_MFEM_P1_GridFunction, SaveLoadRoundTrip_Mixed2D)
  {
    // Create a mixed mesh with triangles and quads, >= 16 elements
    Mesh mesh = Mesh<Context::Local>::Builder()
      .initialize(2)
      .nodes(25)  // 5x5 grid of nodes
      .vertex({0.0, 0.0}).vertex({0.25, 0.0}).vertex({0.5, 0.0}).vertex({0.75, 0.0}).vertex({1.0, 0.0})
      .vertex({0.0, 0.25}).vertex({0.25, 0.25}).vertex({0.5, 0.25}).vertex({0.75, 0.25}).vertex({1.0, 0.25})
      .vertex({0.0, 0.5}).vertex({0.25, 0.5}).vertex({0.5, 0.5}).vertex({0.75, 0.5}).vertex({1.0, 0.5})
      .vertex({0.0, 0.75}).vertex({0.25, 0.75}).vertex({0.5, 0.75}).vertex({0.75, 0.75}).vertex({1.0, 0.75})
      .vertex({0.0, 1.0}).vertex({0.25, 1.0}).vertex({0.5, 1.0}).vertex({0.75, 1.0}).vertex({1.0, 1.0})
      // First row: triangles
      .polytope(Polytope::Type::Triangle, {0, 1, 5})
      .polytope(Polytope::Type::Triangle, {1, 6, 5})
      .polytope(Polytope::Type::Triangle, {1, 2, 6})
      .polytope(Polytope::Type::Triangle, {2, 7, 6})
      // Second row: quads
      .polytope(Polytope::Type::Quadrilateral, {5, 6, 11, 10})
      .polytope(Polytope::Type::Quadrilateral, {6, 7, 12, 11})
      .polytope(Polytope::Type::Quadrilateral, {7, 8, 13, 12})
      .polytope(Polytope::Type::Quadrilateral, {8, 9, 14, 13})
      // Third row: triangles
      .polytope(Polytope::Type::Triangle, {10, 11, 15})
      .polytope(Polytope::Type::Triangle, {11, 16, 15})
      .polytope(Polytope::Type::Triangle, {11, 12, 16})
      .polytope(Polytope::Type::Triangle, {12, 17, 16})
      // Fourth row: quads
      .polytope(Polytope::Type::Quadrilateral, {15, 16, 21, 20})
      .polytope(Polytope::Type::Quadrilateral, {16, 17, 22, 21})
      .polytope(Polytope::Type::Quadrilateral, {17, 18, 23, 22})
      .polytope(Polytope::Type::Quadrilateral, {18, 19, 24, 23})
      .finalize();
    
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    ASSERT_GE(mesh.getCellCount(), 16u);
    
    P1 fes(mesh);
    GridFunction gf(fes);
    
    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y(); 
    });
    gf.project(func);
    
    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);
    
    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P1");
    
    ss.clear();
    ss.seekg(0);
    
    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);
    
    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test file-based I/O for P1 GridFunction
   */
  TEST(Rodin_IO_MFEM_P1_GridFunction, SaveLoadFile)
  {
    // Create a simple 2D mesh
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    
    P1 fes(mesh);
    GridFunction gf(fes);
    
    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y(); 
    });
    gf.project(func);
    
    // Save to file
    const std::string filename = "/tmp/test_p1.gf";
    {
      std::ofstream ofs(filename);
      GridFunctionPrinter<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> printer(gf);
      printer.print(ofs);
    }
    
    // Load from file
    GridFunction gf_loaded(fes);
    {
      std::ifstream ifs(filename);
      GridFunctionLoader<FileFormat::MFEM, P1<Real>, Math::Vector<Real>> loader(gf_loaded);
      loader.load(ifs);
    }
    
    // Verify the data matches
    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
    
    // Clean up
    std::remove(filename.c_str());
  }
}
