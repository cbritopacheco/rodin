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
   * @brief Test saving and loading H1 GridFunction with degree 1 (scalar)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree1_Scalar)
  {
    // Create a simple 2D mesh
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });

    // Compute connectivity required for H1 spaces
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Create H1 space of degree 1
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction gf(fes);

    // Set the grid function to a simple linear function
    RealFunction linear_func([](const Geometry::Point& p) { return p.x() + 2.0 * p.y(); });
    gf.project(linear_func);

    // Save to stringstream
    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    // Verify the header content
    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");

    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P1");

    std::getline(ss, line);
    EXPECT_EQ(line, "VDim: 1");

    std::getline(ss, line);
    EXPECT_EQ(line, "Ordering: 0");

    // Reset stream position
    ss.clear();
    ss.seekg(0);

    // Load into a new grid function
    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    // Verify the data matches
    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test saving and loading H1 GridFunction with degree 2 (scalar)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree2_Scalar)
  {
    // Create a simple 2D mesh
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });

    // Compute connectivity required for H1 spaces
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Create H1 space of degree 2
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Set the grid function to a quadratic function
    RealFunction quad_func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y(); 
    });
    gf.project(quad_func);

    // Save to stringstream
    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    // Verify the header content
    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");

    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P2");

    std::getline(ss, line);
    EXPECT_EQ(line, "VDim: 1");

    std::getline(ss, line);
    EXPECT_EQ(line, "Ordering: 0");

    // Reset stream position
    ss.clear();
    ss.seekg(0);

    // Load into a new grid function
    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    // Verify the data matches
    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test saving and loading H1 GridFunction with degree 3 (scalar)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree3_Scalar)
  {
    // Create a simple 2D mesh
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });

    // Compute connectivity required for H1 spaces
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Create H1 space of degree 3
    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    // Set the grid function to a cubic function
    RealFunction cubic_func([](const Geometry::Point& p) { 
      return p.x() * p.x() * p.x() + p.y() * p.y() * p.y(); 
    });
    gf.project(cubic_func);

    // Save to stringstream
    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<3, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    // Verify the header content
    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");

    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P3");

    std::getline(ss, line);
    EXPECT_EQ(line, "VDim: 1");

    std::getline(ss, line);
    EXPECT_EQ(line, "Ordering: 0");

    // Reset stream position
    ss.clear();
    ss.seekg(0);

    // Load into a new grid function
    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<3, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    // Verify the data matches (relaxed tolerance for higher degree projections)
    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-5);
    }
  }

  /**
   * @brief Test saving H1 GridFunction to file and loading it back
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadFile_H1_Degree2_Scalar)
  {
    // Create a simple 2D mesh
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });

    // Compute connectivity required for H1 spaces
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // Create H1 space of degree 2
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    // Set the grid function to a quadratic function
    RealFunction quad_func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y(); 
    });
    gf.project(quad_func);

    // Save to file
    const std::string filename = "/tmp/test_h1_degree2.gf";
    {
      std::ofstream ofs(filename);
      GridFunctionPrinter<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> printer(gf);
      printer.print(ofs);
    }

    // Load from file
    GridFunction gf_loaded(fes);
    {
      std::ifstream ifs(filename);
      GridFunctionLoader<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> loader(gf_loaded);
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

  /**
   * @brief Test H1 degree 4 on triangle mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree4_Triangle)
  {
    // Create 2D triangle mesh with 32 elements (4x4 grid = 32 triangles)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    // Use a polynomial of degree 4
    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 4) + std::pow(p.y(), 4) + p.x() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    // Verify header
    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P4");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 5 on triangle mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree5_Triangle)
  {
    // Create 2D triangle mesh with 32 elements
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 5>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 5) + std::pow(p.y(), 5) + p.x() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<5, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P5");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<5, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 6 on triangle mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree6_Triangle)
  {
    // Create 2D triangle mesh with 32 elements
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 6) + std::pow(p.y(), 6) + p.x() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<6, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P6");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<6, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 4 on quadrilateral mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree4_Quadrilateral)
  {
    // Create 2D quad mesh with 9 elements (3x3 grid of NODES)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 9u);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 4) + std::pow(p.y(), 4) + p.x() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P4");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 5 on quadrilateral mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree5_Quadrilateral)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 9u);

    H1 fes(std::integral_constant<size_t, 5>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 5) + std::pow(p.y(), 5) + p.x() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<5, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P5");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<5, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 6 on quadrilateral mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree6_Quadrilateral)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 9u);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 6) + std::pow(p.y(), 6) + p.x() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<6, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P6");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<6, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 4 on tetrahedron mesh (3D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree4_Tetrahedron)
  {
    // Create 3D tet mesh with >= 16 elements (3x3x3 grid = 162 tets)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 4) + std::pow(p.y(), 4) + std::pow(p.z(), 4); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P4");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 5 on tetrahedron mesh (3D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree5_Tetrahedron)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 5>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 5) + std::pow(p.y(), 5) + std::pow(p.z(), 5); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<5, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P5");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<5, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 6 on tetrahedron mesh (3D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree6_Tetrahedron)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 6) + std::pow(p.y(), 6) + std::pow(p.z(), 6); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<6, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P6");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<6, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 4 on mixed 2D mesh (triangles and quads)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree4_Mixed2D)
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

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 4) + std::pow(p.y(), 4) + p.x() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P4");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 5 on mixed 2D mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree5_Mixed2D)
  {
    // Reuse the same mixed mesh structure
    Mesh mesh = Mesh<Context::Local>::Builder()
      .initialize(2)
      .nodes(25)
      .vertex({0.0, 0.0}).vertex({0.25, 0.0}).vertex({0.5, 0.0}).vertex({0.75, 0.0}).vertex({1.0, 0.0})
      .vertex({0.0, 0.25}).vertex({0.25, 0.25}).vertex({0.5, 0.25}).vertex({0.75, 0.25}).vertex({1.0, 0.25})
      .vertex({0.0, 0.5}).vertex({0.25, 0.5}).vertex({0.5, 0.5}).vertex({0.75, 0.5}).vertex({1.0, 0.5})
      .vertex({0.0, 0.75}).vertex({0.25, 0.75}).vertex({0.5, 0.75}).vertex({0.75, 0.75}).vertex({1.0, 0.75})
      .vertex({0.0, 1.0}).vertex({0.25, 1.0}).vertex({0.5, 1.0}).vertex({0.75, 1.0}).vertex({1.0, 1.0})
      .polytope(Polytope::Type::Triangle, {0, 1, 5})
      .polytope(Polytope::Type::Triangle, {1, 6, 5})
      .polytope(Polytope::Type::Triangle, {1, 2, 6})
      .polytope(Polytope::Type::Triangle, {2, 7, 6})
      .polytope(Polytope::Type::Quadrilateral, {5, 6, 11, 10})
      .polytope(Polytope::Type::Quadrilateral, {6, 7, 12, 11})
      .polytope(Polytope::Type::Quadrilateral, {7, 8, 13, 12})
      .polytope(Polytope::Type::Quadrilateral, {8, 9, 14, 13})
      .polytope(Polytope::Type::Triangle, {10, 11, 15})
      .polytope(Polytope::Type::Triangle, {11, 16, 15})
      .polytope(Polytope::Type::Triangle, {11, 12, 16})
      .polytope(Polytope::Type::Triangle, {12, 17, 16})
      .polytope(Polytope::Type::Quadrilateral, {15, 16, 21, 20})
      .polytope(Polytope::Type::Quadrilateral, {16, 17, 22, 21})
      .polytope(Polytope::Type::Quadrilateral, {17, 18, 23, 22})
      .polytope(Polytope::Type::Quadrilateral, {18, 19, 24, 23})
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 5>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 5) + std::pow(p.y(), 5) + p.x() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<5, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P5");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<5, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 6 on mixed 2D mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree6_Mixed2D)
  {
    Mesh mesh = Mesh<Context::Local>::Builder()
      .initialize(2)
      .nodes(25)
      .vertex({0.0, 0.0}).vertex({0.25, 0.0}).vertex({0.5, 0.0}).vertex({0.75, 0.0}).vertex({1.0, 0.0})
      .vertex({0.0, 0.25}).vertex({0.25, 0.25}).vertex({0.5, 0.25}).vertex({0.75, 0.25}).vertex({1.0, 0.25})
      .vertex({0.0, 0.5}).vertex({0.25, 0.5}).vertex({0.5, 0.5}).vertex({0.75, 0.5}).vertex({1.0, 0.5})
      .vertex({0.0, 0.75}).vertex({0.25, 0.75}).vertex({0.5, 0.75}).vertex({0.75, 0.75}).vertex({1.0, 0.75})
      .vertex({0.0, 1.0}).vertex({0.25, 1.0}).vertex({0.5, 1.0}).vertex({0.75, 1.0}).vertex({1.0, 1.0})
      .polytope(Polytope::Type::Triangle, {0, 1, 5})
      .polytope(Polytope::Type::Triangle, {1, 6, 5})
      .polytope(Polytope::Type::Triangle, {1, 2, 6})
      .polytope(Polytope::Type::Triangle, {2, 7, 6})
      .polytope(Polytope::Type::Quadrilateral, {5, 6, 11, 10})
      .polytope(Polytope::Type::Quadrilateral, {6, 7, 12, 11})
      .polytope(Polytope::Type::Quadrilateral, {7, 8, 13, 12})
      .polytope(Polytope::Type::Quadrilateral, {8, 9, 14, 13})
      .polytope(Polytope::Type::Triangle, {10, 11, 15})
      .polytope(Polytope::Type::Triangle, {11, 16, 15})
      .polytope(Polytope::Type::Triangle, {11, 12, 16})
      .polytope(Polytope::Type::Triangle, {12, 17, 16})
      .polytope(Polytope::Type::Quadrilateral, {15, 16, 21, 20})
      .polytope(Polytope::Type::Quadrilateral, {16, 17, 22, 21})
      .polytope(Polytope::Type::Quadrilateral, {17, 18, 23, 22})
      .polytope(Polytope::Type::Quadrilateral, {18, 19, 24, 23})
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 6>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 6) + std::pow(p.y(), 6) + p.x() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<6, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P6");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<6, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 1 on quadrilateral mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree1_Quadrilateral)
  {
    // Create 2D quad mesh with 16 elements (4x4 grid)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 9u);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P1");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 2 on quadrilateral mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree2_Quadrilateral)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 9u);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P2");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 1 on tetrahedron mesh (3D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree1_Tetrahedron)
  {
    // Create 3D tet mesh with >= 16 elements (3x3x3 grid = 162 tets)
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y() + 3.0 * p.z(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P1");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 2 on tetrahedron mesh (3D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree2_Tetrahedron)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, { 3, 3, 3 });

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y() + p.z() * p.z(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P2");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree1_Triangle_Large)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_EQ(mesh.getCellCount(), 18u);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P1");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 2 on triangle mesh with >= 16 elements
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree2_Triangle_Large)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P2");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 1 on mixed 2D mesh (triangles and quads)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree1_Mixed2D)
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

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P1");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 2 on mixed 2D mesh
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree2_Mixed2D)
  {
    // Reuse the same mixed mesh structure
    Mesh mesh = Mesh<Context::Local>::Builder()
      .initialize(2)
      .nodes(25)
      .vertex({0.0, 0.0}).vertex({0.25, 0.0}).vertex({0.5, 0.0}).vertex({0.75, 0.0}).vertex({1.0, 0.0})
      .vertex({0.0, 0.25}).vertex({0.25, 0.25}).vertex({0.5, 0.25}).vertex({0.75, 0.25}).vertex({1.0, 0.25})
      .vertex({0.0, 0.5}).vertex({0.25, 0.5}).vertex({0.5, 0.5}).vertex({0.75, 0.5}).vertex({1.0, 0.5})
      .vertex({0.0, 0.75}).vertex({0.25, 0.75}).vertex({0.5, 0.75}).vertex({0.75, 0.75}).vertex({1.0, 0.75})
      .vertex({0.0, 1.0}).vertex({0.25, 1.0}).vertex({0.5, 1.0}).vertex({0.75, 1.0}).vertex({1.0, 1.0})
      .polytope(Polytope::Type::Triangle, {0, 1, 5})
      .polytope(Polytope::Type::Triangle, {1, 6, 5})
      .polytope(Polytope::Type::Triangle, {1, 2, 6})
      .polytope(Polytope::Type::Triangle, {2, 7, 6})
      .polytope(Polytope::Type::Quadrilateral, {5, 6, 11, 10})
      .polytope(Polytope::Type::Quadrilateral, {6, 7, 12, 11})
      .polytope(Polytope::Type::Quadrilateral, {7, 8, 13, 12})
      .polytope(Polytope::Type::Quadrilateral, {8, 9, 14, 13})
      .polytope(Polytope::Type::Triangle, {10, 11, 15})
      .polytope(Polytope::Type::Triangle, {11, 16, 15})
      .polytope(Polytope::Type::Triangle, {11, 12, 16})
      .polytope(Polytope::Type::Triangle, {12, 17, 16})
      .polytope(Polytope::Type::Quadrilateral, {15, 16, 21, 20})
      .polytope(Polytope::Type::Quadrilateral, {16, 17, 22, 21})
      .polytope(Polytope::Type::Quadrilateral, {17, 18, 23, 22})
      .polytope(Polytope::Type::Quadrilateral, {18, 19, 24, 23})
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_2D_P2");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 1 on segment mesh (1D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree1_Segment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, { 16 });

    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 15u);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return 2.0 * p.x() + 1.0; 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_1D_P1");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 2 on segment mesh (1D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree2_Segment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, { 16 });

    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 15u);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() * p.x() + 2.0 * p.x() + 1.0; 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_1D_P2");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 3 on segment mesh (1D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree3_Segment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, { 16 });

    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 15u);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 3) + p.x() * p.x() + p.x(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<3, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_1D_P3");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<3, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 4 on segment mesh (1D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree4_Segment)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, { 16 });

    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 15u);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 4) + std::pow(p.x(), 3) + p.x(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_1D_P4");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 1 on wedge mesh (3D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree1_Wedge)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 3, 3, 3 });

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() + 2.0 * p.y() + 3.0 * p.z(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P1");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 2 on wedge mesh (3D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree2_Wedge)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 3, 3, 3 });

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return p.x() * p.x() + p.y() * p.y() + p.z() * p.z(); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P2");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<2, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
    }
  }

  /**
   * @brief Test H1 degree 3 on wedge mesh (3D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree3_Wedge)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 3, 3, 3 });

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 3>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 3) + std::pow(p.y(), 3) + std::pow(p.z(), 3); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<3, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P3");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<3, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }

  /**
   * @brief Test H1 degree 4 on wedge mesh (3D)
   */
  TEST(Rodin_IO_MFEM_H1_GridFunction, SaveLoadRoundTrip_H1_Degree4_Wedge)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 3, 3, 3 });

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GE(mesh.getCellCount(), 16u);

    H1 fes(std::integral_constant<size_t, 4>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p) { 
      return std::pow(p.x(), 4) + std::pow(p.y(), 4) + std::pow(p.z(), 4); 
    });
    gf.project(func);

    std::stringstream ss;
    GridFunctionPrinter<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> printer(gf);
    printer.print(ss);

    std::string line;
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementSpace");
    std::getline(ss, line);
    EXPECT_EQ(line, "FiniteElementCollection: H1_3D_P4");

    ss.clear();
    ss.seekg(0);

    GridFunction gf_loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, H1<4, Real>, Math::Vector<Real>> loader(gf_loaded);
    loader.load(ss);

    ASSERT_EQ(gf.getSize(), gf_loaded.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-4);
    }
  }
}
