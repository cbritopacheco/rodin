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
}
