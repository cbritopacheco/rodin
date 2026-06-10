/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/path.hpp>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Rodin/Configure.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/IO.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    boost::filesystem::path mfemH1FixturePath(const std::string& filename)
    {
      boost::filesystem::path path(RODIN_RESOURCES_DIR);
      path /= "tests";
      path /= "mfem";
      path /= "h1";
      path /= filename;
      return path;
    }

    std::string mfemH1FixtureName(
        const std::string& geometry,
        size_t order,
        const std::string& field,
        const std::string& ordering,
        const std::string& extension)
    {
      return "h1_" + geometry + "_p" + std::to_string(order)
        + "_" + field + "_mfem_" + ordering + "." + extension;
    }

    Mesh<Context::Local> loadMFEMH1FixtureMesh(const std::string& filename)
    {
      boost::filesystem::ifstream in(mfemH1FixturePath(filename));
      if (!in)
        throw std::runtime_error(std::string("Could not open fixture mesh: ") + filename);

      Mesh mesh;
      MeshLoader<FileFormat::MFEM, Context::Local> loader(mesh);
      loader.load(in);

      for (size_t d = mesh.getDimension(); d > 0; --d)
        mesh.getConnectivity().compute(d, d - 1);

      return mesh;
    }

    Real mfemH1FixtureScalar(const Geometry::Point& p)
    {
      const Real x = p.getDimension() > 0 ? p(0) : 0.0;
      const Real y = p.getDimension() > 1 ? p(1) : 0.0;
      const Real z = p.getDimension() > 2 ? p(2) : 0.0;

      return 1.0 + 2.0 * x - 3.0 * y + 0.5 * z
        + x * x + 0.25 * y * y - 0.75 * z * z
        + x * y - 0.4 * x * z + 0.2 * y * z;
    }

    auto mfemH1FixtureScalarFunction()
    {
      return RealFunction([](const Geometry::Point& p)
      {
        return mfemH1FixtureScalar(p);
      });
    }

    auto mfemH1FixtureVectorFunction()
    {
      return VectorFunction{
        [](const Geometry::Point& p)
        {
          return mfemH1FixtureScalar(p);
        },
        [](const Geometry::Point& p)
        {
          const Real x = p.getDimension() > 0 ? p(0) : 0.0;
          const Real y = p.getDimension() > 1 ? p(1) : 0.0;
          const Real z = p.getDimension() > 2 ? p(2) : 0.0;
          return -2.0 + x + 4.0 * y - z + x * z + 0.5 * y * y;
        },
        [](const Geometry::Point& p)
        {
          const Real x = p.getDimension() > 0 ? p(0) : 0.0;
          const Real y = p.getDimension() > 1 ? p(1) : 0.0;
          const Real z = p.getDimension() > 2 ? p(2) : 0.0;
          return 3.0 - x * y + z * z + 0.25 * x * x;
        }
      };
    }

    std::vector<Real> readMFEMGridFunctionValues(std::istream& in)
    {
      std::string line;
      bool sawOrdering = false;

      while (std::getline(in, line))
      {
        if (!sawOrdering)
        {
          sawOrdering = line.rfind("Ordering:", 0) == 0;
          continue;
        }

        if (line.find_first_not_of(" \t\r\n") == std::string::npos)
          break;
      }

      std::vector<Real> values;
      Real value;
      while (in >> value)
        values.push_back(value);
      return values;
    }

    std::vector<Real> readMFEMGridFunctionValues(const boost::filesystem::path& path)
    {
      boost::filesystem::ifstream in(path);
      if (!in)
        throw std::runtime_error("Could not open fixture grid function: " + path.string());
      return readMFEMGridFunctionValues(in);
    }

    std::vector<Real> readMFEMGridFunctionValues(const std::string& contents)
    {
      std::istringstream in(contents);
      return readMFEMGridFunctionValues(in);
    }

    template <class GF1, class GF2>
    void expectGridFunctionNear(const GF1& actual, const GF2& expected, Real tolerance)
    {
      ASSERT_EQ(actual.getSize(), expected.getSize());
      for (Index i = 0; i < static_cast<Index>(actual.getSize()); ++i)
        EXPECT_NEAR(actual[i], expected[i], tolerance) << "global dof " << i;
    }

    template <class FES>
    void loadMFEMGridFunctionFixture(
        const std::string& gfFilename,
        GridFunction<FES, Math::Vector<Real>>& gf)
    {
      boost::filesystem::ifstream in(mfemH1FixturePath(gfFilename));
      ASSERT_TRUE(static_cast<bool>(in)) << gfFilename;
      GridFunctionLoader<FileFormat::MFEM, FES, Math::Vector<Real>> loader(gf);
      loader.load(in);
    }

    template <class FES>
    void roundTripMFEMGridFunctionAndCompare(
        const GridFunction<FES, Math::Vector<Real>>& gf,
        Real tolerance)
    {
      std::stringstream out;
      GridFunctionPrinter<FileFormat::MFEM, FES, Math::Vector<Real>> printer(gf);
      printer.print(out);

      GridFunction<FES, Math::Vector<Real>> reloaded(gf.getFiniteElementSpace());
      GridFunctionLoader<FileFormat::MFEM, FES, Math::Vector<Real>> loader(reloaded);
      loader.load(out);

      expectGridFunctionNear(reloaded, gf, tolerance);
    }

    template <class FES>
    void loadMFEMScalarFixtureAndCompare(
        const std::string& gfFilename,
        const FES& fes,
        Real tolerance)
    {
      GridFunction<FES, Math::Vector<Real>> gf(fes);
      loadMFEMGridFunctionFixture(gfFilename, gf);

      GridFunction expected(fes);
      expected.project(mfemH1FixtureScalarFunction());

      expectGridFunctionNear(gf, expected, tolerance);
    }

    template <size_t K>
    void loadMFEMScalarFixtureCaseAndCompare(
        const std::string& geometry,
        Real tolerance)
    {
      const auto meshFilename =
        mfemH1FixtureName(geometry, K, "scalar", "by_nodes", "mesh");
      const auto gfFilename =
        mfemH1FixtureName(geometry, K, "scalar", "by_nodes", "gf");

      Mesh mesh = loadMFEMH1FixtureMesh(meshFilename);
      H1 fes(std::integral_constant<size_t, K>{}, mesh);

      loadMFEMScalarFixtureAndCompare(gfFilename, fes, tolerance);
    }

    template <size_t K>
    void roundTripMFEMScalarFixtureCaseAndCompare(
        const std::string& geometry,
        Real tolerance)
    {
      const auto meshFilename =
        mfemH1FixtureName(geometry, K, "scalar", "by_nodes", "mesh");
      const auto gfFilename =
        mfemH1FixtureName(geometry, K, "scalar", "by_nodes", "gf");

      Mesh mesh = loadMFEMH1FixtureMesh(meshFilename);
      H1<K, Real> fes(std::integral_constant<size_t, K>{}, mesh);

      GridFunction<H1<K, Real>, Math::Vector<Real>> gf(fes);
      loadMFEMGridFunctionFixture(gfFilename, gf);
      roundTripMFEMGridFunctionAndCompare(gf, tolerance);
    }

    template <size_t K>
    void loadMFEMScalarFileFixtureFamilyAndCompare()
    {
      loadMFEMScalarFixtureCaseAndCompare<K>("triangle", 1e-9);
      loadMFEMScalarFixtureCaseAndCompare<K>("mixed2d", 1e-9);
      loadMFEMScalarFixtureCaseAndCompare<K>("tetrahedron", 1e-9);
      loadMFEMScalarFixtureCaseAndCompare<K>("hexahedron", 1e-9);
      loadMFEMScalarFixtureCaseAndCompare<K>("wedge", 1e-9);
    }

    template <size_t K>
    void roundTripMFEMScalarFileFixtureFamilyAndCompare()
    {
      roundTripMFEMScalarFixtureCaseAndCompare<K>("triangle", 1e-9);
      roundTripMFEMScalarFixtureCaseAndCompare<K>("mixed2d", 1e-9);
      roundTripMFEMScalarFixtureCaseAndCompare<K>("tetrahedron", 1e-9);
      roundTripMFEMScalarFixtureCaseAndCompare<K>("hexahedron", 1e-9);
      roundTripMFEMScalarFixtureCaseAndCompare<K>("wedge", 1e-9);
    }

    Mesh<Context::Local> makeUniform3DMeshForH1RoundTrip(Polytope::Type geometry)
    {
      return LocalMesh::UniformGrid(geometry, { 3, 3, 3 });
    }

    std::string geometryName(Polytope::Type geometry)
    {
      switch (geometry)
      {
        case Polytope::Type::Tetrahedron: return "Tetrahedron";
        case Polytope::Type::Hexahedron:  return "Hexahedron";
        case Polytope::Type::Pyramid:     return "Pyramid";
        case Polytope::Type::Wedge:       return "Wedge";
        default:                          return "Unknown";
      }
    }

    template <size_t K>
    void loadMFEMVectorFixtureAndCompare(
        const std::string& meshFilename,
        const std::string& gfFilename,
        Real tolerance)
    {
      Mesh mesh = loadMFEMH1FixtureMesh(meshFilename);

      using FES = H1<K, Math::SpatialVector<Real>>;
      FES fes(std::integral_constant<size_t, K>{}, mesh, 3);

      GridFunction<FES, Math::Vector<Real>> gf(fes);
      loadMFEMGridFunctionFixture(gfFilename, gf);

      GridFunction expected(fes);
      expected.project(mfemH1FixtureVectorFunction());
      expectGridFunctionNear(gf, expected, tolerance);
    }

    template <size_t K>
    void loadAndPrintMFEMVectorFixtureAndCompare(
        const std::string& meshFilename,
        const std::string& gfFilename,
        const std::string& expectedGfFilename,
        Real tolerance)
    {
      Mesh mesh = loadMFEMH1FixtureMesh(meshFilename);

      using FES = H1<K, Math::SpatialVector<Real>>;
      FES fes(std::integral_constant<size_t, K>{}, mesh, 3);

      GridFunction<FES, Math::Vector<Real>> gf(fes);
      loadMFEMGridFunctionFixture(gfFilename, gf);

      GridFunction expected(fes);
      expected.project(mfemH1FixtureVectorFunction());
      expectGridFunctionNear(gf, expected, tolerance);

      std::stringstream out;
      GridFunctionPrinter<FileFormat::MFEM, FES, Math::Vector<Real>> printer(gf);
      printer.print(out);

      EXPECT_NE(out.str().find("Ordering: 0"), std::string::npos);

      const auto expectedValues = readMFEMGridFunctionValues(mfemH1FixturePath(expectedGfFilename));
      const auto actualValues = readMFEMGridFunctionValues(out.str());

      ASSERT_EQ(actualValues.size(), expectedValues.size());
      for (size_t i = 0; i < actualValues.size(); ++i)
      {
        EXPECT_NEAR(actualValues[i], expectedValues[i], tolerance)
          << "MFEM stream value " << i;
      }
    }

    template <size_t K>
    void loadMFEMVectorFixtureCaseAndCompare(
        const std::string& geometry,
        const std::string& ordering,
        Real tolerance)
    {
      const auto meshFilename =
        mfemH1FixtureName(geometry, K, "vector", ordering, "mesh");
      const auto gfFilename =
        mfemH1FixtureName(geometry, K, "vector", ordering, "gf");

      loadMFEMVectorFixtureAndCompare<K>(
        meshFilename, gfFilename, tolerance);
    }

    template <size_t K>
    void roundTripMFEMVectorFixtureCaseAndCompare(
        const std::string& geometry,
        const std::string& ordering,
        Real tolerance)
    {
      const auto meshFilename =
        mfemH1FixtureName(geometry, K, "vector", ordering, "mesh");
      const auto gfFilename =
        mfemH1FixtureName(geometry, K, "vector", ordering, "gf");

      Mesh mesh = loadMFEMH1FixtureMesh(meshFilename);

      using FES = H1<K, Math::SpatialVector<Real>>;
      FES fes(std::integral_constant<size_t, K>{}, mesh, 3);

      GridFunction<FES, Math::Vector<Real>> gf(fes);
      loadMFEMGridFunctionFixture(gfFilename, gf);
      roundTripMFEMGridFunctionAndCompare(gf, tolerance);
    }

    template <size_t K>
    void loadAndPrintMFEMVectorFixtureCaseAndCompare(
        const std::string& geometry,
        Real tolerance)
    {
      const auto meshFilename =
        mfemH1FixtureName(geometry, K, "vector", "by_vdim", "mesh");
      const auto gfFilename =
        mfemH1FixtureName(geometry, K, "vector", "by_vdim", "gf");
      const auto expectedGfFilename =
        mfemH1FixtureName(geometry, K, "vector", "by_nodes", "gf");

      loadAndPrintMFEMVectorFixtureAndCompare<K>(
        meshFilename, gfFilename, expectedGfFilename, tolerance);
    }

    template <size_t K>
    void loadMFEMVectorFileFixtureFamilyAndCompare(const std::string& ordering)
    {
      loadMFEMVectorFixtureCaseAndCompare<K>("triangle", ordering, 1e-9);
      loadMFEMVectorFixtureCaseAndCompare<K>("mixed2d", ordering, 1e-9);
      loadMFEMVectorFixtureCaseAndCompare<K>("tetrahedron", ordering, 1e-9);
      loadMFEMVectorFixtureCaseAndCompare<K>("hexahedron", ordering, 1e-9);
      loadMFEMVectorFixtureCaseAndCompare<K>("wedge", ordering, 1e-9);
    }

    template <size_t K>
    void roundTripMFEMVectorFileFixtureFamilyAndCompare(const std::string& ordering)
    {
      roundTripMFEMVectorFixtureCaseAndCompare<K>("triangle", ordering, 1e-9);
      roundTripMFEMVectorFixtureCaseAndCompare<K>("mixed2d", ordering, 1e-9);
      roundTripMFEMVectorFixtureCaseAndCompare<K>("tetrahedron", ordering, 1e-9);
      roundTripMFEMVectorFixtureCaseAndCompare<K>("hexahedron", ordering, 1e-9);
      roundTripMFEMVectorFixtureCaseAndCompare<K>("wedge", ordering, 1e-9);
    }

    template <size_t K>
    void loadAndPrintMFEMVectorFileFixtureFamilyAndCompare()
    {
      loadAndPrintMFEMVectorFixtureCaseAndCompare<K>("triangle", 1e-9);
      loadAndPrintMFEMVectorFixtureCaseAndCompare<K>("mixed2d", 1e-9);
      loadAndPrintMFEMVectorFixtureCaseAndCompare<K>("tetrahedron", 1e-9);
      loadAndPrintMFEMVectorFixtureCaseAndCompare<K>("hexahedron", 1e-9);
      loadAndPrintMFEMVectorFixtureCaseAndCompare<K>("wedge", 1e-9);
    }
  }

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
    // EXPECT_EQ(line, "Ordering: 0");

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
    // EXPECT_EQ(line, "Ordering: 0");

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
    // EXPECT_EQ(line, "Ordering: 0");

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

  class H1Degree2RoundTrip3DGeometryCoverage
    : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(H1Degree2RoundTrip3DGeometryCoverage, SaveLoadRoundTrip)
  {
    const auto geometry = GetParam();
    Mesh mesh = makeUniform3DMeshForH1RoundTrip(geometry);

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    ASSERT_GT(mesh.getCellCount(), 0u);

    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);

    RealFunction func([](const Geometry::Point& p)
    {
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
      EXPECT_NEAR(gf[i], gf_loaded[i], 1e-10);
  }

  INSTANTIATE_TEST_SUITE_P(
    All3DPolytopes,
    H1Degree2RoundTrip3DGeometryCoverage,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info)
    {
      return geometryName(info.param);
    });

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

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMStringFixture_H1_Degree1_Triangle_Scalar_ByNodes)
  {
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    GridFunction gf(fes);

    const std::string fixture = R"(FiniteElementSpace
FiniteElementCollection: H1_2D_P1
VDim: 1
Ordering: 0

1
2
3
)";

    std::istringstream in(fixture);
    GridFunctionLoader<FileFormat::MFEM, H1<1, Real>, Math::Vector<Real>> loader(gf);
    loader.load(in);

    GridFunction expected(fes);
    expected.project(RealFunction([](const Geometry::Point& p)
    {
      return 1.0 + p.x() + 2.0 * p.y();
    }));

    expectGridFunctionNear(gf, expected, 1e-12);
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree1_Scalar_ByNodes)
  {
    loadMFEMScalarFileFixtureFamilyAndCompare<1>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree2_Scalar_ByNodes)
  {
    loadMFEMScalarFileFixtureFamilyAndCompare<2>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree3_Scalar_ByNodes)
  {
    loadMFEMScalarFileFixtureFamilyAndCompare<3>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree4_Scalar_ByNodes)
  {
    loadMFEMScalarFileFixtureFamilyAndCompare<4>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree5_Scalar_ByNodes)
  {
    loadMFEMScalarFileFixtureFamilyAndCompare<5>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree6_Scalar_ByNodes)
  {
    loadMFEMScalarFileFixtureFamilyAndCompare<6>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree1_Scalar)
  {
    roundTripMFEMScalarFileFixtureFamilyAndCompare<1>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree2_Scalar)
  {
    roundTripMFEMScalarFileFixtureFamilyAndCompare<2>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree3_Scalar)
  {
    roundTripMFEMScalarFileFixtureFamilyAndCompare<3>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree4_Scalar)
  {
    roundTripMFEMScalarFileFixtureFamilyAndCompare<4>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree5_Scalar)
  {
    roundTripMFEMScalarFileFixtureFamilyAndCompare<5>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree6_Scalar)
  {
    roundTripMFEMScalarFileFixtureFamilyAndCompare<6>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree1_Vector_ByNodes)
  {
    loadMFEMVectorFileFixtureFamilyAndCompare<1>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree2_Vector_ByNodes)
  {
    loadMFEMVectorFileFixtureFamilyAndCompare<2>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree3_Vector_ByNodes)
  {
    loadMFEMVectorFileFixtureFamilyAndCompare<3>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree4_Vector_ByNodes)
  {
    loadMFEMVectorFileFixtureFamilyAndCompare<4>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree5_Vector_ByNodes)
  {
    loadMFEMVectorFileFixtureFamilyAndCompare<5>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadMFEMFileFixtures_H1_Degree6_Vector_ByNodes)
  {
    loadMFEMVectorFileFixtureFamilyAndCompare<6>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree1_Vector_ByNodes)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<1>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree2_Vector_ByNodes)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<2>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree3_Vector_ByNodes)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<3>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree4_Vector_ByNodes)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<4>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree5_Vector_ByNodes)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<5>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree6_Vector_ByNodes)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<6>("by_nodes");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadByVDimPrintByNodesMFEMFileFixtures_H1_Degree1_Vector)
  {
    loadAndPrintMFEMVectorFileFixtureFamilyAndCompare<1>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadByVDimPrintByNodesMFEMFileFixtures_H1_Degree2_Vector)
  {
    loadAndPrintMFEMVectorFileFixtureFamilyAndCompare<2>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadByVDimPrintByNodesMFEMFileFixtures_H1_Degree3_Vector)
  {
    loadAndPrintMFEMVectorFileFixtureFamilyAndCompare<3>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadByVDimPrintByNodesMFEMFileFixtures_H1_Degree4_Vector)
  {
    loadAndPrintMFEMVectorFileFixtureFamilyAndCompare<4>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadByVDimPrintByNodesMFEMFileFixtures_H1_Degree5_Vector)
  {
    loadAndPrintMFEMVectorFileFixtureFamilyAndCompare<5>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, LoadByVDimPrintByNodesMFEMFileFixtures_H1_Degree6_Vector)
  {
    loadAndPrintMFEMVectorFileFixtureFamilyAndCompare<6>();
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree1_Vector_ByVDim)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<1>("by_vdim");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree2_Vector_ByVDim)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<2>("by_vdim");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree3_Vector_ByVDim)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<3>("by_vdim");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree4_Vector_ByVDim)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<4>("by_vdim");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree5_Vector_ByVDim)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<5>("by_vdim");
  }

  TEST(Rodin_IO_MFEM_H1_GridFunction, RoundTripMFEMFileFixtures_H1_Degree6_Vector_ByVDim)
  {
    roundTripMFEMVectorFileFixtureFamilyAndCompare<6>("by_vdim");
  }
}
