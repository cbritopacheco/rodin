/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <sstream>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/IO.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::IO;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    Mesh<Context::Local> makePointMesh()
    {
      Mesh<Context::Local>::Builder builder;
      return builder
        .initialize(1)
        .nodes(1)
        .vertex({ 0.0 })
        .finalize();
    }

    Mesh<Context::Local> makeMesh(Polytope::Type geometry)
    {
      switch (geometry)
      {
        case Polytope::Type::Point:
          return makePointMesh();
        case Polytope::Type::Segment:
          return LocalMesh::UniformGrid(geometry, { 5 });
        case Polytope::Type::Triangle:
        case Polytope::Type::Quadrilateral:
          return LocalMesh::UniformGrid(geometry, { 3, 3 });
        case Polytope::Type::Tetrahedron:
        case Polytope::Type::Pyramid:
        case Polytope::Type::Hexahedron:
        case Polytope::Type::Wedge:
          return LocalMesh::UniformGrid(geometry, { 2, 2, 2 });
      }
      assert(false);
      return {};
    }

    std::string geometryName(Polytope::Type geometry)
    {
      switch (geometry)
      {
        case Polytope::Type::Point:         return "Point";
        case Polytope::Type::Segment:       return "Segment";
        case Polytope::Type::Triangle:      return "Triangle";
        case Polytope::Type::Quadrilateral: return "Quadrilateral";
        case Polytope::Type::Tetrahedron:   return "Tetrahedron";
        case Polytope::Type::Pyramid:       return "Pyramid";
        case Polytope::Type::Hexahedron:    return "Hexahedron";
        case Polytope::Type::Wedge:         return "Wedge";
      }
      return "Unknown";
    }
  }

  class P0MFEMGridFunctionCoverage : public ::testing::TestWithParam<Polytope::Type>
  {};

  /// @brief Verifies scalar round trip for P0 MFEM grid function coverage by checking tolerance-based numerical results, exact expected values.
  TEST_P(P0MFEMGridFunctionCoverage, ScalarRoundTrip)
  {
    Mesh mesh = makeMesh(GetParam());
    P0 fes(mesh);
    GridFunction gf(fes);

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); ++i)
      gf[i] = 1.25 + static_cast<Real>(i);

    std::stringstream stream;
    GridFunctionPrinter<FileFormat::MFEM, P0<Real>, Math::Vector<Real>> printer(gf);
    printer.print(stream);

    GridFunction loaded(fes);
    GridFunctionLoader<FileFormat::MFEM, P0<Real>, Math::Vector<Real>> loader(loaded);
    loader.load(stream);

    ASSERT_EQ(loaded.getSize(), gf.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); ++i)
      EXPECT_DOUBLE_EQ(loaded[i], gf[i]) << "dof " << i;
  }

  /// @brief Instantiates P0 MFEM Grid Function Coverage over the All Supported Geometries parameter coverage.
  INSTANTIATE_TEST_SUITE_P(
      AllSupportedGeometries,
      P0MFEMGridFunctionCoverage,
      ::testing::Values(
        Polytope::Type::Point,
        Polytope::Type::Segment,
        Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral,
        Polytope::Type::Tetrahedron,
        Polytope::Type::Pyramid,
        Polytope::Type::Hexahedron,
        Polytope::Type::Wedge),
      [](const ::testing::TestParamInfo<Polytope::Type>& info)
      {
        return geometryName(info.param);
      });

  /// @brief Verifies scalar multiline load for IO MFEM P0 grid function by checking tolerance-based numerical results.
  TEST(Rodin_IO_MFEM_P0_GridFunction, ScalarMultilineLoad)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, { 4 });

    P0 fes(mesh);
    GridFunction gf(fes);

    std::stringstream stream;
    stream
      << "FiniteElementSpace\n"
      << "FiniteElementCollection: L2_1D_P0\n"
      << "VDim: 1\n"
      << "Ordering: 1\n\n"
      << "1 2\n"
      << "3\n";

    GridFunctionLoader<FileFormat::MFEM, P0<Real>, Math::Vector<Real>> loader(gf);
    loader.load(stream);

    EXPECT_DOUBLE_EQ(gf[0], 1);
    EXPECT_DOUBLE_EQ(gf[1], 2);
    EXPECT_DOUBLE_EQ(gf[2], 3);
  }
}
