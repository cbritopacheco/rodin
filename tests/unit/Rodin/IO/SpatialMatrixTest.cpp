/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Matrix field persistence preserves entries and rectangular shape. */
#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <sstream>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/IO.h"
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::IO;
namespace
{
  template <class FES>
  void checkPersistence(const FES& fes)
  {
    SCOPED_TRACE(typeid(FES).name());
    GridFunction field(fes), loaded(fes);
    field = [&](const Point& p) {
      typename FormLanguage::Traits<FES>::RangeType value(
        fes.getRows(), fes.getColumns());
      for (size_t r = 0; r < fes.getRows(); ++r)
        for (size_t c = 0; c < fes.getColumns(); ++c)
        {
          using Scalar = typename FormLanguage::Traits<FES>::ScalarType;
          value(r, c) = Scalar(1 + 7 * r + c + p.x());
          if constexpr (std::is_same_v<Scalar, Complex>)
            value(r, c) += Complex(0, 2 + r - Real(c));
        }
      return value;
    };
    const std::string filename = "/tmp/rodin_matrix_field.h5";
    field.save(filename, FileFormat::HDF5);
    loaded.load(filename, FileFormat::HDF5);
    EXPECT_LE((field.getData() - loaded.getData()).norm(), 1e-12);
    std::remove(filename.c_str());
  }
  template <class FES>
  void checkVertexFormat(const FES& fes, FileFormat format = FileFormat::MEDIT)
  {
    SCOPED_TRACE(typeid(FES).name());
    GridFunction field(fes), loaded(fes);
    field = [&](const Point& p) {
      Math::SpatialMatrix<Real> value(fes.getRows(), fes.getColumns());
      for (size_t r = 0; r < fes.getRows(); ++r)
        for (size_t c = 0; c < fes.getColumns(); ++c)
          value(r, c) = 1 + 7 * r + c + p.x();
      return value;
    };
    std::stringstream stream;
    if (format == FileFormat::MEDIT)
    {
      GridFunctionPrinter<FileFormat::MEDIT, FES, Math::Vector<Real>>(field).print(
        stream);
      GridFunctionLoader<FileFormat::MEDIT, FES, Math::Vector<Real>>(loaded).load(stream);
    }
    else
    {
      GridFunctionPrinter<FileFormat::MFEM, FES, Math::Vector<Real>>(field).print(stream);
      GridFunctionLoader<FileFormat::MFEM, FES, Math::Vector<Real>>(loaded).load(stream);
    }
    EXPECT_LE((field.getData() - loaded.getData()).norm(), 1e-10);
  }
}
TEST(SpatialMatrixIO, HDF5AllSpacesRealAndComplex)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  mesh.getConnectivity().compute(2, 1);
  for (auto shape : {std::pair<size_t, size_t>{2, 3}, {3, 3}})
  {
    auto [r, c] = shape;
    checkPersistence(P0(mesh, r, c));
    checkPersistence(P0g(mesh, r, c));
    checkPersistence(P1(mesh, r, c));
    checkPersistence(H1(std::integral_constant<size_t, 1>{}, mesh, r, c));
    checkPersistence(H1(std::integral_constant<size_t, 2>{}, mesh, r, c));
    checkPersistence(H1(std::integral_constant<size_t, 3>{}, mesh, r, c));
    using Matrix = Math::SpatialMatrix<Complex>;
    checkPersistence(P0<Matrix>(mesh, r, c));
    checkPersistence(P0g<Matrix, LocalMesh>(mesh, r, c));
    checkPersistence(P1<Matrix>(mesh, r, c));
    checkPersistence(H1<1, Matrix>(std::integral_constant<size_t, 1>{}, mesh, r, c));
    checkPersistence(H1<2, Matrix>(std::integral_constant<size_t, 2>{}, mesh, r, c));
    checkPersistence(H1<3, Matrix>(std::integral_constant<size_t, 3>{}, mesh, r, c));
  }
}
TEST(SpatialMatrixIO, MEDITRectangularVertexFields)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  mesh.getConnectivity().compute(2, 1);
  checkVertexFormat(P1(mesh, 2, 3));
  checkVertexFormat(H1(std::integral_constant<size_t, 1>{}, mesh, 2, 3));
  checkVertexFormat(H1(std::integral_constant<size_t, 2>{}, mesh, 2, 3));
  checkVertexFormat(H1(std::integral_constant<size_t, 3>{}, mesh, 2, 3));
}

TEST(SpatialMatrixIO, MFEMUsesScalarNodePermutations)
{
  for (auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Wedge, Polytope::Type::Hexahedron})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const size_t D = Polytope::Traits(geometry).getDimension();
    auto mesh = D == 1 ? LocalMesh::UniformGrid(geometry, {2})
      : D == 2         ? LocalMesh::UniformGrid(geometry, {2, 2})
                       : LocalMesh::UniformGrid(geometry, {2, 2, 2});
    for (size_t d = 1; d <= D; ++d)
      for (size_t lower = 0; lower < d; ++lower)
        mesh.getConnectivity().compute(d, lower);
    checkVertexFormat(P0(mesh, 2, 3), FileFormat::MFEM);
    checkVertexFormat(P1(mesh, 2, 3), FileFormat::MFEM);
    checkVertexFormat(
      H1(std::integral_constant<size_t, 1>{}, mesh, 2, 3), FileFormat::MFEM);
    checkVertexFormat(
      H1(std::integral_constant<size_t, 2>{}, mesh, 2, 3), FileFormat::MFEM);
    checkVertexFormat(
      H1(std::integral_constant<size_t, 3>{}, mesh, 2, 3), FileFormat::MFEM);
  }
}

TEST(SpatialMatrixIO, XDMFShapeAndNonsymmetricEntries)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  for (auto shape : {std::pair<size_t, size_t>{2, 3}, {3, 3}})
  {
    auto [r, c] = shape;
    P1 fes(mesh, r, c);
    GridFunction field(fes);
    Math::SpatialMatrix<Real> value(r, c);
    for (size_t i = 0; i < r; ++i)
      for (size_t j = 0; j < c; ++j)
        value(i, j) = 1 + 7 * i + j;
    field = value;
    const boost::filesystem::path dir = "/tmp/rodin_matrix_xdmf";
    boost::filesystem::remove_all(dir);
    boost::filesystem::create_directories(dir);
    {
      XDMF output(dir / "field");
      output.setMesh(mesh);
      output.add("matrix", field, XDMF::Center::Node);
      output.write(0);
      output.close();
    }
    std::ifstream input((dir / "field.xdmf").string());
    ASSERT_TRUE(input.good());
    const std::string xml((std::istreambuf_iterator<char>(input)), {});
    EXPECT_NE(xml.find(r == 3 ? "AttributeType=\"Tensor\"" : "AttributeType=\"Matrix\""),
      std::string::npos);
    EXPECT_NE(xml.find("Dimensions=\"" + std::to_string(mesh.getVertexCount()) + " " +
                std::to_string(r) + " " + std::to_string(c) + "\""),
      std::string::npos);
    boost::filesystem::remove_all(dir);
  }
}
