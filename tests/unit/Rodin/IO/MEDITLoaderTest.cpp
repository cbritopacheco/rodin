/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Analytic vertex-field loading and transactional MEDIT input validation.
#include <gtest/gtest.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/IO.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::IO;

namespace
{
  /// @brief Independent affine entries distinguish every component and coordinate axis.
  template <class Scalar, size_t Rank, bool Quadratic>
  void checkLoaderContracts()
  {
    using Range = std::conditional_t<Rank == 0, Scalar,
      std::conditional_t<Rank == 1, Math::SpatialVector<Scalar>,
        Math::SpatialMatrix<Scalar>>>;
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -3);
      else
        return Scalar(2);
    }();
    for (const auto geometry :
      {Polytope::Type::Triangle, Polytope::Type::Point, Polytope::Type::Segment,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Pyramid, Polytope::Type::Wedge, Polytope::Type::Hexahedron})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = [&] {
        if (dimension == 0)
          return LocalMesh::Builder().initialize(1).nodes(1).vertex({0}).finalize();
        return LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
      }();
      for (size_t d = 1; d <= dimension; ++d)
      {
        mesh.getConnectivity().compute(dimension, d);
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      }
      const size_t vdim = Rank == 0 ? 1 : Rank == 1 ? mesh.getSpaceDimension() : 6;
      auto space = [&] {
        if constexpr (Quadratic)
        {
          if constexpr (Rank == 0)
            return H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh);
          else if constexpr (Rank == 1)
            return H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh, vdim);
          else
            return H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh, 2, 3);
        }
        else
        {
          if constexpr (Rank == 0)
            return P1<Range>(mesh);
          else if constexpr (Rank == 1)
            return P1<Range>(mesh, vdim);
          else
            return P1<Range>(mesh, 2, 3);
        }
      }();
      const size_t spaceDimension = mesh.getSpaceDimension();
      auto entry = [factor, spaceDimension](
                     const Point& point, size_t component) -> Scalar {
        const size_t row = Rank == 2 ? component / 3 : component;
        const size_t column = Rank == 2 ? component % 3 : 0;
        Real value = 1 + 3 * row + column;
        for (size_t axis = 0; axis < spaceDimension; ++axis)
          value += Real((axis + 1) * (column + 1)) * point(axis);
        return factor * Scalar(value);
      };
      auto header = [&](bool inlineValues) {
        std::ostringstream stream;
        stream << "\nMeshVersionFormatted" << (inlineValues ? " 2\n" : "\n2\n\n")
               << "Dimension" << (inlineValues ? " " : "\n") << mesh.getSpaceDimension()
               << "\n\nSolAtVertices\n"
               << mesh.getVertexCount() << '\n';
        if constexpr (Rank == 2)
          stream << "6 1 1 1 1 1 1\n";
        else
          stream << "1 " << (vdim == 1 ? 1 : 2) << '\n';
        return stream.str();
      };
      std::ostringstream payload;
      payload << std::setprecision(std::numeric_limits<Real>::max_digits10);
      for (auto vertex = mesh.getVertex(); vertex; ++vertex)
      {
        const Point point(*vertex, Math::SpatialPoint{}, vertex->getCoordinates());
        for (size_t component = 0; component < vdim; ++component)
          payload << entry(point, component) << ' ';
        payload << '\n';
      }
      GridFunction field(space);
      using FES = std::decay_t<decltype(space)>;
      GridFunctionLoader<FileFormat::MEDIT, FES, Math::Vector<Scalar>> loader(field);
      auto verify = [&] {
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const auto& formula =
            QF::PolytopeQuadratureFormula::get(3, cell->getGeometry());
          for (size_t sample = 0; sample < formula.getSize(); ++sample)
          {
            const Point point(*cell, formula.getPoint(sample));
            const auto value = field.getValue(point);
            for (size_t component = 0; component < vdim; ++component)
            {
              Scalar actual;
              if constexpr (Rank == 0)
                actual = value;
              else if constexpr (Rank == 1)
                actual = value[component];
              else
                actual = value(component / 3, component % 3);
              EXPECT_LT(std::abs(actual - entry(point, component)), 1e-10);
            }
          }
        }
      };
      for (const bool inlineValues : {false, true})
      {
        std::istringstream input(header(inlineValues) + payload.str() + "\nEnd\n");
        loader.load(input);
        verify();
      }
      if (geometry != Polytope::Type::Triangle)
        continue;
      const std::string valid = header(true) + payload.str() + "\nEnd\n";
      auto replace = [](std::string value, const std::string& before,
                       const std::string& after) {
        const auto position = value.find(before);
        EXPECT_NE(position, std::string::npos);
        if (position != std::string::npos)
          value.replace(position, before.size(), after);
        return value;
      };
      std::vector<std::string> invalid;
      if constexpr (Rank == 2)
      {
        invalid = {replace(valid, "Formatted 2", "Formatted 7"),
          replace(valid, "Dimension 2", "Dimension 3"),
          replace(valid, "Vertices\n4", "Vertices\n3"),
          replace(valid, "6 1 1 1 1 1 1", "5 1 1 1 1 1"),
          replace(valid, "6 1 1 1 1 1 1", "6 1 1 2 1 1 1"), header(true) + "1 2\nEnd\n"};
      }
      else
      {
        const std::string type = vdim == 1 ? "1 1\n" : "1 2\n";
        invalid = {"MeshVersionFormatted\ninvalid\n",
          "MeshVersionFormatted 2\nDimension\ninvalid\n",
          replace(valid, "Dimension 2", "Dimension 3"),
          replace(valid, "SolAtVertices", "SolAtEdges"),
          replace(valid, "Vertices\n4", "Vertices\ninvalid"),
          replace(valid, "Vertices\n4", "Vertices\n3"), replace(valid, type, "invalid\n"),
          replace(valid, type, "2 1\n"),
          replace(valid, type, vdim == 1 ? "1 2\n" : "1 1\n"),
          replace(valid, type, "1 3\n"), header(true) + "1\nEnd\n"};
      }
      for (const auto& text : invalid)
      {
        SCOPED_TRACE(text);
        field.getData().setConstant(Scalar(-17));
        const Math::Vector<Scalar> before = field.getData();
        std::string firstMessage;
        for (size_t repeat = 0; repeat < 2; ++repeat)
        {
          bool rejected = false;
          try
          {
            std::istringstream input(text);
            loader.load(input);
          }
          catch (const Alert::Exception& error)
          {
            rejected = true;
            const std::string message = error.what();
            EXPECT_NE(message.find("member function"), std::string::npos);
            if (repeat == 0)
              firstMessage = message;
            else
              EXPECT_EQ(message, firstMessage);
          }
          EXPECT_TRUE(rejected);
          EXPECT_EQ((field.getData() - before).norm(), 0);
        }
        std::istringstream recovery(valid);
        loader.load(recovery);
        verify();
      }
    }
  }
}

TEST(MEDITLoader, P1RealScalar)
{
  checkLoaderContracts<Real, 0, false>();
}
TEST(MEDITLoader, P1ComplexScalar)
{
  checkLoaderContracts<Complex, 0, false>();
}
TEST(MEDITLoader, P1RealVector)
{
  checkLoaderContracts<Real, 1, false>();
}
TEST(MEDITLoader, P1ComplexVector)
{
  checkLoaderContracts<Complex, 1, false>();
}
TEST(MEDITLoader, H1RealScalar)
{
  checkLoaderContracts<Real, 0, true>();
}
TEST(MEDITLoader, H1ComplexScalar)
{
  checkLoaderContracts<Complex, 0, true>();
}
TEST(MEDITLoader, H1RealVector)
{
  checkLoaderContracts<Real, 1, true>();
}
TEST(MEDITLoader, H1ComplexVector)
{
  checkLoaderContracts<Complex, 1, true>();
}
TEST(MEDITLoader, P1RealMatrix)
{
  checkLoaderContracts<Real, 2, false>();
}
TEST(MEDITLoader, P1ComplexMatrix)
{
  checkLoaderContracts<Complex, 2, false>();
}
TEST(MEDITLoader, H1RealMatrix)
{
  checkLoaderContracts<Real, 2, true>();
}
TEST(MEDITLoader, H1ComplexMatrix)
{
  checkLoaderContracts<Complex, 2, true>();
}

/// @brief Nonspatial vector components are declared as independent real fields.
TEST(MEDITInterchange, ExplicitComponentLayout)
{
  const auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  P1<Math::SpatialVector<Real>> space(mesh, 3);
  GridFunction field(space);
  field = VectorFunction(1., 2., 3.);
  std::ostringstream output;
  GridFunctionPrinter<FileFormat::MEDIT, decltype(space), Math::Vector<Real>>(field)
    .print(output);
  std::istringstream external(output.str());
  std::string keyword;
  size_t version = 0, dimension = 0, vertices = 0, fields = 0;
  ASSERT_TRUE(bool(external >> keyword >> version));
  EXPECT_EQ(keyword, "MeshVersionFormatted");
  ASSERT_TRUE(bool(external >> keyword >> dimension));
  EXPECT_EQ(dimension, 2u);
  ASSERT_TRUE(bool(external >> keyword >> vertices >> fields));
  EXPECT_EQ(keyword, "SolAtVertices");
  EXPECT_EQ(fields, 3u);
  for (size_t fieldIndex = 0; fieldIndex < fields; ++fieldIndex)
  {
    size_t type = 0;
    ASSERT_TRUE(bool(external >> type));
    EXPECT_EQ(type, 1u);
  }
  for (size_t vertex = 0; vertex < vertices; ++vertex)
  {
    for (size_t component = 0; component < fields; ++component)
    {
      Real value = 0;
      ASSERT_TRUE(bool(external >> value));
      EXPECT_EQ(value, Real(component + 1));
    }
  }
  GridFunction restored(space);
  std::istringstream input(output.str());
  GridFunctionLoader<FileFormat::MEDIT, decltype(space), Math::Vector<Real>>(restored)
    .load(input);
  EXPECT_EQ(restored.getData(), field.getData());
}

/// @brief Unsupported complex output is rejected before any format bytes are written.
TEST(MEDITInterchange, RejectsComplexOutput)
{
  const auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  P1<Complex> space(mesh);
  GridFunction field(space);
  field = ComplexFunction(Complex(1, 2));
  std::ostringstream output;
  output << "retained";
  EXPECT_THROW(
    (GridFunctionPrinter<FileFormat::MEDIT, decltype(space), Math::Vector<Complex>>(field)
        .print(output)),
    Alert::Exception);
  EXPECT_EQ(output.str(), "retained");
}
