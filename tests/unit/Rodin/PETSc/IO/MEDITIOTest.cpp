/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Vertex-value, row-layout and complex-rejection oracles for PETSc MEDIT output.
#include <gtest/gtest.h>
#include <petsc.h>
#include <memory>
#include <sstream>
#include <string>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>
#include <Rodin/PETSc/IO.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::IO;

namespace
{
  /// @brief Real output preserves vertex values; complex output is rejected transactionally.
  template <bool Quadratic, size_t Rank, class Scalar = PetscScalar>
  void checkVertexRows()
  {
    using Range = std::conditional_t<Rank == 0, Scalar,
      std::conditional_t<Rank == 1, Math::SpatialVector<Scalar>,
        Math::SpatialMatrix<Scalar>>>;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.scale(Real(0.5));
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    auto space = [&] {
      if constexpr (Quadratic)
      {
        if constexpr (Rank == 0)
          return H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh);
        else if constexpr (Rank == 1)
          return H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh, 2);
        else
          return H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh, 2, 3);
      }
      else
      {
        if constexpr (Rank == 0)
          return P1<Range>(mesh);
        else if constexpr (Rank == 1)
          return P1<Range>(mesh, 2);
        else
          return P1<Range>(mesh, 2, 3);
      }
    }();
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -2);
      else
        return Scalar(1);
    }();
    auto entry = [factor](const Point& p, size_t row, size_t column) -> Scalar {
      return factor *
        Scalar(1 + 3 * row + column + Real(0.125) * p.x() +
          (Real(0.25) + Real(0.125) * column) * p.y());
    };
    auto value = [&](const Point& p) -> Range {
      if constexpr (Rank == 0)
        return entry(p, 0, 0);
      else if constexpr (Rank == 1)
      {
        Range result(2);
        for (size_t row = 0; row < 2; ++row)
          result[row] = entry(p, row, 0);
        return result;
      }
      else
      {
        Range result(2, 3);
        for (size_t row = 0; row < 2; ++row)
        {
          for (size_t column = 0; column < 3; ++column)
            result(row, column) = entry(p, row, column);
        }
        return result;
      }
    };
    PETSc::Variational::GridFunction field(space);
    field = value;
    using FES = std::decay_t<decltype(space)>;
    using Base = GridFunctionPrinterBase<FileFormat::MEDIT, FES, Vec>;
    std::unique_ptr<Base> printer(
      new GridFunctionPrinter<FileFormat::MEDIT, FES, Vec>(field));
    const size_t components = Rank == 0 ? 1 : Rank == 1 ? 2 : 6;
    for (size_t repeat = 0; repeat < 2; ++repeat)
    {
      std::stringstream stream;
      if constexpr (std::is_same_v<Scalar, Complex>)
      {
        stream << "retained";
        EXPECT_THROW(printer->print(stream), Alert::Exception);
        EXPECT_EQ(stream.str(), "retained");
        continue;
      }
      printer->print(stream);
      std::string keyword;
      size_t version, dimension, vertices, fields;
      ASSERT_TRUE(stream >> keyword >> version);
      EXPECT_EQ(keyword, "MeshVersionFormatted");
      EXPECT_EQ(version, 2u);
      ASSERT_TRUE(stream >> keyword >> dimension);
      EXPECT_EQ(keyword, "Dimension");
      EXPECT_EQ(dimension, 2u);
      ASSERT_TRUE(stream >> keyword >> vertices >> fields);
      EXPECT_EQ(keyword, "SolAtVertices");
      EXPECT_EQ(vertices, mesh.getVertexCount());
      EXPECT_EQ(fields, Rank == 2 ? components : 1u);
      for (size_t component = 0; component < fields; ++component)
      {
        size_t type;
        ASSERT_TRUE(stream >> type);
        EXPECT_EQ(type, Rank == 1 ? 2u : 1u);
      }
      std::string line;
      std::getline(stream, line);
      for (auto vertex = mesh.getVertex(); vertex; ++vertex)
      {
        ASSERT_TRUE(std::getline(stream, line));
        std::stringstream row(line);
        const Point point(*vertex, Math::SpatialPoint{}, vertex->getCoordinates());
        for (size_t component = 0; component < components; ++component)
        {
          Scalar actual;
          ASSERT_TRUE(row >> actual);
          const size_t r = Rank == 2 ? component / 3 : component;
          const size_t c = Rank == 2 ? component % 3 : 0;
          EXPECT_LT(std::abs(actual - entry(point, r, c)), 1e-12);
        }
        row >> std::ws;
        EXPECT_TRUE(row.eof());
      }
      ASSERT_TRUE(stream >> keyword);
      EXPECT_EQ(keyword, "End");
      stream >> std::ws;
      EXPECT_TRUE(stream.eof());
    }
  }
}

TEST(PETScMEDIT, LinearScalar)
{
  checkVertexRows<false, 0>();
}
TEST(PETScMEDIT, LinearVector)
{
  checkVertexRows<false, 1>();
}
TEST(PETScMEDIT, LinearMatrix)
{
  checkVertexRows<false, 2>();
}
TEST(PETScMEDIT, QuadraticScalar)
{
  checkVertexRows<true, 0>();
}
TEST(PETScMEDIT, QuadraticVector)
{
  checkVertexRows<true, 1>();
}
TEST(PETScMEDIT, QuadraticMatrix)
{
  checkVertexRows<true, 2>();
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
