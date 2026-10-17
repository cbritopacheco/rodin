/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2024.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "InterfaceTestUtils.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief Average of a constant function across any face should equal the constant.
  TEST(Rodin_Variational_Average, ConstantFunction_AverageEqualsConstant)
  {
    // Average of a constant function across any face should equal the constant
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    RealFunction f(42.0);
    auto avg_f = Average(f);

    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(avg_f.getValue(p), 42.0, 1e-10);
  }

  /// @brief Average of a constant vector function should equal the same vector.
  TEST(Rodin_Variational_Average, ConstantVectorFunction_AverageEqualsConstant)
  {
    // Average of a constant vector function should equal the same vector
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    VectorFunction f{3.0, 7.0};
    auto avg_f = Average(f);

    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    auto result = avg_f.getValue(p);
    EXPECT_NEAR(result(0), 3.0, 1e-10);
    EXPECT_NEAR(result(1), 7.0, 1e-10);
  }

  /// @brief Average of a continuous P1 grid function should equal the function value.
  TEST(Rodin_Variational_Average, GridFunction_ContinuousFunction_AverageEqualsValue)
  {
    // Average of a continuous P1 grid function should equal the function value
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    GridFunction gf(fes);
    RealFunction f(5.0);
    gf.project(f);
    auto avg_gf = Average(gf);

    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(avg_gf.getValue(p), 5.0, 1e-10);
  }

  /// @brief Verifies copy construction for variational average by checking tolerance-based numerical results, true predicates, copy semantics.
  TEST(Rodin_Variational_Average, CopyConstruction)
  {
    RealFunction f(9.0);
    auto avg_f = Average(f);
    auto avg_copy = avg_f;

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(avg_copy.getValue(p), 9.0, 1e-10);
  }

  /// @brief Verifies move construction for variational average by checking tolerance-based numerical results, true predicates, move semantics.
  TEST(Rodin_Variational_Average, MoveConstruction)
  {
    RealFunction f(9.0);
    auto avg_f = Average(f);
    auto avg_moved = std::move(avg_f);

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(avg_moved.getValue(p), 9.0, 1e-10);
  }

  /// @brief Verifies polymorphic copy for variational average by checking tolerance-based numerical results, true predicates, copy semantics.
  TEST(Rodin_Variational_Average, PolymorphicCopy)
  {
    RealFunction f(11.0);
    auto avg_f = Average(f);
    std::unique_ptr<decltype(avg_f)> copy(avg_f.copy());

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(copy->getValue(p), 11.0, 1e-10);
  }

  /// @brief For a constant function, Jump=0 and Average=constant.
  TEST(Rodin_Variational_Average, JumpAndAverage_Consistency)
  {
    // For a constant function, Jump=0 and Average=constant
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    RealFunction f(10.0);
    auto avg_f = Average(f);
    auto jump_f = Jump(f);

    auto [found, p] = findInteriorFacePoint(mesh);
    ASSERT_TRUE(found) << "No interior face found in mesh";

    EXPECT_NEAR(jump_f.getValue(p), 0.0, 1e-10);
    EXPECT_NEAR(avg_f.getValue(p), 10.0, 1e-10);
  }
}

/// @brief Quadrature-associated interface values retain distinct traces and P0 basis values.
TEST(Rodin_Variational_Average, DiscontinuousIntegrationPointTracesAllGeometries)
{
  const auto check = []<class Scalar>() {
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Wedge, Polytope::Type::Hexahedron})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const auto dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {3})
        : dimension == 2         ? LocalMesh::UniformGrid(geometry, {3, 3})
                                 : LocalMesh::UniformGrid(geometry, {3, 3, 3});
      mesh.getConnectivity().compute(dimension - 1, dimension);
      P0<Scalar> fes(mesh);
      GridFunction field(fes);
      for (Index cell = 0; cell < mesh.getCellCount(); ++cell)
      {
        if constexpr (std::is_same_v<Scalar, Complex>)
          field[cell] = Complex(Real(2 + 2 * cell), Real(3 + cell));
        else
          field[cell] = Real(2 + 2 * cell);
      }
      auto average = Average(field);
      auto jump = Jump(field);
      TrialFunction trial(fes);
      auto averageBasis = Average(trial);
      auto jumpBasis = Jump(trial);
      bool checked = false;
      for (auto face = mesh.getFace(); face; ++face)
      {
        const auto& incidence = mesh.getConnectivity().getIncidence(
          {dimension - 1, dimension}, face->getIndex());
        if (incidence.size() != 2)
          continue;
        const Index first = *incidence.begin();
        const Index second = *std::next(incidence.begin());
        const auto& quadrature =
          QF::PolytopeQuadratureFormula::get(2, face->getGeometry());
        const auto& points = face->getQuadrature(quadrature);
        for (size_t q = 0; q < points.getSize(); ++q)
        {
          const IntegrationPoint point(points.getPoint(q), &quadrature, q);
          const Scalar firstValue = field[first];
          const Scalar secondValue = field[second];
          EXPECT_EQ(average.getValue(point), Real(0.5) * (firstValue + secondValue));
          EXPECT_EQ(jump.getValue(point), firstValue - secondValue);
          EXPECT_EQ(average.getValue(point.getPoint()), average.getValue(point));
          EXPECT_EQ(jump.getValue(point.getPoint()), jump.getValue(point));
          EXPECT_EQ(&averageBasis.setIntegrationPoint(point), &averageBasis);
          EXPECT_EQ(&jumpBasis.setIntegrationPoint(point), &jumpBasis);
          EXPECT_EQ(&averageBasis.getIntegrationPoint(), &point);
          EXPECT_EQ(&jumpBasis.getIntegrationPoint(), &point);
          EXPECT_EQ(averageBasis.getBasis(0), Scalar(1));
          EXPECT_EQ(jumpBasis.getBasis(0), Scalar(0));
        }
        checked = true;
        break;
      }
      EXPECT_TRUE(checked);
    }
  };
  check.template operator()<Real>();
  check.template operator()<Complex>();
}
