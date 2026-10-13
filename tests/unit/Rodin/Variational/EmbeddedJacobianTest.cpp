/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/**
 * @brief Verifies ambient Jacobian columns against analytic tangent projectors.
 *
 * On an affine embedded cell with tangent projector P, the field
 * @f$u_c=(c+1)x_{s-1}@f$ has derivative
 * @f$J(u)_{cj}=(c+1)P_{s-1,j}@f$. Complex cases have a nonreal multiplier.
 * Grid-function evaluation at direct, unindexed and indexed points is checked
 * independently of shape-function reconstruction. Matrix dimensions and every
 * ambient derivative are checked for P1 and H1 degrees one through six.
 * Constant P0g fields additionally require correctly shaped zero Jacobians
 * and scalar gradients, with no mutation of their borrowed operands.
 */
TEST(EmbeddedJacobian, AmbientColumnsAndTangentAffineFields)
{
  const auto checkMesh = [](const LocalMesh& mesh, const Math::SpatialMatrix<Real>& P) {
    const size_t dim = mesh.getDimension();
    const size_t ambient = mesh.getSpaceDimension();
    const auto checkSpace = [&]<class Space>(const Space& space, bool constant = false) {
      using Scalar = typename FormLanguage::Traits<Space>::ScalarType;
      const size_t components = space.getVectorDimension();
      const Scalar multiplier = []() {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(1, Real(0.5));
        else
          return Scalar(1);
      }();
      GridFunction field(space);
      field = VectorFunction(components, [=](const Point& point) {
        Math::SpatialVector<Scalar> value(components);
        for (size_t c = 0; c < components; ++c)
        {
          value(c) = multiplier * Real(c + 1) *
            (constant ? Real(1) : point.getPhysicalCoordinates()(ambient - 1));
        }
        return value;
      });
      TrialFunction u(space);
      auto shape = Jacobian(u);
      auto jacobian = Jacobian(field);
      // Affine reproduction differs only by interpolation and summation rounding.
      constexpr Real AffineTolerance = 1e-11;
      const auto checkValue = [&](const auto& value) {
        ASSERT_EQ(value.rows(), components);
        ASSERT_EQ(value.cols(), ambient);
        for (size_t c = 0; c < components; ++c)
        {
          for (size_t j = 0; j < ambient; ++j)
          {
            const Scalar expected =
              constant ? Scalar(0) : multiplier * Real(c + 1) * P(ambient - 1, j);
            EXPECT_LE(std::abs(value(c, j) - expected), AffineTolerance);
          }
        }
      };
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto& formula = QF::PolytopeQuadratureFormula::get(4, cell->getGeometry());
        const PolytopeQuadrature quadrature(*cell, formula);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& point = quadrature.getPoint(qp);
          checkValue(jacobian.getValue(point));
          checkValue(jacobian.getValue(IntegrationPoint(point)));
          const IntegrationPoint ip(point, &formula, qp);
          checkValue(jacobian.getValue(ip));
          shape.setIntegrationPoint(ip);
          Math::SpatialMatrix<Scalar> reconstructed(components, ambient);
          reconstructed.setZero();
          for (size_t local = 0; local < u.getDOFs(*cell); ++local)
          {
            const auto basis = shape.getBasis(local);
            ASSERT_EQ(basis.rows(), components);
            ASSERT_EQ(basis.cols(), ambient);
            reconstructed +=
              field[space.getGlobalIndex({dim, cell->getIndex()}, local)] * basis;
          }
          checkValue(reconstructed);
        }
      }
    };
    const auto checkScalar = [&]<class Scalar>() {
      for (size_t components = 1; components <= 3; ++components)
      {
        SCOPED_TRACE(components);
        using Range = Math::SpatialVector<Scalar>;
        checkSpace(P0g<Range, LocalMesh>(mesh, components), true);
        checkSpace(P1<Range>(mesh, components));
        checkSpace(H1<1, Range>(std::integral_constant<size_t, 1>{}, mesh, components));
        checkSpace(H1<2, Range>(std::integral_constant<size_t, 2>{}, mesh, components));
        checkSpace(H1<3, Range>(std::integral_constant<size_t, 3>{}, mesh, components));
        checkSpace(H1<4, Range>(std::integral_constant<size_t, 4>{}, mesh, components));
        checkSpace(H1<5, Range>(std::integral_constant<size_t, 5>{}, mesh, components));
        checkSpace(H1<6, Range>(std::integral_constant<size_t, 6>{}, mesh, components));
      }
      P0g<Scalar, LocalMesh> scalarSpace(mesh);
      GridFunction scalarField(scalarSpace);
      scalarField = Scalar(1);
      TrialFunction scalarTrial(scalarSpace);
      auto gradient = Grad(scalarField);
      auto shapeGradient = Grad(scalarTrial);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto& formula = QF::PolytopeQuadratureFormula::get(4, cell->getGeometry());
        const PolytopeQuadrature quadrature(*cell, formula);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& point = quadrature.getPoint(qp);
          const IntegrationPoint ip(point, &formula, qp);
          const auto pointValue = gradient.getValue(point);
          const auto indexedValue = gradient.getValue(ip);
          ASSERT_EQ(pointValue.size(), ambient);
          ASSERT_EQ(indexedValue.size(), ambient);
          EXPECT_EQ(pointValue.norm(), 0);
          EXPECT_EQ(indexedValue.norm(), 0);
          shapeGradient.setIntegrationPoint(ip);
          const auto basis = shapeGradient.getBasis(0);
          ASSERT_EQ(basis.size(), ambient);
          EXPECT_EQ(basis.norm(), 0);
        }
      }
    };
    checkScalar.template operator()<Real>();
    checkScalar.template operator()<Complex>();
  };
  for (const size_t ambient : {size_t(2), size_t(3)})
  {
    auto builder = LocalMesh::Builder();
    builder.initialize(ambient).nodes(2);
    Math::SpatialPoint zero(ambient), one(ambient);
    zero.setZero();
    one.setConstant(1);
    builder.vertex(zero).vertex(one);
    auto mesh = builder.polytope(Polytope::Type::Segment, {0, 1}).finalize();
    Math::SpatialMatrix<Real> projector(ambient, ambient);
    projector.setConstant(Real(1) / ambient);
    checkMesh(mesh, projector);
  }
  for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral})
  {
    auto builder = LocalMesh::Builder();
    builder.initialize(3).nodes(4);
    builder.vertex({0, 0, 0}).vertex({1, 0, 1}).vertex({0, 1, 1}).vertex({1, 1, 2});
    if (geometry == Polytope::Type::Triangle)
      builder.polytope(geometry, {0, 1, 2});
    else
      builder.polytope(geometry, {0, 1, 3, 2});
    auto mesh = builder.finalize();
    mesh.getConnectivity().compute(2, 1);
    Math::SpatialMatrix<Real> projector(3, 3);
    const Math::SpatialVector<Real> normal{-1, -1, 1};
    for (size_t i = 0; i < 3; ++i)
    {
      for (size_t j = 0; j < 3; ++j)
        projector(i, j) = Real(i == j) - normal(i) * normal(j) / 3;
    }
    checkMesh(mesh, projector);
  }
  auto pointBuilder = LocalMesh::Builder();
  pointBuilder.initialize(3).nodes(1).vertex({0, 0, 0});
  Math::SpatialMatrix<Real> zero(3, 3);
  zero.setZero();
  checkMesh(pointBuilder.finalize(), zero);
}
