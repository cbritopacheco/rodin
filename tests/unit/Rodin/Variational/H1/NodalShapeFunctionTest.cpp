/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Exact pointwise nodal identities for scalar, vector and matrix H1 bases.
 */
#include <gtest/gtest.h>
#include <cmath>
#include "Rodin/Variational/H1.h"

namespace Rodin::Tests::Unit
{
  template <class Scalar>
  class H1NodalShapeFunctionTest : public ::testing::Test
  {};

  using H1NodalScalars = ::testing::Types<Real, Complex>;
  TYPED_TEST_SUITE(H1NodalShapeFunctionTest, H1NodalScalars);

  /**
   * @brief Checks exact cardinality against node indices rather than a modal oracle.
   * @tparam Degree Polynomial degree.
   * @tparam Scalar Scalar coefficient type.
   */
  template <size_t Degree, class Scalar>
  void checkNodalShapeFunctions()
  {
    using namespace Geometry;
    using namespace Variational;
    using G = Polytope::Type;
    for (const auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
           G::Pyramid, G::Hexahedron, G::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      const auto counts = Array<size_t>::Constant(dimension, 2);
      auto mesh = LocalMesh::UniformGrid(geometry, counts);
      for (size_t d = 1; d <= dimension; ++d)
      {
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      }
      const auto cell = *mesh.getCell();
      const H1Element<Degree, Scalar> scalar(geometry);
      H1<Degree, Math::SpatialVector<Scalar>> vectorSpace(
        std::integral_constant<size_t, Degree>{}, mesh, 2);
      H1<Degree, Math::SpatialMatrix<Scalar>> matrixSpace(
        std::integral_constant<size_t, Degree>{}, mesh, 2, 2);
      // Explicit scalar type also exercises complex-valued pointwise bindings.
      H1<Degree, Scalar> typedScalarSpace(std::integral_constant<size_t, Degree>{}, mesh);
      TrialFunction scalarTrial(typedScalarSpace);
      TrialFunction vectorTrial(vectorSpace);
      TrialFunction matrixTrial(matrixSpace);
      for (size_t node = 0; node < scalar.getCount(); ++node)
      {
        SCOPED_TRACE(node);
        const Point point(cell, scalar.getNode(node));
        const IntegrationPoint sample(point);
        scalarTrial.setIntegrationPoint(sample);
        vectorTrial.setIntegrationPoint(sample);
        matrixTrial.setIntegrationPoint(sample);
        for (size_t basis = 0; basis < scalar.getCount(); ++basis)
        {
          const Scalar expected(basis == node);
          EXPECT_EQ(scalarTrial.getBasis(basis), expected);
          for (size_t component = 0; component < 2; ++component)
          {
            const auto& value = vectorTrial.getBasis(2 * basis + component);
            for (size_t entry = 0; entry < 2; ++entry)
              EXPECT_EQ(value(entry), entry == component ? expected : Scalar(0));
          }
          for (size_t component = 0; component < 4; ++component)
          {
            const auto& value = matrixTrial.getBasis(4 * basis + component);
            for (size_t entry = 0; entry < 4; ++entry)
            {
              EXPECT_EQ(
                value(entry / 2, entry % 2), entry == component ? expected : Scalar(0));
            }
          }
        }
      }
    }
  }

  /// @brief Node lookup uses exact reference coordinates and ignores only trailing padding.
  TEST(H1NodalShapeFunction, ExactReferenceNodeLookup)
  {
    using namespace Geometry;
    using namespace Variational;
    const H1Element<2, Real> element(Polytope::Type::Triangle);
    const auto& node = element.getNode(0);
    ASSERT_EQ(element.getNodeIndex(node), Optional<size_t>(0));
    auto nearby = node;
    nearby.x() = std::nextafter(node.x(), node.x() + Real(1));
    EXPECT_FALSE(element.getNodeIndex(nearby));
    EXPECT_FALSE(element.getNodeIndex(Math::SpatialPoint{Real(0)}));
    EXPECT_EQ(element.getNodeIndex(Math::SpatialPoint{{node.x(), node.y(), 3}}),
      Optional<size_t>(0));
    const H1Element<2, Real> point(Polytope::Type::Point);
    EXPECT_EQ(point.getNodeIndex(Math::SpatialPoint{{3, 2, 1}}), Optional<size_t>(0));
  }

  /// @brief Quadratic pointwise bases retain exact structural zeros and unit nodal values.
  TYPED_TEST(H1NodalShapeFunctionTest, QuadraticCardinality)
  {
    checkNodalShapeFunctions<2, TypeParam>();
  }

  /// @brief Cubic pointwise bases retain the same cardinality across reference geometries.
  TYPED_TEST(H1NodalShapeFunctionTest, CubicCardinality)
  {
    checkNodalShapeFunctions<3, TypeParam>();
  }
}
