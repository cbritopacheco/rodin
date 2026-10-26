/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Exact matrix-valued H1 basis binding and reference-association checks.
 */
#include <gtest/gtest.h>
#include <type_traits>
#include <utility>
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/QF/QuadratureFormula.h"
#include "Rodin/Variational/H1.h"
#include "../../../../convergence/Convergence.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Checks matrix component bases without using a cached-table oracle.
   *
   * For the scalar basis @f$\phi_a@f$ and logical component unit matrix
   * @f$E_c@f$, the matrix basis is @f$\Phi_{a,c}=\phi_a E_c@f$. Every active
   * matrix entry is compared exactly against direct scalar basis evaluation.
   * Formula/sample transitions and unproven sampling metadata must preserve
   * that identity. Point is an evaluation case, not a refinement hierarchy.
   *
   * @tparam Case Degree constant and real or complex scalar type.
   */
  template <class Case>
  class H1MatrixShapeFunctionTest : public ::testing::Test
  {};

  using H1MatrixShapeCases =
    ::testing::Types<std::pair<std::integral_constant<size_t, 1>, Real>,
      std::pair<std::integral_constant<size_t, 2>, Real>,
      std::pair<std::integral_constant<size_t, 3>, Real>,
      std::pair<std::integral_constant<size_t, 4>, Real>,
      std::pair<std::integral_constant<size_t, 5>, Real>,
      std::pair<std::integral_constant<size_t, 6>, Real>,
      std::pair<std::integral_constant<size_t, 1>, Complex>,
      std::pair<std::integral_constant<size_t, 2>, Complex>,
      std::pair<std::integral_constant<size_t, 3>, Complex>,
      std::pair<std::integral_constant<size_t, 4>, Complex>,
      std::pair<std::integral_constant<size_t, 5>, Complex>,
      std::pair<std::integral_constant<size_t, 6>, Complex>>;
  TYPED_TEST_SUITE(H1MatrixShapeFunctionTest, H1MatrixShapeCases);

  TYPED_TEST(H1MatrixShapeFunctionTest, ReferenceBindingAndValueLifetimes)
  {
    using namespace Geometry;
    using namespace Variational;
    using Order = typename TypeParam::first_type;
    using Scalar = typename TypeParam::second_type;
    using Matrix = Math::SpatialMatrix<Scalar>;
    using G = Polytope::Type;
    for (const auto geometry : {G::Point, G::Segment, G::Triangle, G::Quadrilateral,
           G::Tetrahedron, G::Pyramid, G::Hexahedron, G::Wedge})
    {
      SCOPED_TRACE(int(geometry));
      auto mesh = [&]() {
        if (geometry == G::Point)
          return LocalMesh::Builder().initialize(3).nodes(1).vertex({0, 0, 0}).finalize();
        return Convergence::UniformGrid(geometry).makeMesh(2);
      }();
      for (size_t d = 1; d <= mesh.getDimension(); ++d)
      {
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      }
      const auto cell = *mesh.getCell();
      const auto& first = QF::QuadratureFormula::get(2, geometry);
      const auto& second = QF::QuadratureFormula::get(4, geometry);
      const PolytopeQuadrature firstSamples(cell, first), secondSamples(cell, second);
      // The single borrowed cursor outlives every expression created below.
      IntegrationPoint cursor(firstSamples, 0);
      for (size_t rows = 1; rows <= 3; ++rows)
      {
        for (size_t columns = 1; columns <= 3; ++columns)
        {
          SCOPED_TRACE(::testing::Message() << rows << "x" << columns);
          H1<Order::value, Matrix> space(Order{}, mesh, rows, columns);
          const auto& element =
            space.getFiniteElement(cell.getDimension(), cell.getIndex());
          const auto& scalar = element.getScalarElement();
          const size_t components = rows * columns;
          const auto& mutableFormula =
            QF::QuadratureFormula::get<QF::GaussLegendre>(geometry, size_t{2});
          const PolytopeQuadrature savedSamples(cell, mutableFormula);
          TrialFunction original(space);
          const auto check = [&](auto& unknown) {
            auto& shape = unknown.getDerived();
            shape.setIntegrationPoint(cursor);
            ASSERT_EQ(shape.getDOFs(cell), scalar.getCount() * components);
            const auto& coordinates = cursor.getPoint().getReferenceCoordinates();
            for (size_t a = 0; a < scalar.getCount(); ++a)
            {
              const Scalar value = scalar.getBasis(a)(coordinates);
              for (size_t c = 0; c < components; ++c)
              {
                const auto& actual = shape.getBasis(a * components + c);
                ASSERT_EQ(actual.rows(), rows);
                ASSERT_EQ(actual.cols(), columns);
                for (size_t entry = 0; entry < components; ++entry)
                {
                  ASSERT_EQ(actual(entry / columns, entry % columns),
                    entry == c ? value : Scalar(0))
                    << "basis=" << a << " component=" << c << " entry=" << entry;
                }
              }
            }
          };
          for (const auto* samples : {&firstSamples, &secondSamples, &firstSamples})
          {
            for (const size_t qp : {size_t(0), samples->getSize() - 1, size_t(0)})
            {
              cursor = IntegrationPoint(*samples, qp);
              check(original);
              check(original);
            }
          }
          TrialFunction copied(original);
          TrialFunction moved(std::move(copied));
          cursor = IntegrationPoint(secondSamples, secondSamples.getSize() - 1);
          check(moved);
          cursor = IntegrationPoint(firstSamples, 0);
          check(original);

          // Formula/index provenance alone must not select the origin's sample.
          const Point direct(cell, second.getPoint(second.getSize() - 1));
          cursor = IntegrationPoint(direct, &first, 0);
          ASSERT_EQ(cursor.getReferenceQuadratureFormula(), nullptr);
          check(original);
          cursor = IntegrationPoint(direct);
          check(original);

          cursor = IntegrationPoint(savedSamples, 0);
          check(original);
          ASSERT_EQ(cursor.getReferenceQuadratureFormula(), &mutableFormula);
          check(original);
          const auto& replacement =
            QF::QuadratureFormula::get<QF::GaussLegendre>(geometry, size_t{4});
          const PolytopeQuadrature rebuiltSamples(cell, replacement);
          cursor = IntegrationPoint(rebuiltSamples, rebuiltSamples.getSize() - 1);
          ASSERT_EQ(cursor.getReferenceQuadratureFormula(), &replacement);
          check(original);
          cursor = IntegrationPoint(firstSamples, 0);
          check(original);
        }
      }
    }
  }

  TYPED_TEST(H1MatrixShapeFunctionTest, BorrowsReadyTablesWithoutFillingMisses)
  {
    using namespace Geometry;
    using namespace Variational;
    using Order = typename TypeParam::first_type;
    using Scalar = typename TypeParam::second_type;
    using Matrix = Math::SpatialMatrix<Scalar>;
    using G = Polytope::Type;
    // One more identity than the current eight-entry reference-table store.
    // This checks the present eviction protocol, not a public capacity promise.
    constexpr size_t EvictionWorkingSet = 9;
    for (const auto geometry : {G::Point, G::Segment, G::Triangle, G::Quadrilateral,
           G::Tetrahedron, G::Pyramid, G::Hexahedron, G::Wedge})
    {
      SCOPED_TRACE(int(geometry));
      auto mesh = [&]() {
        if (geometry == G::Point)
          return LocalMesh::Builder().initialize(3).nodes(1).vertex({0, 0, 0}).finalize();
        return Convergence::UniformGrid(geometry).makeMesh(2);
      }();
      for (size_t d = 1; d <= mesh.getDimension(); ++d)
      {
        for (size_t lower = 0; lower < d; ++lower)
          mesh.getConnectivity().compute(d, lower);
      }
      const auto cell = *mesh.getCell();
      H1<Order::value, Matrix> space(Order{}, mesh, 2, 3);
      const auto& scalar =
        space.getFiniteElement(cell.getDimension(), cell.getIndex()).getScalarElement();
      const auto& formula =
        QF::QuadratureFormula::get<QF::GaussLegendre>(geometry, size_t{4});
      const auto& identity = formula;
      const PolytopeQuadrature samples(cell, formula);
      IntegrationPoint cursor(samples, 0);
      const auto check = [&](auto& unknown) {
        auto& shape = unknown.getDerived();
        shape.setIntegrationPoint(cursor);
        const auto& coordinates = cursor.getPoint().getReferenceCoordinates();
        for (size_t a = 0; a < scalar.getCount(); ++a)
        {
          const Scalar value = scalar.getBasis(a)(coordinates);
          for (size_t component = 0; component < 6; ++component)
          {
            const auto& actual = shape.getBasis(6 * a + component);
            ASSERT_EQ(actual.rows(), 2);
            ASSERT_EQ(actual.cols(), 3);
            for (size_t entry = 0; entry < 6; ++entry)
            {
              ASSERT_EQ(
                actual(entry / 3, entry % 3), entry == component ? value : Scalar(0));
            }
          }
        }
      };

      ASSERT_FALSE(scalar.findTabulation(identity));
      TrialFunction first(space);
      check(first);
      // A pointwise value request must not build all samples and derivatives.
      ASSERT_FALSE(scalar.findTabulation(identity));

      const auto* ready = &scalar.getTabulation(formula);
      ASSERT_EQ(&scalar.findTabulation(identity)->get(), ready);
      TrialFunction warmed(space);
      for (const size_t qp : {samples.getSize() - 1, size_t(0)})
      {
        cursor = IntegrationPoint(samples, qp);
        check(warmed);
        ASSERT_TRUE(scalar.findTabulation(identity));
        ASSERT_EQ(&scalar.findTabulation(identity)->get(), ready);
      }

      // Existing tables still cannot authorize unrelated coordinates or index.
      const Point direct(cell, formula.getPoint(formula.getSize() - 1));
      cursor = IntegrationPoint(direct, &formula, 0);
      ASSERT_EQ(cursor.getReferenceQuadratureFormula(), nullptr);
      check(warmed);

      for (size_t i = 0; i < EvictionWorkingSet; ++i)
      {
        const auto& other =
          QF::QuadratureFormula::get<QF::GaussLegendre>(geometry, i + 5);
        scalar.getTabulation(other);
      }
      ASSERT_FALSE(scalar.findTabulation(identity));
      TrialFunction evicted(space);
      cursor = IntegrationPoint(samples, samples.getSize() - 1);
      check(evicted);
      ASSERT_FALSE(scalar.findTabulation(identity));
      // A previously bound expression must not retain a borrow across eviction.
      check(warmed);
      ASSERT_FALSE(scalar.findTabulation(identity));
    }
  }
}
