/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Exact matrix-basis lifecycle checks across finite-element families.
 */
#include <gtest/gtest.h>
#include <type_traits>
#include <utility>
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/P0.h"
#include "Rodin/Variational/P0g.h"
#include "Rodin/Variational/P1.h"
#include "../../../convergence/Convergence.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Checks basis support after changing an existing matrix-space definition.
   *
   * The formula identity, sample index, geometry and space address are held
   * fixed while matrix rows and columns change. In particular, the transition
   * @f$2\times3\to3\times2@f$ retains the local basis count but changes every
   * returned matrix shape. For scalar bases @f$\phi_a@f$, every logical entry
   * of @f$\Phi_{a,c}=\phi_a E_c@f$ is checked against direct scalar evaluation.
   * No tolerance or coordinate-based cache identity is used.
   *
   * @tparam Case Polynomial degree constant and matrix-valued space type.
   */
  template <class Case>
  class MatrixShapeLifecycleTest : public ::testing::Test
  {};

  /// @brief Polynomial degree and matrix-space type used by one typed case.
  template <size_t K, class FES>
  using MatrixShapeCase = std::pair<std::integral_constant<size_t, K>, FES>;

  using RealMatrix = Math::SpatialMatrix<Real>;
  using ComplexMatrix = Math::SpatialMatrix<Complex>;
  using MatrixShapeCases =
    ::testing::Types<MatrixShapeCase<0, Variational::P0<RealMatrix>>,
      MatrixShapeCase<0, Variational::P0g<RealMatrix, Geometry::LocalMesh>>,
      MatrixShapeCase<1, Variational::P1<RealMatrix>>,
      MatrixShapeCase<1, Variational::H1<1, RealMatrix>>,
      MatrixShapeCase<2, Variational::H1<2, RealMatrix>>,
      MatrixShapeCase<3, Variational::H1<3, RealMatrix>>,
      MatrixShapeCase<4, Variational::H1<4, RealMatrix>>,
      MatrixShapeCase<5, Variational::H1<5, RealMatrix>>,
      MatrixShapeCase<6, Variational::H1<6, RealMatrix>>,
      MatrixShapeCase<0, Variational::P0<ComplexMatrix>>,
      MatrixShapeCase<0, Variational::P0g<ComplexMatrix, Geometry::LocalMesh>>,
      MatrixShapeCase<1, Variational::P1<ComplexMatrix>>,
      MatrixShapeCase<1, Variational::H1<1, ComplexMatrix>>,
      MatrixShapeCase<2, Variational::H1<2, ComplexMatrix>>,
      MatrixShapeCase<3, Variational::H1<3, ComplexMatrix>>,
      MatrixShapeCase<4, Variational::H1<4, ComplexMatrix>>,
      MatrixShapeCase<5, Variational::H1<5, ComplexMatrix>>,
      MatrixShapeCase<6, Variational::H1<6, ComplexMatrix>>>;
  TYPED_TEST_SUITE(MatrixShapeLifecycleTest, MatrixShapeCases);

  TYPED_TEST(MatrixShapeLifecycleTest, SpaceAssignmentRefreshesLogicalSupport)
  {
    using namespace Geometry;
    using namespace Variational;
    using Order = typename TypeParam::first_type;
    using FES = typename TypeParam::second_type;
    using Scalar = typename FormLanguage::Traits<FES>::ScalarType;
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
      const auto makeSpace = [&](size_t rows, size_t columns) {
        if constexpr (std::is_constructible_v<FES, Order, const LocalMesh&, size_t,
                        size_t>)
          return FES(Order{}, mesh, rows, columns);
        else
          return FES(mesh, rows, columns);
      };
      const auto cell = *mesh.getCell();
      const auto& formula = QF::PolytopeQuadratureFormula::get(4, geometry);
      const PolytopeQuadrature quadrature(cell, formula);
      const IntegrationPoint cursor(quadrature, formula.getSize() - 1);
      auto space = makeSpace(2, 3);
      TrialFunction original(space);
      original.getSolution().getData().setZero();
      const auto check = [&](auto& unknown) {
        auto& shape = unknown.getDerived();
        shape.setIntegrationPoint(cursor);
        const auto& element =
          space.getFiniteElement(cell.getDimension(), cell.getIndex());
        const auto& scalar = element.getScalarElement();
        const size_t rows = space.getRows(), columns = space.getColumns();
        const size_t components = rows * columns;
        ASSERT_EQ(&shape.getFiniteElementSpace(), &space);
        ASSERT_EQ(shape.getDOFs(cell), scalar.getCount() * components);
        for (size_t a = 0; a < scalar.getCount(); ++a)
        {
          const Scalar value =
            scalar.getBasis(a)(cursor.getPoint().getReferenceCoordinates());
          for (size_t c = 0; c < components; ++c)
          {
            const auto& actual = shape.getBasis(a * components + c);
            ASSERT_EQ(actual.rows(), rows);
            ASSERT_EQ(actual.cols(), columns);
            for (size_t entry = 0; entry < components; ++entry)
            {
              ASSERT_EQ(
                actual(entry / columns, entry % columns), entry == c ? value : Scalar(0));
            }
          }
        }
      };
      check(original);
      size_t transition = 0;
      for (const auto& [rows, columns] : {std::pair<size_t, size_t>{3, 2}, {1, 2}, {2, 3},
             {3, 3}, {1, 1}, {1, 3}, {3, 1}, {2, 2}, {2, 1}, {2, 3}})
      {
        SCOPED_TRACE(::testing::Message() << rows << "x" << columns);
        auto replacement = makeSpace(rows, columns);
        if (transition++ % 2 == 0)
          space = replacement;
        else
          space = std::move(replacement);
        ASSERT_EQ(space.getRows(), rows);
        ASSERT_EQ(space.getColumns(), columns);
        // Coefficients do not determine basis values, but the copied unknowns
        // must still carry a valid vector for the newly assigned space.
        original.getSolution().getData().resize(space.getSize());
        original.getSolution().getData().setZero();
        ASSERT_EQ(cursor.getReferenceQuadratureFormula(), &formula);
        check(original);
        check(original);
        TrialFunction copied(original);
        check(copied);
        TrialFunction moved(std::move(copied));
        check(moved);
      }
    }
  }
}
