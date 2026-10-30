/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Essential elimination on low/high-order and mixed assembled systems.
 *
 * The affine solution satisfies @f$-\Delta u+u=u@f$ with its exact trace.
 * The mixed-field check uses a one-way mass coupling, deliberately omitting
 * its transpose block. Both tests independently solve the eliminated system.
 */
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Geometry and polynomial degree for the manufactured system.
  class SparseElimination
    : public ::testing::TestWithParam<std::pair<Polytope::Type, size_t>>
  {
    public:
      /**
       * @brief Checks the affine reaction-diffusion problem and a mixed system.
       * @tparam Degree Dedicated P1 or high-order H1 degree.
       * @param mixedGeometry Use a conforming triangle/quadrilateral mesh in two dimensions.
       */
      template <size_t Degree>
      void check(bool mixedGeometry = false)
      {
        const auto geometry = GetParam().first;
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        Array<size_t> shape(dimension);
        shape.setConstant(3);
        auto mesh = [&]() {
          if (!mixedGeometry)
            return LocalMesh::UniformGrid(geometry, shape);
          auto builder = LocalMesh::Builder();
          builder.initialize(2).nodes(6);
          builder.vertex({0, 0})
            .vertex({1, 0})
            .vertex({2, 0})
            .vertex({0, 1})
            .vertex({1, 1})
            .vertex({2, 1});
          builder.polytope(Polytope::Type::Triangle, {0, 1, 3})
            .polytope(Polytope::Type::Triangle, {1, 4, 3})
            .polytope(Polytope::Type::Quadrilateral, {1, 2, 5, 4});
          return builder.finalize();
        }();
        mesh.scale(Real(0.5));
        for (size_t high = 1; high <= dimension; ++high)
        {
          for (size_t low = 0; low < high; ++low)
            mesh.getConnectivity().compute(high, low);
        }
        mesh.getConnectivity().compute(dimension - 1, dimension);
        auto space = [&]() {
          if constexpr (Degree == 1)
            return P1<Real>(mesh);
          else
            return H1<Degree, Real>(std::integral_constant<size_t, Degree>{}, mesh);
        }();
        TrialFunction u(space);
        TestFunction v(space);
        RealFunction exact([](const Point& point) {
          Real value = 1;
          for (size_t axis = 0; axis < point.getCoordinates().size(); ++axis)
            value += Real(axis + 1) * point.getCoordinates()[axis];
          return value;
        });
        GridFunction coefficients(space);
        coefficients = exact;
        auto stiffness = Integral(Grad(u), Grad(v));
        auto mass = Integral(u, v);
        auto force = Integral(exact, v);
        // Rational pyramid bases need the established manufactured-test rule.
        constexpr size_t PyramidQuadratureOrder = 12;
        const size_t quadratureOrder =
          geometry == Polytope::Type::Pyramid ? PyramidQuadratureOrder : 2 * Degree;
        stiffness.setOrder(quadratureOrder);
        mass.setOrder(quadratureOrder);
        force.setOrder(quadratureOrder);
        BilinearForm form(u, v);
        form = stiffness + mass;
        form.assemble();
        LinearForm load(v);
        load = force;
        load.assemble();
        auto boundary = DirichletBC(u, exact);
        boundary.assemble();
        const auto& fixed =
          std::get<typename decltype(boundary)::ValueDOFs>(boundary.getDOFs());
        using Matrix = Math::SparseMatrix<Real>;
        using Vector = Math::Vector<Real>;
        using System = Math::LinearSystem<Matrix, Vector>;
        System system;
        system.getOperator() = form.getOperator();
        system.getVector() = load.getVector();
        system.eliminate(fixed);
        const Math::Matrix<Real> dense = system.getOperator();
        const Vector solution = dense.fullPivLu().solve(system.getVector());
        EXPECT_LT((solution - coefficients.getData()).norm(), 2e-9);
        EXPECT_LT(
          (system.getOperator() * coefficients.getData() - system.getVector()).norm(),
          2e-10);

        H1<1, Real> low(std::integral_constant<size_t, 1>{}, mesh);
        TrialFunction p(low);
        TestFunction q(low);
        BilinearForm diagonal(p, q);
        diagonal = Integral(Grad(p), Grad(q)) + Integral(p, q);
        diagonal.assemble();
        auto mixedMass = Integral(p, v);
        mixedMass.setOrder(quadratureOrder);
        BilinearForm coupling(p, v);
        coupling = mixedMass;
        coupling.assemble();
        GridFunction lowExact(low);
        lowExact = exact;
        const Index highCount = form.getOperator().rows();
        const Index lowCount = diagonal.getOperator().rows();
        std::vector<Math::SparseTriplet<Real>> entries;
        for (Index column = 0; column < highCount; ++column)
        {
          for (Matrix::InnerIterator it(form.getOperator(), column); it; ++it)
            entries.emplace_back(it.row(), column, it.value());
        }
        for (Index column = 0; column < lowCount; ++column)
        {
          for (Matrix::InnerIterator it(diagonal.getOperator(), column); it; ++it)
            entries.emplace_back(highCount + it.row(), highCount + column, it.value());
          for (Matrix::InnerIterator it(coupling.getOperator(), column); it; ++it)
            entries.emplace_back(it.row(), highCount + column, -Real(0.3) * it.value());
        }
        Matrix original(highCount + lowCount, highCount + lowCount);
        original.setFromTriplets(entries.begin(), entries.end());
        Vector mixedExact(highCount + lowCount);
        mixedExact << coefficients.getData(), lowExact.getData();
        auto lowBoundary = DirichletBC(p, exact);
        lowBoundary.assemble();
        IndexMap<Real> mixedFixed = fixed;
        for (const auto& [index, value] :
          std::get<typename decltype(lowBoundary)::ValueDOFs>(lowBoundary.getDOFs()))
          mixedFixed.emplace(highCount + index, value);
        Vector prescribed = Vector::Zero(mixedExact.size());
        Eigen::DiagonalMatrix<Real, Eigen::Dynamic> free(mixedExact.size());
        free.diagonal().setOnes();
        for (const auto& [index, value] : mixedFixed)
        {
          prescribed[index] = value;
          free.diagonal()[index] = 0;
        }
        Math::Matrix<Real> expected = free * original * free;
        Vector expectedLoad =
          free * (original * mixedExact - original * prescribed) + prescribed;
        for (const auto& [index, value] : mixedFixed)
          expected(index, index) = 1;
        System mixed;
        mixed.getOperator() = original;
        mixed.getVector() = original * mixedExact;
        mixed.eliminate(mixedFixed);
        EXPECT_LT((Math::Matrix<Real>(mixed.getOperator()) - expected).norm(), 1e-11);
        EXPECT_LT((mixed.getVector() - expectedLoad).norm(), 1e-10);
        const Vector mixedSolution =
          Math::Matrix<Real>(mixed.getOperator()).fullPivLu().solve(mixed.getVector());
        EXPECT_LT((mixedSolution - mixedExact).norm(), 3e-9);
        EXPECT_LT((mixed.getOperator() * mixedExact - mixed.getVector()).norm(), 2e-10);
      }
  };
}

/// @brief Manufactured affine PDE and asymmetric P1/H1 mixed-field elimination.
TEST_P(SparseElimination, AffineReactionDiffusionAndOneWayMixedCoupling)
{
  if (GetParam().second == 1)
    check<1>();
  else
    check<3>();
  if (GetParam().first == Polytope::Type::Quadrilateral)
  {
    if (GetParam().second == 1)
      check<1>(true);
    else
      check<3>(true);
  }
}

INSTANTIATE_TEST_SUITE_P(AllGeometries, SparseElimination,
  ::testing::Values(std::pair{Polytope::Type::Segment, size_t(1)},
    std::pair{Polytope::Type::Triangle, size_t(1)},
    std::pair{Polytope::Type::Quadrilateral, size_t(1)},
    std::pair{Polytope::Type::Tetrahedron, size_t(1)},
    std::pair{Polytope::Type::Hexahedron, size_t(1)},
    std::pair{Polytope::Type::Pyramid, size_t(1)},
    std::pair{Polytope::Type::Wedge, size_t(1)},
    std::pair{Polytope::Type::Segment, size_t(3)},
    std::pair{Polytope::Type::Triangle, size_t(3)},
    std::pair{Polytope::Type::Quadrilateral, size_t(3)},
    std::pair{Polytope::Type::Tetrahedron, size_t(3)},
    std::pair{Polytope::Type::Hexahedron, size_t(3)},
    std::pair{Polytope::Type::Pyramid, size_t(3)},
    std::pair{Polytope::Type::Wedge, size_t(3)}));
