/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Distributed physical-derivative state regressions.
 */
#include <gtest/gtest.h>
#include <boost/mpi.hpp>
#include "Rodin/MPI/Variational.h"
#include "../../../../convergence/MPIConvergence.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  boost::mpi::environment* environment;
  boost::mpi::communicator* world;

  /**
   * @brief Checks rank-local cache invalidation with unchanged logical entities.
   *
   * Spaces are constructed collectively by the caller. Empty ranks perform no
   * evaluation; nonempty ranks scale their owned and ghost coordinates by the
   * same binary factor. Geometry updates and basis binding are local operations,
   * so no collective is permitted inside this branch. Every component must
   * equal a fresh expression and exactly half its previous value.
   *
   * @param mesh Distributed mesh whose local geometry is modified.
   * @param space Native scalar or vector P1 space.
   */
  template <class FES>
  void checkState(Mesh<Context::MPI>& mesh, FES& space)
  {
    const auto iterator = mesh.getCell();
    if (!iterator)
      return;
    const auto cell = *iterator;
    TestFunction u(space);
    const auto& formula = QF::PolytopeQuadratureFormula::get(4, cell.getGeometry());
    const auto exercise = [&](auto& cached, auto makeFresh) {
      using Value = std::remove_cvref_t<decltype(cached.getBasis(size_t(0)))>;
      using Scalar = typename FormLanguage::Traits<FES>::ScalarType;
      const PolytopeQuadrature before(cell, formula);
      const IntegrationPoint first(before, 0);
      cached.setIntegrationPoint(first);
      std::vector<Value> original;
      for (size_t a = 0; a < cached.getDOFs(cell); ++a)
        original.push_back(cached.getBasis(a));
      const auto identity = cell.getTransformation().getCacheIdentity();
      mesh.scale(2);
      const PolytopeQuadrature after(cell, formula);
      EXPECT_NE(cell.getTransformation().getCacheIdentity(), identity);
      const IntegrationPoint second(after, 0);
      cached.setIntegrationPoint(second);
      auto fresh = makeFresh(u);
      fresh.setIntegrationPoint(second);
      for (size_t a = 0; a < original.size(); ++a)
      {
        const auto actual = cached.getBasis(a);
        const auto expected = fresh.getBasis(a);
        if constexpr (FormLanguage::IsMatrixRange<Value>::Value)
        {
          EXPECT_EQ(actual.rows(), expected.rows());
          EXPECT_EQ(actual.cols(), expected.cols());
          for (Eigen::Index i = 0; i < actual.rows(); ++i)
          {
            for (Eigen::Index j = 0; j < actual.cols(); ++j)
            {
              EXPECT_EQ(actual(i, j), expected(i, j));
              EXPECT_EQ(Scalar(2) * actual(i, j), original[a](i, j));
            }
          }
        }
        else if constexpr (FormLanguage::IsVectorRange<Value>::Value)
        {
          EXPECT_EQ(actual.size(), expected.size());
          for (size_t i = 0; i < actual.size(); ++i)
          {
            EXPECT_EQ(actual(i), expected(i));
            EXPECT_EQ(Scalar(2) * actual(i), original[a](i));
          }
        }
        else
        {
          EXPECT_EQ(actual, expected);
          EXPECT_EQ(Scalar(2) * actual, original[a]);
        }
      }
      mesh.scale(Real(0.5));
    };
    if constexpr (FormLanguage::IsVectorRange<
                    typename FormLanguage::Traits<FES>::RangeType>::Value)
    {
      auto jacobian = Jacobian(u);
      exercise(jacobian, [](auto& unknown) { return Jacobian(unknown); });
      auto divergence = Div(u);
      exercise(divergence, [](auto& unknown) { return Div(unknown); });
    }
    else
    {
      auto gradient = Grad(u);
      exercise(gradient, [](auto& unknown) { return Grad(unknown); });
    }
  }
}

/// @brief Real/complex native P1 derivatives rebind over all distributed cell families.
TEST(MPI_DerivativeState, GeometryRebindAllGeometries)
{
  const Context::MPI context(*environment, *world);
  using G = Polytope::Type;
  for (const auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
         G::Pyramid, G::Hexahedron, G::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = Tests::Convergence::DistributedUniformGrid(context, geometry).makeMesh(2);
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    const auto fields = [&]<class Scalar>() {
      P1<Scalar, decltype(mesh)> scalar(mesh);
      P1<Math::SpatialVector<Scalar>, decltype(mesh)> vector(mesh, dimension);
      checkState(mesh, scalar);
      checkState(mesh, vector);
    };
    fields.template operator()<Real>();
    fields.template operator()<Complex>();
  }
}

/**
 * @brief Certifies H1 scalar/vector geometry rebinding through every supported order.
 *
 * Every rank constructs the real/complex spaces collectively before any
 * empty-holder branch. Binding, reference association and geometry updates
 * remain rank-local. Orders one through six and all seven cell families use
 * the same exact dilation/fresh-expression checks as native P1.
 */
TEST(MPI_DerivativeState, H1GeometryRebindAllOrdersAllGeometries)
{
  const Context::MPI context(*environment, *world);
  using G = Polytope::Type;
  for (const auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
         G::Pyramid, G::Hexahedron, G::Wedge})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    auto mesh = Tests::Convergence::DistributedUniformGrid(context, geometry).makeMesh(2);
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    const auto fields = [&]<class Scalar>() {
      const auto degree = [&]<size_t K>(std::integral_constant<size_t, K> order) {
        SCOPED_TRACE(K);
        H1<K, Scalar, decltype(mesh)> scalar(order, mesh);
        H1<K, Math::SpatialVector<Scalar>, decltype(mesh)> vector(order, mesh, dimension);
        checkState(mesh, scalar);
        checkState(mesh, vector);
      };
      degree(std::integral_constant<size_t, 1>{});
      degree(std::integral_constant<size_t, 2>{});
      degree(std::integral_constant<size_t, 3>{});
      degree(std::integral_constant<size_t, 4>{});
      degree(std::integral_constant<size_t, 5>{});
      degree(std::integral_constant<size_t, 6>{});
    };
    fields.template operator()<Real>();
    fields.template operator()<Complex>();
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment mpiEnvironment(argc, argv);
  boost::mpi::communicator mpiWorld;
  environment = &mpiEnvironment;
  world = &mpiWorld;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
