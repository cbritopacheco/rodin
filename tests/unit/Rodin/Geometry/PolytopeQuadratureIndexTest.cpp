/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/PolytopeQuadrature.h>
#include <Rodin/QF/Centroid.h>
#include <Rodin/QF/GenericPolytopeQuadrature.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  namespace
  {
    Mesh<Context::Local> makeTriangleMesh()
    {
      return Mesh<Context::Local>::Builder()
        .initialize(2)
        .nodes(3)
        .vertex({0, 0})
        .vertex({1, 0})
        .vertex({0, 1})
        .polytope(Polytope::Type::Triangle, {0, 1, 2})
        .finalize();
    }
  }

  TEST(Geometry_PolytopeQuadratureIndex, CachesSingleFormulaPerPolytope)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    QF::Centroid qf(Polytope::Type::Triangle);
    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    size_t factoryCalls = 0;
    const auto factory = [&]()
    {
      factoryCalls++;
      return std::make_unique<PolytopeQuadrature>(polytope, qf);
    };

    const auto& q1 = index.get({2, 0}, mesh.getPolytopeCount(2), qf, factory);
    const auto& q2 = index.get({2, 0}, mesh.getPolytopeCount(2), qf, factory);

    EXPECT_EQ(&q1, &q2);
    EXPECT_EQ(factoryCalls, 1);
  }

  TEST(Geometry_PolytopeQuadratureIndex, DistinguishesQuadratureFormulaKeys)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    QF::Centroid qf1(Polytope::Type::Triangle);
    const auto& qf2 = QF::GenericPolytopeQuadrature::get(2, Polytope::Type::Triangle);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    const auto& qc = index.get(
      {2, 0}, mesh.getPolytopeCount(2), qf1,
      [&]() { return std::make_unique<PolytopeQuadrature>(polytope, qf1); });
    const auto& qg = index.get(
      {2, 0}, mesh.getPolytopeCount(2), qf2,
      [&]() { return std::make_unique<PolytopeQuadrature>(polytope, qf2); });

    EXPECT_NE(&qc, &qg);
    EXPECT_EQ(&qc.getQuadratureFormula(), &qf1);
    EXPECT_EQ(&qg.getQuadratureFormula(), &qf2);
  }

  TEST(Geometry_PolytopeQuadratureIndex, ClearDropsCachedEntries)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    QF::Centroid qf(Polytope::Type::Triangle);
    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    size_t factoryCalls = 0;
    const auto factory = [&]()
    {
      factoryCalls++;
      return std::make_unique<PolytopeQuadrature>(polytope, qf);
    };

    index.get({2, 0}, mesh.getPolytopeCount(2), qf, factory);
    index.clear();
    index.resize(2, mesh.getPolytopeCount(2));
    index.get({2, 0}, mesh.getPolytopeCount(2), qf, factory);

    EXPECT_EQ(factoryCalls, 2);
  }

  TEST(Geometry_PolytopeQuadratureIndex, ThrowsOnInvalidDimension)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    QF::Centroid qf(Polytope::Type::Triangle);
    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());

    EXPECT_THROW(
      index.get(
        {3, 0},
        1,
        qf,
        [&]() { return std::make_unique<PolytopeQuadrature>(polytope, qf); }),
      std::out_of_range);
  }
}
