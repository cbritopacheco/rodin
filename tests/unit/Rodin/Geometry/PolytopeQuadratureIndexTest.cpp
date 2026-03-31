/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <vector>

#include <Rodin/Alert/Exception.h>
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
    constexpr size_t kMaxQuadraturesPerPolytope = 8;

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
      Alert::Exception);
  }

  TEST(Geometry_PolytopeQuadratureIndex, ThrowsWhenPerPolytopeCapacityExceeded)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    std::vector<std::unique_ptr<QF::Centroid>> quadratureFormulas;
    // Reserve one extra slot for the expected-to-fail insertion attempt.
    quadratureFormulas.reserve(kMaxQuadraturesPerPolytope + 1);

    // Cache keying uses quadrature formula object identity (address), so this test
    // deliberately creates distinct formula objects to fill each slot.
    for (size_t i = 0; i < kMaxQuadraturesPerPolytope; ++i)
    {
      quadratureFormulas.emplace_back(std::make_unique<QF::Centroid>(Polytope::Type::Triangle));
      EXPECT_NO_THROW(index.get(
        {2, 0},
        mesh.getPolytopeCount(2),
        *quadratureFormulas.back(),
        [&]()
        {
          return std::make_unique<PolytopeQuadrature>(polytope, *quadratureFormulas.back());
        }));
    }

    quadratureFormulas.emplace_back(std::make_unique<QF::Centroid>(Polytope::Type::Triangle));
    EXPECT_THROW(
      index.get(
        {2, 0},
        mesh.getPolytopeCount(2),
        *quadratureFormulas.back(),
        [&]()
        {
          return std::make_unique<PolytopeQuadrature>(polytope, *quadratureFormulas.back());
        }),
      Alert::Exception);
  }
}
