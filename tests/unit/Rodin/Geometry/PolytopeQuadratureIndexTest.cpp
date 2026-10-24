/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <array>
#include <atomic>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include <Rodin/Alert/Exception.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/PolytopeQuadrature.h>
#include <Rodin/QF/Centroid.h>
#include <Rodin/QF/GaussLegendre.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include "../../../convergence/Convergence.h"

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

  /// @brief Supplies equal centroid samples under distinct factory keys for index tests.
  class IndexedCentroid final : public QF::QuadratureFormulaBase
  {
    public:
      /**
       * @param geometry Reference-cell geometry.
       * @param tag Distinguishes numerically identical inputs to the tested index.
       */
      IndexedCentroid(Polytope::Type geometry, [[maybe_unused]] size_t tag)
        : m_formula(geometry)
      {}

      /// @returns One centroid sample.
      size_t getSize() const override
      {
        return m_formula.getSize();
      }

      /**
       * @param index Sample index.
       * @returns Reference measure.
       */
      Real getWeight(size_t index) const override
      {
        return m_formula.getWeight(index);
      }

      /**
       * @param index Sample index.
       * @returns Reference centroid.
       */
      const Math::SpatialPoint& getPoint(size_t index) const override
      {
        return m_formula.getPoint(index);
      }

    private:
      QF::Centroid m_formula;
  };

  /// @brief Rejected point indices preserve the mapped point and reference formula.
  TEST(Geometry_PolytopeQuadrature, CheckedPointBoundsUseMemberContext)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);
    const auto& formula =
      QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);
    const PolytopeQuadrature quadrature(polytope, formula);
    ASSERT_EQ(quadrature.getSize(), 1u);
    const auto* stored = &quadrature.getPoint(0);
    const auto coordinates = stored->getPhysicalCoordinates();
    for (size_t index : {size_t{1}, std::numeric_limits<size_t>::max()})
    {
      try
      {
        quadrature.getPoint(index);
        ADD_FAILURE() << "Out-of-range quadrature points must be rejected.";
      }
      catch (const Alert::Exception& error)
      {
        EXPECT_NE(std::string(error.what()).find("getPoint"), std::string::npos);
        EXPECT_NE(
          std::string(error.what()).find("In member function"), std::string::npos);
      }
      catch (...)
      {
        ADD_FAILURE() << "Checked point access must raise an Alert exception.";
      }
    }
    EXPECT_EQ(quadrature.getSize(), 1u);
    EXPECT_EQ(&quadrature.getQuadratureFormula(), &formula);
    EXPECT_EQ(&quadrature.getPoint(0), stored);
    EXPECT_EQ((stored->getPhysicalCoordinates() - coordinates).norm(), 0);
  }

  class QuadratureVersionTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  /**
   * @brief Formula assignment and storage reuse create distinct stable entries.
   *
   * Retained mapped samples keep their original coordinates until mesh flush.
   * New requests use the current formula version, even at the same address.
   * More than sixteen versions exercise the full index beyond its fast lookup.
   */
  TEST_P(QuadratureVersionTest, FormulaVersionsRetainBorrowedSamples)
  {
    auto mesh = GetParam() == Polytope::Type::Point
      ? LocalMesh::Builder().initialize(3).nodes(1).vertex({0, 0, 0}).finalize()
      : Convergence::UniformGrid(GetParam()).makeMesh(3);
    const auto cell = *mesh.getCell();
    std::vector<const PolytopeQuadrature*> retained;
    std::vector<const QF::QuadratureFormula*> formulas;
    std::vector<Math::SpatialPoint> coordinates;
    for (size_t version = 0; version < 20; ++version)
    {
      const auto& formula =
        QF::QuadratureFormula::get<QF::GaussLegendre>(GetParam(), version + 2);
      const auto& quadrature = cell.getQuadrature(formula);
      EXPECT_EQ(&quadrature, &cell.getQuadrature(formula));
      EXPECT_EQ(&quadrature.getQuadratureFormula(), &formula);
      ASSERT_EQ(quadrature.getSize(), formula.getSize());
      for (size_t q = 0; q < quadrature.getSize(); ++q)
      {
        EXPECT_EQ(
          (quadrature.getPoint(q).getReferenceCoordinates() - formula.getPoint(q)).norm(),
          0);
      }
      for (const auto* previous : retained)
        EXPECT_NE(&quadrature, previous);
      retained.push_back(&quadrature);
      formulas.push_back(&formula);
      coordinates.push_back(formula.getPoint(0));
      for (size_t old = 0; old < retained.size(); ++old)
      {
        EXPECT_EQ(&retained[old]->getQuadratureFormula(), formulas[old]);
        EXPECT_EQ(
          (retained[old]->getPoint(0).getReferenceCoordinates() - coordinates[old])
            .norm(),
          0);
      }
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, QuadratureVersionTest,
    ::testing::Values(Polytope::Type::Point, Polytope::Type::Segment,
      Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge));

  /// @brief Verifies caches single formula per polytope for geometry polytope quadrature index by checking exact expected values.
  TEST(Geometry_PolytopeQuadratureIndex, CachesSingleFormulaPerPolytope)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    size_t factoryCalls = 0;
    const auto factory = [&]() {
      ++factoryCalls;
      return std::make_unique<PolytopeQuadrature>(polytope, qf);
    };

    const auto& q1 = index.get({2, 0}, mesh.getPolytopeCount(2), qf, factory);
    const auto& q2 = index.get({2, 0}, mesh.getPolytopeCount(2), qf, factory);

    EXPECT_EQ(&q1, &q2);
    EXPECT_EQ(factoryCalls, 1u);
    EXPECT_EQ(&q1.getQuadratureFormula(), &qf);
  }

  /// @brief Verifies distinguishes quadrature formula keys for geometry polytope quadrature index by checking exact expected values.
  TEST(Geometry_PolytopeQuadratureIndex, DistinguishesQuadratureFormulaKeys)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    const auto& qf1 = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);
    const auto& qf2 = QF::PolytopeQuadratureFormula::get(2, Polytope::Type::Triangle);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    const auto& qc = index.get({2, 0}, mesh.getPolytopeCount(2), qf1,
      [&]() { return std::make_unique<PolytopeQuadrature>(polytope, qf1); });

    const auto& qg = index.get({2, 0}, mesh.getPolytopeCount(2), qf2,
      [&]() { return std::make_unique<PolytopeQuadrature>(polytope, qf2); });

    EXPECT_NE(&qc, &qg);
    EXPECT_EQ(&qc.getQuadratureFormula(), &qf1);
    EXPECT_EQ(&qg.getQuadratureFormula(), &qf2);
  }

  /// @brief Verifies clear drops cached entries for geometry polytope quadrature index by checking exact expected values.
  TEST(Geometry_PolytopeQuadratureIndex, ClearDropsCachedEntries)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    size_t factoryCalls = 0;
    const auto factory = [&]() {
      ++factoryCalls;
      return std::make_unique<PolytopeQuadrature>(polytope, qf);
    };

    index.get({2, 0}, mesh.getPolytopeCount(2), qf, factory);
    index.clear();
    index.resize(2, mesh.getPolytopeCount(2));
    const auto& q2 = index.get({2, 0}, mesh.getPolytopeCount(2), qf, factory);
    const auto& q3 = index.get({2, 0}, mesh.getPolytopeCount(2), qf, factory);

    EXPECT_EQ(factoryCalls, 2u);
    EXPECT_EQ(&q2, &q3);
    EXPECT_EQ(&q2.getQuadratureFormula(), &qf);
  }

  /// @brief Verifies mesh move constructor drops stale entries for geometry polytope quadrature index by checking exact expected values, move semantics.
  TEST(Geometry_PolytopeQuadratureIndex, MeshMoveConstructorDropsStaleEntries)
  {
    auto mesh = makeTriangleMesh();
    const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);

    mesh.getQuadrature(2, 0, qf);

    Mesh<Context::Local> moved(std::move(mesh));
    const auto& quadrature = moved.getQuadrature(2, 0, qf);
    const auto& point = quadrature.getPoint(0);

    EXPECT_EQ(&point.getPolytope().getMesh(), &static_cast<const MeshBase&>(moved));
  }

  /// @brief Verifies mesh move assignment drops stale entries for geometry polytope quadrature index by checking exact expected values, move semantics.
  TEST(Geometry_PolytopeQuadratureIndex, MeshMoveAssignmentDropsStaleEntries)
  {
    auto mesh = makeTriangleMesh();
    const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);

    mesh.getQuadrature(2, 0, qf);

    Mesh<Context::Local> moved;
    moved = std::move(mesh);
    const auto& quadrature = moved.getQuadrature(2, 0, qf);
    const auto& point = quadrature.getPoint(0);

    EXPECT_EQ(&point.getPolytope().getMesh(), &static_cast<const MeshBase&>(moved));
  }

  /// @brief Verifies throws on invalid dimension for geometry polytope quadrature index by checking exception behavior.
  TEST(Geometry_PolytopeQuadratureIndex, ThrowsOnInvalidDimension)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());

    EXPECT_THROW(index.get({3, 0}, 1, qf,
                   [&]() { return std::make_unique<PolytopeQuadrature>(polytope, qf); }),
      Alert::Exception);
  }

  /// @brief Verifies throws on out of range polytope index for geometry polytope quadrature index by checking exception behavior.
  TEST(Geometry_PolytopeQuadratureIndex, ThrowsOnOutOfRangePolytopeIndex)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    EXPECT_THROW(index.get({2, mesh.getPolytopeCount(2)}, mesh.getPolytopeCount(2), qf,
                   [&]() { return std::make_unique<PolytopeQuadrature>(polytope, qf); }),
      Alert::Exception);
  }

  /// @brief Verifies exceeding capacity does not throw for geometry polytope quadrature index by checking no-throw behavior.
  TEST(Geometry_PolytopeQuadratureIndex, ExceedingCapacityDoesNotThrow)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    std::vector<const QF::QuadratureFormula*> quadratureFormulas;
    quadratureFormulas.reserve(kMaxQuadraturesPerPolytope + 1);

    for (size_t i = 0; i < kMaxQuadraturesPerPolytope + 1; ++i)
    {
      quadratureFormulas.emplace_back(&QF::QuadratureFormula::get<IndexedCentroid>(
        Polytope::Type::Triangle, quadratureFormulas.size()));

      EXPECT_NO_THROW(
        index.get({2, 0}, mesh.getPolytopeCount(2), *quadratureFormulas.back(), [&]() {
          return std::make_unique<PolytopeQuadrature>(
            polytope, *quadratureFormulas.back());
        }));
    }
  }

  /// @brief Verifies formula entries remain stable beyond the former bounded-cache capacity.
  TEST(Geometry_PolytopeQuadratureIndex, RetainsFormulaBeyondFormerCapacity)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    std::vector<const QF::QuadratureFormula*> quadratureFormulas;
    quadratureFormulas.reserve(kMaxQuadraturesPerPolytope + 1);

    size_t factoryCalls = 0;

    // Fill the slot with distinct formula identities.
    for (size_t i = 0; i < kMaxQuadraturesPerPolytope; ++i)
    {
      quadratureFormulas.emplace_back(&QF::QuadratureFormula::get<IndexedCentroid>(
        Polytope::Type::Triangle, quadratureFormulas.size()));

      index.get({2, 0}, mesh.getPolytopeCount(2), *quadratureFormulas.back(), [&]() {
        ++factoryCalls;
        return std::make_unique<PolytopeQuadrature>(polytope, *quadratureFormulas.back());
      });
    }

    EXPECT_EQ(factoryCalls, kMaxQuadraturesPerPolytope);

    // Insert one more distinct formula.
    quadratureFormulas.emplace_back(&QF::QuadratureFormula::get<IndexedCentroid>(
      Polytope::Type::Triangle, quadratureFormulas.size()));

    index.get({2, 0}, mesh.getPolytopeCount(2), *quadratureFormulas.back(), [&]() {
      ++factoryCalls;
      return std::make_unique<PolytopeQuadrature>(polytope, *quadratureFormulas.back());
    });

    EXPECT_EQ(factoryCalls, kMaxQuadraturesPerPolytope + 1);

    // Stable formula storage retains the first quadrature until clear().
    index.get({2, 0}, mesh.getPolytopeCount(2), *quadratureFormulas.front(), [&]() {
      ++factoryCalls;
      return std::make_unique<PolytopeQuadrature>(polytope, *quadratureFormulas.front());
    });

    EXPECT_EQ(factoryCalls, kMaxQuadraturesPerPolytope + 1);
  }

  /// @brief Verifies immediate repeated access after registering several formulas.
  TEST(Geometry_PolytopeQuadratureIndex, PreservesImmediateRepeatAcrossFormulas)
  {
    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());
    index.resize(2, mesh.getPolytopeCount(2));

    std::vector<const QF::QuadratureFormula*> quadratureFormulas;
    quadratureFormulas.reserve(kMaxQuadraturesPerPolytope + 1);

    size_t factoryCalls = 0;

    for (size_t i = 0; i < kMaxQuadraturesPerPolytope + 1; ++i)
    {
      quadratureFormulas.emplace_back(&QF::QuadratureFormula::get<IndexedCentroid>(
        Polytope::Type::Triangle, quadratureFormulas.size()));

      index.get({2, 0}, mesh.getPolytopeCount(2), *quadratureFormulas.back(), [&]() {
        ++factoryCalls;
        return std::make_unique<PolytopeQuadrature>(polytope, *quadratureFormulas.back());
      });
    }

    EXPECT_EQ(factoryCalls, kMaxQuadraturesPerPolytope + 1);

    // Immediate repeat of the most recently inserted formula should not rebuild.
    const auto& q1 =
      index.get({2, 0}, mesh.getPolytopeCount(2), *quadratureFormulas.back(), [&]() {
        ++factoryCalls;
        return std::make_unique<PolytopeQuadrature>(polytope, *quadratureFormulas.back());
      });

    const auto& q2 =
      index.get({2, 0}, mesh.getPolytopeCount(2), *quadratureFormulas.back(), [&]() {
        ++factoryCalls;
        return std::make_unique<PolytopeQuadrature>(polytope, *quadratureFormulas.back());
      });

    EXPECT_EQ(&q1, &q2);
    EXPECT_EQ(factoryCalls, kMaxQuadraturesPerPolytope + 1);
  }

  /// @brief Verifies that concurrent formula lookup constructs each cached entry once.
  TEST(Geometry_PolytopeQuadratureIndex, ConstructsDistinctFormulasConcurrently)
  {
    constexpr size_t formulaCount = 8;
    constexpr size_t repetitions = 64;

    auto mesh = makeTriangleMesh();
    const auto polytope = *mesh.getPolytope(2, 0);

    PolytopeQuadratureIndex index;
    index.initialize(mesh.getDimension());

    std::array<const QF::QuadratureFormula*, formulaCount> formulas;
    std::array<std::atomic<size_t>, formulaCount> factoryCalls{};
    for (size_t i = 0; i < formulas.size(); ++i)
    {
      formulas[i] =
        &QF::QuadratureFormula::get<IndexedCentroid>(Polytope::Type::Triangle, i);
    }

    std::vector<std::thread> threads;
    threads.reserve(formulaCount);
    for (size_t i = 0; i < formulaCount; ++i)
    {
      threads.emplace_back([&, i]() {
        for (size_t repetition = 0; repetition < repetitions; ++repetition)
        {
          const auto& quadrature =
            index.get({2, 0}, mesh.getPolytopeCount(2), *formulas[i], [&]() {
              ++factoryCalls[i];
              return std::make_unique<PolytopeQuadrature>(polytope, *formulas[i]);
            });
          EXPECT_EQ(&quadrature.getQuadratureFormula(), formulas[i]);
        }
      });
    }

    for (auto& thread : threads)
      thread.join();

    for (const auto& calls : factoryCalls)
      EXPECT_EQ(calls.load(), 1u);
  }
}
