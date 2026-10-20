/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <array>
#include <atomic>
#include <barrier>
#include <thread>
#include <type_traits>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/QF/Centroid.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/QF/GaussLobatto.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/QF/TensorProduct.h"
#include "Rodin/Geometry/PolytopeQuadrature.h"
#include "Rodin/Variational/P1/P1Element.h"
#include "Rodin/Variational/IntegrationPoint.h"

using namespace Rodin;

namespace
{
  using G = Geometry::Polytope::Type;

  template <class Formula>
  concept CanCacheMappedSamples = requires(const Geometry::Polytope& cell,
    const Formula& formula) { cell.getQuadrature(formula); };

  template <class Formula>
  concept CanCacheBasis = requires(const Variational::RealP1Element& element,
    const Formula& formula) { element.getTabulation(formula); };

  static_assert(!std::is_constructible_v<Geometry::PolytopeQuadrature,
                const Geometry::Polytope&, const QF::GaussLegendre&>);
  static_assert(!std::is_constructible_v<Variational::IntegrationPoint,
                const Geometry::Point&, const QF::GaussLegendre*, size_t>);
  static_assert(!CanCacheMappedSamples<QF::GaussLegendre>);
  static_assert(!CanCacheBasis<QF::GaussLegendre>);
  static_assert(CanCacheMappedSamples<QF::QuadratureFormula>);
  static_assert(CanCacheBasis<QF::QuadratureFormula>);
  static_assert(!std::is_copy_constructible_v<QF::QuadratureFormula>);
  static_assert(!std::is_move_constructible_v<QF::QuadratureFormula>);
  static_assert(!std::is_copy_assignable_v<QF::QuadratureFormula>);
  static_assert(!std::is_move_assignable_v<QF::QuadratureFormula>);

  /// @brief Generator with controlled publication failures and sample observations.
  class FactoryFormula final : public QF::QuadratureFormulaBase
  {
    public:
      inline static std::atomic<size_t> constructions{0};
      inline static std::atomic<size_t> reads{0};
      inline static std::atomic<bool> fail{false};
      inline static std::atomic<bool> failSample{false};
      inline static std::atomic<size_t> nextTag{1000};
      inline static std::barrier<>* gate = nullptr;

      /**
       * @param tag Integer selecting this generator's reference coordinate.
       */
      explicit FactoryFormula(size_t tag)
        : m_point{Real(tag) / 1024}
      {
        ++constructions;
        if (gate)
          gate->arrive_and_wait();
        if (fail.load())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Interrupted test formula generation." << Alert::Raise;
        }
      }

      /// @returns One sample, recording the generator read.
      size_t getSize() const override
      {
        ++reads;
        return 1;
      }

      /**
       * @param index Unused index of the only sample.
       * @returns Unit weight.
       */
      Real getWeight([[maybe_unused]] size_t index) const override
      {
        ++reads;
        return 1;
      }

      /**
       * @param index Unused index of the only sample.
       * @returns Owned reference coordinate.
       */
      const Math::SpatialPoint& getPoint([[maybe_unused]] size_t index) const override
      {
        ++reads;
        if (failSample.load())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Interrupted test sample materialization." << Alert::Raise;
        }
        return m_point;
      }

    private:
      Math::SpatialPoint m_point;
  };
}

/// @brief Equivalent argument qualifications and family factories share publication.
TEST(QuadratureFormula, ReusesNormalizedArguments)
{
  auto geometry = G::Segment;
  size_t count = 4;
  const auto& first = QF::QuadratureFormula::get<QF::GaussLegendre>(geometry, count);
  const auto& second =
    QF::QuadratureFormula::get<QF::GaussLegendre>(G::Segment, size_t{4});
  EXPECT_EQ(&first, &second);
  EXPECT_EQ(&QF::GaussLobatto::get(geometry, count),
    &QF::QuadratureFormula::get<QF::GaussLobatto>(geometry, count));
  EXPECT_EQ(&QF::PolytopeQuadratureFormula::get(count, geometry),
    &QF::QuadratureFormula::get<QF::PolytopeQuadratureFormula>(count, geometry));
}

/// @brief Distinct families and parameter values never alias published objects.
TEST(QuadratureFormula, SeparatesFamiliesAndKeys)
{
  const auto& first = QF::QuadratureFormula::get<QF::GaussLegendre>(G::Point, size_t{2});
  const auto& other = QF::QuadratureFormula::get<QF::GaussLobatto>(G::Point, size_t{2});
  const auto& count = QF::QuadratureFormula::get<QF::GaussLegendre>(G::Point, size_t{3});
  EXPECT_NE(&first, &other);
  EXPECT_NE(&first, &count);
  EXPECT_EQ(first.getWeight(0), other.getWeight(0));
}

/// @brief Tuple pointer ordering distinguishes factor order and repeated factors.
TEST(QuadratureFormula, TensorProductKeysPreserveFactorOrder)
{
  const auto& first =
    QF::QuadratureFormula::get<QF::GaussLegendre>(G::Segment, size_t{2});
  const auto& second = QF::GaussLobatto::get(G::Segment, 3);
  const auto& product =
    QF::QuadratureFormula::get<QF::TensorProduct>(G::Quadrilateral, first, second);
  const auto& reversed =
    QF::QuadratureFormula::get<QF::TensorProduct>(G::Quadrilateral, second, first);
  EXPECT_EQ(&product,
    &QF::QuadratureFormula::get<QF::TensorProduct>(G::Quadrilateral, first, second));
  EXPECT_NE(&product, &reversed);
  const QF::TensorProduct expected(G::Quadrilateral, first, second);
  ASSERT_EQ(product.getSize(), expected.getSize());
  for (size_t q = 0; q < product.getSize(); ++q)
  {
    EXPECT_EQ(product.getWeight(q), expected.getWeight(q));
    EXPECT_EQ((product.getPoint(q) - expected.getPoint(q)).norm(), 0);
  }
  const auto& cube =
    QF::QuadratureFormula::get<QF::TensorProduct>(G::Hexahedron, first, first, first);
  EXPECT_EQ(cube.getSize(), first.getSize() * first.getSize() * first.getSize());
}

/// @brief Reallocation and thread-local eviction preserve formula and point addresses.
TEST(QuadratureFormula, PoolGrowthPreservesBorrowedAddresses)
{
  const auto& first = QF::QuadratureFormula::get<FactoryFormula>(size_t{1});
  const auto* point = &first.getPoint(0);
  for (size_t tag = 2; tag < 64; ++tag)
    QF::QuadratureFormula::get<FactoryFormula>(tag);
  EXPECT_EQ(&first, &QF::QuadratureFormula::get<FactoryFormula>(size_t{1}));
  EXPECT_EQ(point, &first.getPoint(0));
  EXPECT_EQ((*point)(0), Real(1) / 1024);
}

/// @brief Failed construction publishes no formula and a later request can retry.
TEST(QuadratureFormula, RetriesFailedConstruction)
{
  const size_t tag = FactoryFormula::nextTag.fetch_add(1);
  const size_t before = FactoryFormula::constructions.load();
  FactoryFormula::fail = true;
  EXPECT_THROW(QF::QuadratureFormula::get<FactoryFormula>(tag), Alert::Exception);
  FactoryFormula::fail = false;
  const auto& formula = QF::QuadratureFormula::get<FactoryFormula>(tag);
  EXPECT_EQ(FactoryFormula::constructions.load(), before + 2);
  EXPECT_EQ(formula.getPoint(0)(0), Real(tag) / 1024);
}

/// @brief Warm factory and table accesses never query the destroyed generator.
TEST(QuadratureFormula, PublishedSamplesAreIndependentOfGenerator)
{
  const auto& formula = QF::QuadratureFormula::get<FactoryFormula>(size_t{101});
  const size_t reads = FactoryFormula::reads.load();
  const size_t constructions = FactoryFormula::constructions.load();
  Variational::RealP1Element element(G::Segment);
  const auto& table = element.getTabulation(formula);
  EXPECT_EQ(&table, &element.getTabulation(formula));
  EXPECT_EQ(&formula, &QF::QuadratureFormula::get<FactoryFormula>(size_t{101}));
  EXPECT_EQ(FactoryFormula::reads.load(), reads);
  EXPECT_EQ(FactoryFormula::constructions.load(), constructions);
  EXPECT_EQ(table.getBasis(0, 1), formula.getPoint(0)(0));
}

/// @brief Simultaneous cold requests publish one shared object surviving producer exit.
TEST(QuadratureFormula, ConcurrentPublicationSurvivesThreadExit)
{
  constexpr size_t ThreadCount = 8;
  const size_t tag = FactoryFormula::nextTag.fetch_add(1);
  std::barrier start(ThreadCount);
  std::barrier publication(ThreadCount);
  const size_t before = FactoryFormula::constructions.load();
  FactoryFormula::gate = &publication;
  std::array<const QF::QuadratureFormula*, ThreadCount> formulas{};
  std::array<std::thread, ThreadCount> threads;
  for (size_t index = 0; index < ThreadCount; ++index)
  {
    threads[index] = std::thread([&, index] {
      start.arrive_and_wait();
      formulas[index] = &QF::QuadratureFormula::get<FactoryFormula>(tag);
    });
  }
  for (auto& thread : threads)
    thread.join();
  FactoryFormula::gate = nullptr;
  EXPECT_EQ(FactoryFormula::constructions.load(), before + ThreadCount);
  ASSERT_NE(formulas[0], nullptr);
  for (const auto* formula : formulas)
  {
    EXPECT_EQ(formula, formulas[0]);
    EXPECT_EQ(formula->getPoint(0)(0), Real(tag) / 1024);
  }
  EXPECT_EQ(formulas[0], &QF::QuadratureFormula::get<FactoryFormula>(tag));
}

/// @brief Interrupted sample materialization leaves no partially published object.
TEST(QuadratureFormula, RetriesFailedMaterialization)
{
  const size_t tag = FactoryFormula::nextTag.fetch_add(1);
  const size_t before = FactoryFormula::constructions.load();
  FactoryFormula::failSample = true;
  EXPECT_THROW(QF::QuadratureFormula::get<FactoryFormula>(tag), Alert::Exception);
  FactoryFormula::failSample = false;
  const auto& formula = QF::QuadratureFormula::get<FactoryFormula>(tag);
  EXPECT_EQ(FactoryFormula::constructions.load(), before + 2);
  EXPECT_EQ(formula.getPoint(0)(0), Real(tag) / 1024);
}

/// @brief Dispatcher copies retain exact immutable samples after their source dies.
TEST(QuadratureFormula, DispatcherCopiesRetainExistingSamples)
{
  std::unique_ptr<QF::PolytopeQuadratureFormula> copy;
  std::vector<Math::SpatialPoint> points;
  std::vector<Real> weights;
  {
    const QF::PolytopeQuadratureFormula source(8, G::Triangle);
    copy.reset(source.copy());
    ASSERT_EQ(copy->getSize(), source.getSize());
    EXPECT_EQ(&copy->getPoint(0), &source.getPoint(0));
    for (size_t q = 0; q < source.getSize(); ++q)
    {
      points.push_back(source.getPoint(q));
      weights.push_back(source.getWeight(q));
    }
  }
  ASSERT_EQ(copy->getSize(), points.size());
  for (size_t q = 0; q < points.size(); ++q)
  {
    EXPECT_EQ(copy->getWeight(q), weights[q]);
    EXPECT_EQ((copy->getPoint(q) - points[q]).norm(), 0);
  }
}
