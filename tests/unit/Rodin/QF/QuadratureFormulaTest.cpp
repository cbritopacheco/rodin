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
#include "Rodin/QF/GrundmannMoller.h"
#include "Rodin/Adaptation/SWIFT/QualityCovering.h"
#include "Rodin/QF/GaussLobatto.h"
#include "Rodin/QF/QuadratureFormula.h"
#include "QuadratureSelection.h"
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
      inline static std::atomic<size_t> destructions{0};
      inline static std::atomic<const Math::SpatialPoint*> point{nullptr};
      inline static std::atomic<bool> fail{false};
      inline static std::atomic<bool> failSize{false};
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
        point = &m_point;
      }

      FactoryFormula(const FactoryFormula&) = delete;
      FactoryFormula(FactoryFormula&&) = delete;

      ~FactoryFormula() override
      {
        ++destructions;
      }

      /// @returns One sample, recording the generator read.
      size_t getSize() const override
      {
        ++reads;
        if (failSize.load())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Interrupted test formula publication." << Alert::Raise;
        }
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
        return m_point;
      }

    private:
      Math::SpatialPoint m_point;
  };

  /// @brief Strided sample storage whose pointwise accessors record dispatch.
  class SampleViewFormula final : public QF::QuadratureFormulaBase
  {
    public:
      inline static size_t reads = 0;
      inline static size_t destructions = 0;
      inline static bool fail = false;
      inline static const Math::SpatialPoint* points = nullptr;

      explicit SampleViewFormula([[maybe_unused]] size_t tag)
        : m_points{Math::SpatialPoint{Real(0.25)}, Math::SpatialPoint{Real(0.75)}},
          m_weights{Real(0.25), Real(-100), Real(0.75), Real(-200)}
      {
        points = m_points.data();
      }

      SampleViewFormula(const SampleViewFormula&) = delete;
      SampleViewFormula(SampleViewFormula&&) = delete;

      ~SampleViewFormula() override
      {
        ++destructions;
      }

      size_t getSize() const override
      {
        return m_points.size();
      }

      Real getWeight(size_t index) const override
      {
        ++reads;
        return m_weights[2 * index];
      }

      const Math::SpatialPoint& getPoint(size_t index) const override
      {
        ++reads;
        return m_points[index];
      }

      Optional<Samples> getSamples() const override
      {
        if (fail)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Interrupted sample-view publication." << Alert::Raise;
        }
        return Samples{m_points, m_weights.data(), 2};
      }

    private:
      std::array<Math::SpatialPoint, 2> m_points;
      std::array<Real, 4> m_weights;
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
  EXPECT_EQ(&QF::QuadratureFormula::get(count, geometry),
    &QF::QuadratureFormula::get(count, geometry));
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

/// @brief Publication retains the concrete sample storage without copying or moving it.
TEST(QuadratureFormula, PublicationRetainsConcreteSamples)
{
  const size_t destroyed = FactoryFormula::destructions.load();
  const auto& formula = QF::QuadratureFormula::get<FactoryFormula>(size_t{101});
  EXPECT_EQ(&formula.getPoint(0), FactoryFormula::point.load());
  EXPECT_EQ(FactoryFormula::destructions.load(), destroyed);
  const size_t constructions = FactoryFormula::constructions.load();
  Variational::RealP1Element element(G::Segment);
  const auto& table = element.getTabulation(formula);
  const size_t reads = FactoryFormula::reads.load();
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
  const size_t destroyed = FactoryFormula::destructions.load();
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
  EXPECT_EQ(FactoryFormula::destructions.load(), destroyed + ThreadCount - 1);
  ASSERT_NE(formulas[0], nullptr);
  for (const auto* formula : formulas)
  {
    EXPECT_EQ(formula, formulas[0]);
    EXPECT_EQ(formula->getPoint(0)(0), Real(tag) / 1024);
  }
  EXPECT_EQ(formulas[0], &QF::QuadratureFormula::get<FactoryFormula>(tag));
}

/// @brief Interrupted publication destroys its candidate and allows a later retry.
TEST(QuadratureFormula, RetriesFailedPublication)
{
  const size_t tag = FactoryFormula::nextTag.fetch_add(1);
  const size_t before = FactoryFormula::constructions.load();
  const size_t destroyed = FactoryFormula::destructions.load();
  FactoryFormula::failSize = true;
  EXPECT_THROW(QF::QuadratureFormula::get<FactoryFormula>(tag), Alert::Exception);
  EXPECT_EQ(FactoryFormula::destructions.load(), destroyed + 1);
  FactoryFormula::failSize = false;
  const auto& formula = QF::QuadratureFormula::get<FactoryFormula>(tag);
  EXPECT_EQ(FactoryFormula::constructions.load(), before + 2);
  EXPECT_EQ(FactoryFormula::destructions.load(), destroyed + 1);
  EXPECT_EQ(formula.getPoint(0)(0), Real(tag) / 1024);
}

/// @brief Concrete family copies retain exact immutable samples after their source dies.
TEST(QuadratureFormula, ConcreteCopiesRetainExistingSamples)
{
  std::unique_ptr<const QF::QuadratureFormulaBase> copy;
  std::vector<Math::SpatialPoint> points;
  std::vector<Real> weights;
  {
    const auto source = Tests::QF::makeQuadratureFormula(8, G::Triangle);
    copy = Tests::QF::copyQuadratureFormula(*source);
    ASSERT_EQ(copy->getSize(), source->getSize());
    EXPECT_EQ((copy->getPoint(0) - source->getPoint(0)).norm(), 0);
    for (size_t q = 0; q < source->getSize(); ++q)
    {
      points.push_back(source->getPoint(q));
      weights.push_back(source->getWeight(q));
    }
  }
  ASSERT_EQ(copy->getSize(), points.size());
  for (size_t q = 0; q < points.size(); ++q)
  {
    EXPECT_EQ(copy->getWeight(q), weights[q]);
    EXPECT_EQ((copy->getPoint(q) - points[q]).norm(), 0);
  }
}

/// @brief Concrete owned and persistent construction retain identical rules across dispatch boundaries.
TEST(QuadratureFormula, OwnedAndPublishedDispatchAgree)
{
  const auto build = &Tests::QF::makeQuadratureFormula;
  for (const auto geometry : Geometry::Polytope::Types)
  {
    for (const size_t order : {0u, 1u, 2u, 10u, 11u, 15u, 16u, 50u, 51u})
    {
      SCOPED_TRACE(static_cast<size_t>(geometry));
      SCOPED_TRACE(order);
      const auto owned = build(order, geometry);
      std::unique_ptr<const QF::QuadratureFormulaBase> copy;
      {
        const auto source = Tests::QF::makeQuadratureFormula(order, geometry);
        copy = Tests::QF::copyQuadratureFormula(*source);
        EXPECT_EQ((copy->getPoint(0) - source->getPoint(0)).norm(), 0);
      }
      const auto& published = QF::QuadratureFormula::get(order, geometry);
      ASSERT_EQ(owned->getSize(), published.getSize());
      ASSERT_EQ(copy->getSize(), published.getSize());
      for (size_t q = 0; q < published.getSize(); ++q)
      {
        EXPECT_EQ(owned->getWeight(q), published.getWeight(q));
        EXPECT_EQ(copy->getWeight(q), published.getWeight(q));
        const auto& point = published.getPoint(q);
        ASSERT_EQ(owned->getPoint(q).size(), point.size());
        ASSERT_EQ(copy->getPoint(q).size(), point.size());
        for (size_t d = 0; d < point.size(); ++d)
        {
          EXPECT_EQ(owned->getPoint(q)(d), point(d));
          EXPECT_EQ(copy->getPoint(q)(d), point(d));
        }
      }
    }
  }
}

/// @brief Retained strided views avoid both sample copies and pointwise dispatch.
TEST(QuadratureFormula, SampleViewsAvoidCopyAndDispatch)
{
  SampleViewFormula::fail = true;
  EXPECT_THROW(
    QF::QuadratureFormula::get<SampleViewFormula>(size_t{1}), Alert::Exception);
  EXPECT_EQ(SampleViewFormula::destructions, 1);
  SampleViewFormula::fail = false;
  const auto& formula = QF::QuadratureFormula::get<SampleViewFormula>(size_t{1});
  const auto samples = formula.getSamples();
  ASSERT_TRUE(samples);
  EXPECT_EQ(samples->points.data(), SampleViewFormula::points);
  ASSERT_EQ(formula.getSize(), 2);
  EXPECT_EQ(&formula.getPoint(0), SampleViewFormula::points);
  EXPECT_EQ(formula.getPoint(0)(0), Real(0.25));
  EXPECT_EQ(formula.getPoint(1)(0), Real(0.75));
  EXPECT_EQ(formula.getWeight(0), Real(0.25));
  EXPECT_EQ(formula.getWeight(1), Real(0.75));
  Variational::RealP1Element element(G::Segment);
  const auto& table = element.getTabulation(formula);
  EXPECT_EQ(table.getBasis(0, 1), Real(0.25));
  EXPECT_EQ(table.getBasis(1, 1), Real(0.75));
  EXPECT_EQ(SampleViewFormula::reads, 0);
  EXPECT_EQ(SampleViewFormula::destructions, 1);
}

/// @brief Uniform, interleaved and contiguous views preserve every concrete family.
TEST(QuadratureFormula, SampleViewsPreserveConcreteFamilies)
{
  const auto verify = [](const QF::QuadratureFormulaBase& source,
                        const QF::QuadratureFormula& published) {
    ASSERT_TRUE(source.getSamples());
    ASSERT_TRUE(published.getSamples());
    ASSERT_EQ(published.getSize(), source.getSize());
    for (size_t index = 0; index < source.getSize(); ++index)
    {
      EXPECT_EQ(published.getWeight(index), source.getWeight(index));
      EXPECT_EQ((published.getPoint(index) - source.getPoint(index)).norm(), 0);
    }
  };
  for (auto geometry : {G::Point, G::Segment, G::Triangle, G::Quadrilateral,
         G::Tetrahedron, G::Pyramid, G::Hexahedron, G::Wedge})
  {
    const QF::Centroid centroid(geometry);
    verify(centroid, QF::QuadratureFormula::get<QF::Centroid>(geometry));
    const auto dispatcher = Tests::QF::makeQuadratureFormula(4, geometry);
    verify(*dispatcher, QF::QuadratureFormula::get(4, geometry));
    const QF::GaussLobatto lobatto(geometry, 3);
    verify(lobatto, QF::GaussLobatto::get(geometry, 3));
    const Adaptation::SWIFT::WitnessLattice lattice(geometry, 2);
    verify(lattice, Adaptation::SWIFT::WitnessLattice::get(geometry, 2));
    const Adaptation::SWIFT::QualityCovering covering(geometry, 2);
    verify(covering, Adaptation::SWIFT::QualityCovering::get(geometry, 2));
  }
  for (auto geometry : {G::Segment, G::Triangle, G::Tetrahedron})
  {
    const QF::GrundmannMoller formula(2, geometry);
    verify(formula, QF::QuadratureFormula::get<QF::GrundmannMoller>(size_t{2}, geometry));
  }
}

/// @brief Default requests publish the same immutable formula across producer threads.
TEST(QuadratureFormula, DefaultDispatchSurvivesThreadsAndLookupEviction)
{
  constexpr size_t ThreadCount = 4;
  std::array<const QF::QuadratureFormula*, ThreadCount> formulas{};
  std::array<std::thread, ThreadCount> threads;
  std::barrier gate(static_cast<std::ptrdiff_t>(ThreadCount));
  for (size_t i = 0; i < ThreadCount; ++i)
  {
    threads[i] = std::thread([&, i] {
      gate.arrive_and_wait();
      formulas[i] = &QF::QuadratureFormula::get(18, G::Wedge);
      for (size_t order = 1; order <= 12; ++order)
        QF::QuadratureFormula::get(order, G::Triangle);
    });
  }
  for (auto& thread : threads)
    thread.join();
  const auto& expected = QF::QuadratureFormula::get(18, G::Wedge);
  for (const auto* formula : formulas)
  {
    ASSERT_EQ(formula, &expected);
    EXPECT_EQ(formula->getWeight(0), expected.getWeight(0));
    EXPECT_EQ(&formula->getPoint(0), &expected.getPoint(0));
  }
}

/// @brief Unsupported geometries fail before publishing a default lookup entry.
TEST(QuadratureFormula, DefaultDispatchRejectsInvalidGeometry)
{
  EXPECT_THROW(QF::QuadratureFormula::get(2, static_cast<G>(255)), Alert::Exception);
  EXPECT_GT(QF::QuadratureFormula::get(2, G::Triangle).getSize(), 0u);
}
