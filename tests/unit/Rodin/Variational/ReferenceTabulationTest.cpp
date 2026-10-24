/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Certifies element-owned filling and reference-cache publication.
 */
#include <gtest/gtest.h>
#include <utility>
#include <vector>
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Variational/ReferenceTabulationCache.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Exposes exact fill counts for a linear segment's bulk operation.
   *
   * The basis is @f$(1-x,x)@f$. The public bulk operation fills caller-owned numerical storage directly,
   * so dispatch and cache publication can be checked independently of a
   * particular optimized element algorithm.
   */
  class CountingTabulationElement
  {
    public:
      using Store =
        Variational::ReferenceTabulationCache<CountingTabulationElement, Real>;

      /// @brief One linear segment basis function.
      struct Basis
      {
        size_t local;

        /**
         * @brief Selects a basis function.
         * @param index Zero selects the left endpoint, one the right endpoint.
         */
        explicit Basis(size_t index)
          : local(index)
        {}

        /**
         * @brief Evaluates the selected linear basis.
         * @param point Reference segment coordinate.
         * @returns The selected basis value.
         */
        Real operator()(const Math::SpatialPoint& point) const
        {
          return local ? point(0) : 1 - point(0);
        }

        /**
         * @brief Gets the constant first derivative.
         * @tparam Order Derivative order, required to be one in this fixture.
         * @param coordinate Unused coordinate: the segment has only one.
         * @returns A callable returning minus one or one.
         */
        template <size_t Order>
        auto getDerivative([[maybe_unused]] size_t coordinate) const
        {
          static_assert(Order == 1);
          return [index = local]([[maybe_unused]] const Math::SpatialPoint& point) {
            return index ? Real(1) : Real(-1);
          };
        }
      };

      /// @brief Constructs an element with no fills and no requested failure.
      CountingTabulationElement()
        : m_fills(0),
          m_fail(false)
      {}

      /**
       * @brief Gets the segment reference geometry.
       * @returns The segment geometry.
       */
      auto getGeometry() const
      {
        return Geometry::Polytope::Type::Segment;
      }

      /**
       * @brief Gets the basis size.
       * @returns Two endpoint functions.
       */
      size_t getCount() const
      {
        return 2;
      }

      /**
       * @brief Gets one basis descriptor by value.
       * @param index Local endpoint index.
       * @returns The selected immutable descriptor.
       */
      Basis getBasis(size_t index) const
      {
        return Basis(index);
      }

      /**
       * @brief Gets the number of attempted bulk fills.
       * @returns The fill count, including interrupted fills.
       */
      size_t getFills() const
      {
        return m_fills;
      }

      /**
       * @brief Selects whether filling should fail before cache publication.
       * @param fail Whether the test operation raises an Alert exception.
       */
      void setFail(bool fail)
      {
        m_fail = fail;
      }

      /**
       * @brief Fills the table and optionally interrupts its publication.
       * @param table Numerical storage sized and filled by the element.
       * @param formula Reference sample collection.
       */
      void tabulate(Store::Table& table, const QF::QuadratureFormulaBase& formula) const
      {
        ++m_fills;
        table.resize(formula.getSize(), getCount(), 1);
        for (size_t q = 0; q < table.getSize(); ++q)
        {
          const auto& point = formula.getPoint(q);
          for (size_t a = 0; a < table.getCount(); ++a)
          {
            const auto basis = getBasis(a);
            table.getValues(q).data[a] = basis(point);
            table.getGradients(q).data[a] = basis.template getDerivative<1>(0)(point);
          }
        }
        if (m_fail)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Deliberately interrupted reference-table publication" << Alert::Raise;
        }
      }

    private:
      mutable size_t m_fills;
      bool m_fail;
  };

  TEST(ReferenceTabulation, BulkFillRunsOncePerMiss)
  {
    using Store = CountingTabulationElement::Store;
    CountingTabulationElement element;
    const auto& formula = QF::QuadratureFormula::get<QF::PolytopeQuadratureFormula>(
      size_t{3}, element.getGeometry());
    const auto& table = Store::get(element, formula);
    EXPECT_EQ(element.getFills(), 1u);
    EXPECT_EQ(&table, &Store::get(element, formula));
    EXPECT_EQ(element.getFills(), 1u);
    const auto found = Store::find(element, formula);
    ASSERT_TRUE(found);
    EXPECT_EQ(&found->get(), &table);
    for (size_t q = 0; q < formula.getSize(); ++q)
    {
      EXPECT_EQ(table.getBasis(q, 0), 1 - formula.getPoint(q)(0));
      EXPECT_EQ(table.getBasis(q, 1), formula.getPoint(q)(0));
      EXPECT_EQ(table.getDerivative<1>(q, 0, 0), -1);
      EXPECT_EQ(table.getDerivative<1>(q, 1, 0), 1);
    }
  }

  TEST(ReferenceTabulation, EvictionRebuildsThroughElement)
  {
    using Store = CountingTabulationElement::Store;
    CountingTabulationElement element;
    std::vector<const QF::QuadratureFormula*> formulas;
    constexpr size_t CacheCapacity = 8; // The store's documented capacity.
    constexpr size_t FormulaCount = CacheCapacity + 1;
    formulas.reserve(FormulaCount);
    for (size_t i = 0; i < FormulaCount; ++i)
    {
      formulas.push_back(
        &QF::PolytopeQuadratureFormula::get(10 + i, element.getGeometry()));
      Store::get(element, *formulas.back());
    }
    EXPECT_EQ(element.getFills(), FormulaCount);
    EXPECT_FALSE(Store::find(element, *formulas.front()));
    for (size_t i = 1; i < FormulaCount; ++i)
      EXPECT_TRUE(Store::find(element, *formulas[i]));
    const auto& rebuilt = Store::get(element, *formulas.front());
    EXPECT_EQ(element.getFills(), FormulaCount + 1);
    EXPECT_EQ(rebuilt.getDerivative<1>(0, 0, 0), -1);
  }

  TEST(ReferenceTabulation, InterruptedFillIsNotPublished)
  {
    using Store = CountingTabulationElement::Store;
    CountingTabulationElement element;
    const auto& formula = QF::QuadratureFormula::get<QF::PolytopeQuadratureFormula>(
      size_t{40}, element.getGeometry());
    element.setFail(true);
    EXPECT_THROW(Store::get(element, formula), Alert::Exception);
    EXPECT_EQ(element.getFills(), 1u);
    EXPECT_FALSE(Store::find(element, formula));
    element.setFail(false);
    const auto& table = Store::get(element, formula);
    EXPECT_EQ(element.getFills(), 2u);
    EXPECT_TRUE(Store::find(element, formula));
    EXPECT_EQ(table.getDerivative<1>(0, 1, 0), 1);
  }

  /**
   * @brief Failed replacement retires only its selected slot's previous key.
   *
   * Filling a complete cache establishes the round-robin eviction order without
   * using storage addresses. A deliberately interrupted replacement must expose
   * neither the evicted key nor the new key. Surviving entries retain their
   * exact samples. Retrying consumes the next slot, as for every insertion
   * attempt; failed fills do not rewind the existing round-robin policy.
   */
  TEST(ReferenceTabulation, InterruptedEvictionRetiresPreviousIdentity)
  {
    using Store = CountingTabulationElement::Store;
    constexpr size_t CacheCapacity = 8; // The store's documented capacity.
    CountingTabulationElement element;
    std::vector<const QF::QuadratureFormula*> formulas;
    std::vector<Store::Table> snapshots;
    formulas.reserve(CacheCapacity);
    snapshots.reserve(CacheCapacity);
    for (size_t i = 0; i < CacheCapacity; ++i)
    {
      formulas.push_back(
        &QF::PolytopeQuadratureFormula::get(70 + i, element.getGeometry()));
      snapshots.push_back(Store::get(element, *formulas.back()));
    }
    const auto& replacement = QF::QuadratureFormula::get<QF::PolytopeQuadratureFormula>(
      size_t{80}, element.getGeometry());
    element.setFail(true);
    EXPECT_THROW(Store::get(element, replacement), Alert::Exception);
    EXPECT_EQ(element.getFills(), CacheCapacity + 1);
    EXPECT_FALSE(Store::find(element, *formulas.front()));
    EXPECT_FALSE(Store::find(element, replacement));
    for (size_t i = 1; i < CacheCapacity; ++i)
    {
      const auto surviving = Store::find(element, *formulas[i]);
      ASSERT_TRUE(surviving);
      EXPECT_EQ(surviving->get().getValues(), snapshots[i].getValues());
      EXPECT_EQ(surviving->get().getDerivatives(), snapshots[i].getDerivatives());
    }
    element.setFail(false);
    const auto& retried = Store::get(element, replacement);
    const auto published = Store::find(element, replacement);
    ASSERT_TRUE(published);
    EXPECT_EQ(&published->get(), &retried);
    EXPECT_EQ(retried.getSize(), replacement.getSize());
    EXPECT_EQ(element.getFills(), CacheCapacity + 2);
    EXPECT_FALSE(Store::find(element, *formulas[1]));
    for (size_t i = 2; i < CacheCapacity; ++i)
      EXPECT_TRUE(Store::find(element, *formulas[i]));
  }

  TEST(ReferenceTabulation, WarmFormulaReusesPreparedTable)
  {
    using Store = CountingTabulationElement::Store;
    CountingTabulationElement element;
    const auto& formula = QF::PolytopeQuadratureFormula::get(41, element.getGeometry());
    const auto& table = Store::get(element, formula);
    const auto values = table.getValues();
    EXPECT_EQ(&table, &Store::get(element, formula));
    EXPECT_EQ(table.getValues(), values);
    EXPECT_EQ(element.getFills(), 1u);
    const auto& replacement =
      QF::PolytopeQuadratureFormula::get(42, element.getGeometry());
    const auto& replaced = Store::get(element, replacement);
    EXPECT_EQ(replaced.getSize(), replacement.getSize());
    EXPECT_EQ(element.getFills(), 2u);
  }

  TEST(ReferenceTabulation, LocalGeneratorTabulationOwnsItsValues)
  {
    CountingTabulationElement element;
    CountingTabulationElement::Store::Table table;
    {
      QF::GaussLegendre local(Geometry::Polytope::Type::Segment, 2);
      element.tabulate(table, local);
    }
    ASSERT_EQ(table.getSize(), 2u);
    EXPECT_EQ(table.getDerivative<1>(0, 0, 0), -1);
    EXPECT_EQ(table.getDerivative<1>(0, 1, 0), 1);
  }
}

using namespace Rodin;
using namespace Rodin::Variational;

namespace
{
  template <class Table>
  concept HasSecondReferenceDerivative =
    requires(const Table& table) { table.template getDerivative<2>(0, 0, 0); };
  static_assert(!HasSecondReferenceDerivative<ReferenceTabulation<Real>>);
}

TEST(ReferenceTabulation, ControlledStorageExtentsAndOverflow)
{
  ReferenceTabulation<Real> table;
  table.resize(2, 3, 2);
  ASSERT_EQ(table.getValues().size(), 6u);
  ASSERT_EQ(table.getDerivatives().size(), 12u);
  table.getValues(1).data[2] = 7;
  table.getGradients(1).data[5] = 9;
  EXPECT_EQ(table.getBasis(1, 2), 7);
  EXPECT_EQ(table.getDerivative<1>(1, 2, 1), 9);
  const auto saved = table;
  EXPECT_THROW(table.resize(std::numeric_limits<size_t>::max(), 2, 1), Alert::Exception);
  EXPECT_EQ(table.getValues(), saved.getValues());
  EXPECT_EQ(table.getDerivatives(), saved.getDerivatives());
  EXPECT_EQ(table.getSize(), 2u);
  table.resize(1, 4, 0);
  EXPECT_EQ(table.getValues().size(), 4u);
  EXPECT_TRUE(table.getDerivatives().empty());
  EXPECT_TRUE(table.getGradients(0).data.empty());
}

TEST(ReferenceTabulation, MovingPreservesBothStorageContracts)
{
  ReferenceTabulation<Real> source;
  source.resize(2, 3, 2);
  source.getValues(1).data[2] = 7;
  const auto saved = source;
  auto moved = std::move(source);
  EXPECT_EQ(moved.getValues(), saved.getValues());
  EXPECT_EQ(source.getSize(), 0u);
  EXPECT_EQ(source.getCount(), 0u);
  EXPECT_EQ(source.getDimension(), 0u);
  EXPECT_TRUE(source.getValues().empty());
  EXPECT_TRUE(source.getDerivatives().empty());
  source = std::move(moved);
  EXPECT_EQ(source.getValues(), saved.getValues());
  EXPECT_EQ(moved.getSize(), 0u);
  EXPECT_TRUE(moved.getValues().empty());
  moved = source;
  source.resize(0, 4, 3);
  EXPECT_EQ(moved.getValues(), saved.getValues());
  EXPECT_EQ(moved.getDerivatives(), saved.getDerivatives());
  auto& alias = moved;
  moved = std::move(alias);
  EXPECT_EQ(moved.getValues(), saved.getValues());
}
