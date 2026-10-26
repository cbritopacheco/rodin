/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Distance/Base.h"

namespace Rodin::Tests::Unit
{
  class DistanceSelection final : public Distance::Base<DistanceSelection>
  {};
  using Geometry::Attribute;

  /// @brief Region sets are owned copies and setters preserve the derived chaining type.
  TEST(DistanceBase, ExplicitSetsOwnTheirAttributes)
  {
    DistanceSelection model;
    EXPECT_TRUE(model.getInterior().empty());
    EXPECT_TRUE(model.getInterface().empty());
    FlatSet<Attribute> interior{Attribute(3), Attribute(7)};
    FlatSet<Attribute> interface{Attribute(11), Attribute(13)};
    const auto& interiorView = interior;
    const auto& interfaceView = interface;
    EXPECT_EQ(&model.setInterior(interiorView).setInterface(interfaceView), &model);
    interior.clear();
    interface.clear();
    const FlatSet<Attribute> expectedInterior{Attribute(3), Attribute(7)};
    const FlatSet<Attribute> expectedInterface{Attribute(11), Attribute(13)};
    const auto& readonly = model;
    EXPECT_EQ(readonly.getInterior(), expectedInterior);
    EXPECT_EQ(readonly.getInterface(), expectedInterface);
  }

  /// @brief Replacement, deduplication, self-assignment and clearing preserve separate selectors.
  TEST(DistanceBase, ReplacingAndClearingSelectorsAreIndependent)
  {
    DistanceSelection model;
    model.setInterior(Attribute(7), Attribute(3), Attribute(7))
      .setInterface(Attribute(11), Attribute(11));
    const FlatSet<Attribute> originalInterior{Attribute(3), Attribute(7)};
    const FlatSet<Attribute> originalInterface{Attribute(11)};
    EXPECT_EQ(model.getInterior(), originalInterior);
    EXPECT_EQ(model.getInterface(), originalInterface);
    model.setInterior(model.getInterior()).setInterface(model.getInterface());
    EXPECT_EQ(model.getInterior(), originalInterior);
    EXPECT_EQ(model.getInterface(), originalInterface);
    model.setInterior(Attribute(19));
    const FlatSet<Attribute> replacement{Attribute(19)};
    EXPECT_EQ(model.getInterior(), replacement);
    EXPECT_EQ(model.getInterface(), originalInterface);
    const FlatSet<Attribute> empty;
    model.setInterior(empty);
    EXPECT_TRUE(model.getInterior().empty());
    EXPECT_EQ(model.getInterface(), originalInterface);
    model.setInterface(empty);
    EXPECT_TRUE(model.getInterface().empty());
  }

  /// @brief Copying a configured distance model separates subsequent selector changes.
  TEST(DistanceBase, ModelCopiesRetainIndependentRegionSelections)
  {
    DistanceSelection original;
    original.setInterior(Attribute(3), Attribute(5)).setInterface(Attribute(7));
    auto copied = original;
    original.setInterior(Attribute(17)).setInterface(Attribute(19));
    const FlatSet<Attribute> expectedInterior{Attribute(3), Attribute(5)};
    const FlatSet<Attribute> expectedInterface{Attribute(7)};
    EXPECT_EQ(copied.getInterior(), expectedInterior);
    EXPECT_EQ(copied.getInterface(), expectedInterface);
    copied.setInterior(Attribute(23)).setInterface(Attribute(29));
    EXPECT_TRUE(original.getInterior().contains(Attribute(17)));
    EXPECT_FALSE(original.getInterior().contains(Attribute(23)));
    EXPECT_TRUE(original.getInterface().contains(Attribute(19)));
    EXPECT_FALSE(original.getInterface().contains(Attribute(29)));
  }
}
