/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <Rodin/Utility/OptionalReference.h>

#include <type_traits>
#include <utility>

#include <gtest/gtest.h>

namespace Rodin::Tests::Unit
{
  /// @brief Value whose identity and mutations expose reference aliasing.
  struct ReferenceValue
  {
      Integer value;  ///< Value shared through borrowing wrappers.
  };

  /// @brief Engagement and reset preserve the independently owned referent.
  TEST(OptionalReference, EngagementMutationAndReset)
  {
    ReferenceValue value{7};
    Utility::OptionalReference<ReferenceValue> reference;
    EXPECT_FALSE(reference.has_value());
    reference.emplace(value);
    ASSERT_TRUE(reference.has_value());
    EXPECT_EQ(&*reference, &value);
    EXPECT_EQ(reference.operator->(), &value);
    reference->value = 19;
    EXPECT_EQ(value.value, 19);
    (*reference).value += 2;
    EXPECT_EQ(value.value, 21);
    reference.reset();
    EXPECT_FALSE(reference.has_value());
    EXPECT_EQ(value.value, 21);
  }

  /// @brief Copies and moves preserve borrowing while rebinding changes only one wrapper.
  TEST(OptionalReference, CopyMoveAndRebinding)
  {
    ReferenceValue first{3}, second{5};
    Utility::OptionalReference<ReferenceValue> reference(first);
    auto copied(reference);
    EXPECT_EQ(&*copied, &first);
    Utility::OptionalReference<ReferenceValue> moved(std::move(copied));
    copied.reset();
    EXPECT_EQ(&*moved, &first);
    reference.emplace(second);
    EXPECT_EQ(reference.operator->(), &second);
    EXPECT_EQ(moved.operator->(), &first);
    moved->value = 11;
    EXPECT_EQ(first.value, 11);
    EXPECT_EQ(second.value, 5);
    reference = moved;
    EXPECT_EQ(&*reference, &first);
    reference.reset();
    EXPECT_EQ(&*moved, &first);
    EXPECT_EQ(moved->value, 11);
  }

  /// @brief Const referents retain their identity and const access types.
  TEST(OptionalReference, ConstReferent)
  {
    const ReferenceValue value{42};
    Utility::OptionalReference<const ReferenceValue> reference(value);
    static_assert(std::is_same_v<decltype(*reference), const ReferenceValue&>);
    static_assert(
      std::is_same_v<decltype(reference.operator->()), const ReferenceValue*>);
    ASSERT_TRUE(reference.has_value());
    EXPECT_EQ(reference.operator->(), &value);
    EXPECT_EQ(&*reference, &value);
    EXPECT_EQ(reference->value, 42);
  }
}
