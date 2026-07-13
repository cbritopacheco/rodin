#include <gtest/gtest.h>

#include <Rodin/Pair.h>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  //=== Variant Constructor ==================================================

  /// @brief Using deduction guide.
  TEST(Rodin_Pair, VariantConstructor)
  {
    // Using deduction guide.
    Pair p(100, 'a');
    EXPECT_EQ(p.first(), 100);
    EXPECT_EQ(p.second(), 'a');
  }

  //=== Construction =========================================================

  /// @brief Verifies construction for pair by checking exact expected values.
  TEST(Rodin_Pair, Construction)
  {
    Pair<int, std::string> p(42, "Answer");
    EXPECT_EQ(p.first(), 42);
    EXPECT_EQ(p.second(), "Answer");
  }

  //=== Copy and Move Semantics ==============================================

  /// @brief Verifies copy construction for pair by checking exact expected values, copy semantics.
  TEST(Rodin_Pair, CopyConstruction)
  {
    Pair<int, std::string> p1(1, "one");
    Pair<int, std::string> p2(p1);
    EXPECT_EQ(p2.first(), 1);
    EXPECT_EQ(p2.second(), "one");
  }

  /// @brief Verifies move construction for pair by checking exact expected values, move semantics.
  TEST(Rodin_Pair, MoveConstruction)
  {
    Pair<int, std::string> p1(2, "two");
    Pair<int, std::string> p2(std::move(p1));
    EXPECT_EQ(p2.first(), 2);
    EXPECT_EQ(p2.second(), "two");
  }

  /// @brief Verifies copy assignment for pair by checking exact expected values, copy semantics.
  TEST(Rodin_Pair, CopyAssignment)
  {
    Pair<int, char> p1(3, 'c');
    Pair<int, char> p2(0, 'z');
    p2 = p1;
    EXPECT_EQ(p2.first(), 3);
    EXPECT_EQ(p2.second(), 'c');
  }

  /// @brief Verifies move assignment for pair by checking exact expected values, move semantics.
  TEST(Rodin_Pair, MoveAssignment)
  {
    Pair<int, char> p1(7, 'x');
    Pair<int, char> p2(0, 'y');
    p2 = std::move(p1);
    EXPECT_EQ(p2.first(), 7);
    EXPECT_EQ(p2.second(), 'x');
  }

  //=== Accessors ============================================================

  /// @brief Verifies get accessors for pair by checking exact expected values.
  TEST(Rodin_Pair, GetAccessors)
  {
    Pair<int, double> p(5, 3.14);
    // Using Tuple's get<Index>() inherited from Pair.
    EXPECT_EQ(p.get<0>(), 5);
    EXPECT_EQ(p.get<1>(), 3.14);
  }

  /// @brief Verifies first second const for pair by checking exact expected values.
  TEST(Rodin_Pair, FirstSecondConst)
  {
    const Pair<int, char> p(9, 'q');
    EXPECT_EQ(p.first(), 9);
    EXPECT_EQ(p.second(), 'q');
  }
}
