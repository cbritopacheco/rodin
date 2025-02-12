#include <gtest/gtest.h>

#include <Rodin/Tuple.h>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  template <typename T>
  struct IsIntegral
  {
    static constexpr bool Value = std::is_integral_v<T>;
  };

  template <typename T>
  struct AlwaysTrue
  {
    static constexpr bool Value = true;
  };

  template <typename T>
  struct AlwaysFalse
  {
    static constexpr bool Value = false;
  };

  //=== Variant Constructors ==================================================

  TEST(Rodin_Tuple, VariantConstructor)
  {
    {
      Tuple t;
      EXPECT_EQ(t.size(), 0);
    }

    {
      Tuple t{1};
      EXPECT_EQ(t.size(), 1);
      EXPECT_EQ(t.get<0>(), 1);
    }

    {
      Tuple t{'a'};
      EXPECT_EQ(t.size(), 1);
      EXPECT_EQ(t.get<0>(), 'a');
    }

    {
      Tuple t{1, 2};
      EXPECT_EQ(t.size(), 2);
      EXPECT_EQ(t.get<0>(), 1);
      EXPECT_EQ(t.get<1>(), 2);
    }

    {
      Tuple t{'a', 'b'};
      EXPECT_EQ(t.size(), 2);
      EXPECT_EQ(t.get<0>(), 'a');
      EXPECT_EQ(t.get<1>(), 'b');
    }

    {
      Tuple t{1, 'a', 1.0};
      EXPECT_EQ(t.size(), 3);
      EXPECT_EQ(t.get<0>(), 1);
      EXPECT_EQ(t.get<1>(), 'a');
      EXPECT_EQ(t.get<2>(), 1.0);
    }
  }

  //=== Copy and Move Semantics ===============================================

  TEST(Rodin_Tuple, CopyConstructor)
  {
    Tuple t1{1, 'a', 1.0};
    Tuple t2(t1);

    EXPECT_EQ(t2.size(), 3);
    EXPECT_EQ(t2.get<0>(), 1);
    EXPECT_EQ(t2.get<1>(), 'a');
    EXPECT_EQ(t2.get<2>(), 1.0);
  }

  TEST(Rodin_Tuple, CopyAssignment)
  {
    Tuple t1{1, 2, 3};
    Tuple t2{4, 5, 6};
    t2 = t1;

    EXPECT_EQ(t2.size(), 3);
    EXPECT_EQ(t2.get<0>(), 1);
    EXPECT_EQ(t2.get<1>(), 2);
    EXPECT_EQ(t2.get<2>(), 3);
  }

  TEST(Rodin_Tuple, MoveConstructor)
  {
    Tuple t1{7, 8, 9};
    Tuple t2(std::move(t1));

    EXPECT_EQ(t2.size(), 3);
    EXPECT_EQ(t2.get<0>(), 7);
    EXPECT_EQ(t2.get<1>(), 8);
    EXPECT_EQ(t2.get<2>(), 9);
  }

  TEST(Rodin_Tuple, MoveAssignment)
  {
    Tuple t1{7, 8, 9};
    Tuple t2{0, 0, 0};
    t2 = std::move(t1);

    EXPECT_EQ(t2.get<0>(), 7);
    EXPECT_EQ(t2.get<1>(), 8);
    EXPECT_EQ(t2.get<2>(), 9);
  }

  //=== Concatenation =========================================================

  TEST(Rodin_Tuple, ConcatenateNonEmpty)
  {
    Tuple t1{1, 'a', 1.0};
    Tuple t2{2, 'b', 2.0};

    auto r = t1.concatenate(t2);

    EXPECT_EQ(r.size(), t1.size() + t2.size());
    EXPECT_EQ(r.get<0>(), 1);
    EXPECT_EQ(r.get<1>(), 'a');
    EXPECT_EQ(r.get<2>(), 1.0);
    EXPECT_EQ(r.get<3>(), 2);
    EXPECT_EQ(r.get<4>(), 'b');
    EXPECT_EQ(r.get<5>(), 2.0);
  }

  TEST(Rodin_Tuple, ConcatenateWithEmpty)
  {
    Tuple t1{1, 2, 3};
    Tuple empty;

    auto concat1 = t1.concatenate(empty);
    auto concat2 = empty.concatenate(t1);

    EXPECT_EQ(concat1.size(), t1.size());
    EXPECT_EQ(concat2.size(), t1.size());
    EXPECT_EQ(concat1.get<0>(), 1);
    EXPECT_EQ(concat1.get<1>(), 2);
    EXPECT_EQ(concat1.get<2>(), 3);
    EXPECT_EQ(concat2.get<0>(), 1);
    EXPECT_EQ(concat2.get<1>(), 2);
    EXPECT_EQ(concat2.get<2>(), 3);
  }

  TEST(Rodin_Tuple, ConcatenateChain)
  {
    Tuple t1{1};
    Tuple t2{2};
    Tuple t3{3};

    auto concatenated = t1.concatenate(t2).concatenate(t3);
    EXPECT_EQ(concatenated.size(), 3);
    EXPECT_EQ(concatenated.get<0>(), 1);
    EXPECT_EQ(concatenated.get<1>(), 2);
    EXPECT_EQ(concatenated.get<2>(), 3);
  }

  //=== Zip ===================================================================

  TEST(Rodin_Tuple, ZipTwoTuples)
  {
    Tuple t1{1, 'a', 1.0};
    Tuple t2{2, 'b', 2.0};

    auto r = t1.zip(t2);
    EXPECT_EQ(r.size(), t1.size());
    EXPECT_EQ(r.size(), t2.size());
  }

  TEST(Rodin_Tuple, ZipMultipleTuples)
  {
    Tuple t1{1, 2, 3};
    Tuple t2{10, 20, 30};
    Tuple t3{100, 200, 300};

    auto zipped = t1.zip([](auto a, auto b, auto c) { return a + b + c; }, t2, t3);
    EXPECT_EQ(zipped.size(), 3);
    EXPECT_EQ(zipped.get<0>(), 1 + 10 + 100);
    EXPECT_EQ(zipped.get<1>(), 2 + 20 + 200);
    EXPECT_EQ(zipped.get<2>(), 3 + 30 + 300);
  }

  TEST(Rodin_Tuple, ZipCustomLambda)
  {
    Tuple t1{2, 3, 4};
    Tuple t2{5, 6, 7};

    // Multiply corresponding elements using a custom lambda.
    auto zipped = t1.zip([](auto a, auto b) { return a * b; }, t2);
    EXPECT_EQ(zipped.size(), 3);
    EXPECT_EQ(zipped.get<0>(), 2 * 5);
    EXPECT_EQ(zipped.get<1>(), 3 * 6);
    EXPECT_EQ(zipped.get<2>(), 4 * 7);
  }

  //=== Reduce ================================================================

  TEST(Rodin_Tuple, ReduceSum)
  {
    {
      Tuple t{1, 2, 3, 4};
      auto result = t.reduce([](auto a, auto b) { return a + b; });
      EXPECT_EQ(result, 10);
    }

    {
      Tuple t{3, 10, 2, 5};
      auto result = t.reduce([](auto a, auto b) { return a + b; });
      EXPECT_EQ(result, 20);
    }
  }

  TEST(Rodin_Tuple, ReduceProduct)
  {
    Tuple t{2, 3, 4};
    auto result = t.reduce([](auto a, auto b) { return a * b; });
    EXPECT_EQ(result, 24);
  }

  TEST(Rodin_Tuple, ReduceTwoElements)
  {
    Tuple t{5, 10};
    auto result = t.reduce([](auto a, auto b) { return a - b; });
    EXPECT_EQ(result, 5 - 10);
  }

  //=== Map ===================================================================

  TEST(Rodin_Tuple, MapTransform)
  {
    Tuple t{1, 2, 3};
    auto mapped = t.map([](auto x) { return x * 10; });
    EXPECT_EQ(mapped.get<0>(), 10);
    EXPECT_EQ(mapped.get<1>(), 20);
    EXPECT_EQ(mapped.get<2>(), 30);
  }

  TEST(Rodin_Tuple, MapTransformConst)
  {
    const Tuple t{2, 3, 4};
    auto mapped = t.map([](auto x) { return x + 5; });
    EXPECT_EQ(mapped.get<0>(), 7);
    EXPECT_EQ(mapped.get<1>(), 8);
    EXPECT_EQ(mapped.get<2>(), 9);
  }

  TEST(Rodin_Tuple, MapReturnDifferentType)
  {
    Tuple t{1, 2, 3};
    auto mapped = t.map([](auto x) { return std::to_string(x); });
    EXPECT_EQ(mapped.get<0>(), "1");
    EXPECT_EQ(mapped.get<1>(), "2");
    EXPECT_EQ(mapped.get<2>(), "3");
  }

  //=== Filter ================================================================

  TEST(Rodin_Tuple, FilterIntegrals)
  {
    // Mixed types: ints and doubles.
    Tuple t{1, 2.5, 3, 4.5, 5};
    auto filtered = t.filter<IsIntegral>();
    EXPECT_EQ(filtered.size(), 3);
    EXPECT_EQ(filtered.get<0>(), 1);
    EXPECT_EQ(filtered.get<1>(), 3);
    EXPECT_EQ(filtered.get<2>(), 5);
  }

  TEST(Rodin_Tuple, FilterAlwaysTrue)
  {
    Tuple t{10, 20, 30};
    auto filtered = t.filter<AlwaysTrue>();
    EXPECT_EQ(filtered.size(), t.size());
    EXPECT_EQ(filtered.get<0>(), 10);
    EXPECT_EQ(filtered.get<1>(), 20);
    EXPECT_EQ(filtered.get<2>(), 30);
  }

  TEST(Rodin_Tuple, FilterAlwaysFalse)
  {
    Tuple t{10, 20, 30};
    auto filtered = t.filter<AlwaysFalse>();
    EXPECT_EQ(filtered.size(), 0);
  }

  //=== Apply and IApply ======================================================

  TEST(Rodin_Tuple, ApplyModification)
  {
    Tuple t{10, 20, 30};
    t.apply([](auto& x) { x /= 10; });
    EXPECT_EQ(t.get<0>(), 1);
    EXPECT_EQ(t.get<1>(), 2);
    EXPECT_EQ(t.get<2>(), 3);
  }

  TEST(Rodin_Tuple, ApplyConst)
  {
    const Tuple t{100, 200, 300};
    std::vector<int> collected;
    // The const apply should work without modifying t.
    t.apply([&collected](const auto& x) { collected.push_back(x / 100); });
    EXPECT_EQ(collected.size(), 3);
    EXPECT_EQ(collected[0], 1);
    EXPECT_EQ(collected[1], 2);
    EXPECT_EQ(collected[2], 3);
  }

  TEST(Rodin_Tuple, IApplyWithIndex)
  {
    Tuple t{0, 0, 0};
    t.iapply([](std::size_t i, auto& x) { x = static_cast<int>(i * 5); });
    EXPECT_EQ(t.get<0>(), 0);
    EXPECT_EQ(t.get<1>(), 5);
    EXPECT_EQ(t.get<2>(), 10);
  }

  // //=== Product ===============================================================

  TEST(Rodin_Tuple, ProductDefaultPairing)
  {
    Tuple t1{1, 2};
    Tuple t2{'x', 'y'};
    auto productResult = t1.product(t2);
    EXPECT_EQ(productResult.size(), 4);

    // Each element is a Pair (which is a Tuple of two elements).
    auto p0 = productResult.get<0>(); // Expected: Pair(1, 'x')
    auto p1 = productResult.get<1>(); // Expected: Pair(1, 'y')
    auto p2 = productResult.get<2>(); // Expected: Pair(2, 'x')
    auto p3 = productResult.get<3>(); // Expected: Pair(2, 'y')

    EXPECT_EQ(p0.first(), 1);
    EXPECT_EQ(p0.second(), 'x');
    EXPECT_EQ(p1.first(), 1);
    EXPECT_EQ(p1.second(), 'y');
    EXPECT_EQ(p2.first(), 2);
    EXPECT_EQ(p2.second(), 'x');
    EXPECT_EQ(p3.first(), 2);
    EXPECT_EQ(p3.second(), 'y');
  }

  TEST(Rodin_Tuple, ProductCustomFunction)
  {
    Tuple t1{1, 2};
    Tuple t2{10, 20};

    // Custom function returns the sum of the two arguments.
    auto productResult = t1.product([](auto a, auto b) { return a + b; }, t2);
    EXPECT_EQ(productResult.size(), 4);
    EXPECT_EQ(productResult.get<0>(), 11);
    EXPECT_EQ(productResult.get<1>(), 21);
    EXPECT_EQ(productResult.get<2>(), 12);
    EXPECT_EQ(productResult.get<3>(), 22);
  }

  TEST(Rodin_Tuple, ProductWithEmptyTuple)
  {
    Tuple t{1, 2, 3};
    Tuple empty;
    auto productResult = t.product(empty);
    EXPECT_EQ(productResult.size(), 0);
  }

  TEST(Rodin_Tuple, ProductNonSquare)
  {
    Tuple t1{1, 2, 3};
    Tuple t2{'a', 'b'};
    auto productResult = t1.product(t2);
    EXPECT_EQ(productResult.size(), 6);

    // Check first two pairs
    auto p0 = productResult.get<0>(); // Pair(1, 'a')
    auto p1 = productResult.get<1>(); // Pair(1, 'b')
    EXPECT_EQ(p0.first(), 1);
    EXPECT_EQ(p0.second(), 'a');
    EXPECT_EQ(p1.first(), 1);
    EXPECT_EQ(p1.second(), 'b');
  }

  //=== IndexTuple Helper =====================================================

  TEST(Rodin_Tuple, GenerateIndexTuple)
  {
    // Should create a tuple holding indices: 0, 1, 2, 3, 4.
    auto indexTuple = IndexTuple<0, 5>();
    EXPECT_EQ(indexTuple.size(), 5);
    EXPECT_EQ(indexTuple.get<0>(), 0);
    EXPECT_EQ(indexTuple.get<1>(), 1);
    EXPECT_EQ(indexTuple.get<2>(), 2);
    EXPECT_EQ(indexTuple.get<3>(), 3);
    EXPECT_EQ(indexTuple.get<4>(), 4);
  }

  TEST(Rodin_Tuple, GetConstAndNonConst)
  {
    Tuple t{42, 84};
    t.get<0>() = 100;
    EXPECT_EQ(t.get<0>(), 100);
    const Tuple ct{10, 20};
    EXPECT_EQ(ct.get<0>(), 10);
    EXPECT_EQ(ct.get<1>(), 20);
  }
}

