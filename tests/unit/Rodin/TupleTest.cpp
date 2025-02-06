#include <gtest/gtest.h>

#include <Rodin/Tuple.h>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Tuple, SanityTest_VariantConstructor)
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

  TEST(Rodin_Tuple, SanityTest_CopyConstructor)
  {
    Tuple t1{1, 'a', 1.0};
    Tuple t2(t1);

    EXPECT_EQ(t2.size(), 3);
    EXPECT_EQ(t2.get<0>(), 1);
    EXPECT_EQ(t2.get<1>(), 'a');
    EXPECT_EQ(t2.get<2>(), 1.0);
  }

  TEST(Rodin_Tuple, SanityTest_CopyAssignment)
  {
    Tuple t1{1, 2, 3};

    EXPECT_EQ(t1.size(), 3);
    EXPECT_EQ(t1.get<0>(), 1);
    EXPECT_EQ(t1.get<1>(), 2);
    EXPECT_EQ(t1.get<2>(), 3);

    Tuple t2{4, 5, 6};

    EXPECT_EQ(t2.size(), 3);
    EXPECT_EQ(t2.get<0>(), 4);
    EXPECT_EQ(t2.get<1>(), 5);
    EXPECT_EQ(t2.get<2>(), 6);

    t2 = t1;

    EXPECT_EQ(t2.size(), 3);
    EXPECT_EQ(t2.get<0>(), 1);
    EXPECT_EQ(t2.get<1>(), 2);
    EXPECT_EQ(t2.get<2>(), 3);
  }

  TEST(Rodin_Tuple, SanityTest_Concatenate)
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

  TEST(Rodin_Tuple, SanityTest_Zip)
  {
    Tuple t1{1, 'a', 1.0};
    Tuple t2{2, 'b', 2.0};

    auto r = t1.zip(t2);

    EXPECT_EQ(r.size(), t1.size());
    EXPECT_EQ(r.size(), t2.size());
  }

  TEST(Rodin_Tuple, SanityTest_Reduce_Sum)
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

  TEST(Rodin_Tuple, SanityTest_Reduce_Product)
  {
    {
      Tuple t{2, 3, 4};
      auto result = t.reduce([](auto a, auto b) { return a * b; });
      EXPECT_EQ(result, 24);
    }
  }
}

