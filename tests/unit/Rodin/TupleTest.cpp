#include <gtest/gtest.h>

#include <Rodin/Tuple.h>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Tuple, SanityTest_Constructor)
  {
    Tuple t{1, 'a', 1.0};

    EXPECT_EQ(t.size(), 3);
    EXPECT_EQ(t.get<0>(), 1);
    EXPECT_EQ(t.get<1>(), 'a');
    EXPECT_EQ(t.get<2>(), 1.0);
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
}

