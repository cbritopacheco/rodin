/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly/Input.h"
#include "Rodin/FormLanguage/List.h"
#include "Rodin/Types.h"

using namespace Rodin;
using namespace Rodin::Assembly;

// Mock classes for testing
class MockFES
{
  public:
    using ScalarType = Real;
};

class MockTrialFES
{
  public:
    using ScalarType = Real;
};

class MockTestFES
{
  public:
    using ScalarType = Real;
};

class MockLinearFormIntegratorBase
{
  public:
    MockLinearFormIntegratorBase() = default;
};

class MockLocalBilinearFormIntegratorBase
{
  public:
    MockLocalBilinearFormIntegratorBase() = default;
};

class MockGlobalBilinearFormIntegratorBase
{
  public:
    MockGlobalBilinearFormIntegratorBase() = default;
};

// Specialize Traits for our mock classes
namespace Rodin::FormLanguage
{
  template <>
  struct Traits<MockFES>
  {
    using ScalarType = Real;
  };

  template <>
  struct Traits<MockTrialFES>
  {
    using ScalarType = Real;
  };

  template <>
  struct Traits<MockTestFES>
  {
    using ScalarType = Real;
  };
}

class AssemblyInputTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test LinearFormAssemblyInput basic functionality
TEST_F(AssemblyInputTest, LinearFormAssemblyInputBasic)
{
  // Create mock objects
  MockFES fes;
  FormLanguage::List<MockLinearFormIntegratorBase> lfis;
  
  // Create LinearFormAssemblyInput
  LinearFormAssemblyInput<MockFES> input(fes, lfis);
  
  // Test getFES method
  const MockFES& retrievedFES = input.getFES();
  EXPECT_EQ(&retrievedFES, &fes);  // Should be the same object
  
  // Test getLFIs method
  auto& retrievedLFIs = input.getLFIs();
  EXPECT_EQ(&retrievedLFIs, &lfis);  // Should be the same object
}

// Test LinearFormAssemblyInput type traits
TEST_F(AssemblyInputTest, LinearFormAssemblyInputTypes)
{
  using InputType = LinearFormAssemblyInput<MockFES>;
  
  // Test type aliases
  static_assert(std::is_same_v<typename InputType::FESType, MockFES>);
  static_assert(std::is_same_v<typename InputType::ScalarType, Real>);
  
  // Test that LinearFormIntegratorBaseList is properly defined
  using ListType = typename InputType::LinearFormIntegratorBaseList;
  static_assert(std::is_base_of_v<FormLanguage::Base, ListType>);
}

// Test LinearFormAssemblyInput const correctness
TEST_F(AssemblyInputTest, LinearFormAssemblyInputConstCorrectness)
{
  MockFES fes;
  FormLanguage::List<MockLinearFormIntegratorBase> lfis;
  
  LinearFormAssemblyInput<MockFES> input(fes, lfis);
  const LinearFormAssemblyInput<MockFES>& constInput = input;
  
  // Test that const methods work
  const MockFES& retrievedFES = constInput.getFES();
  EXPECT_EQ(&retrievedFES, &fes);
  
  auto& retrievedLFIs = constInput.getLFIs();
  EXPECT_EQ(&retrievedLFIs, &lfis);
}

// Test BilinearFormAssemblyInput basic functionality
TEST_F(AssemblyInputTest, BilinearFormAssemblyInputBasic)
{
  // Create mock objects
  MockTrialFES trialFES;
  MockTestFES testFES;
  FormLanguage::List<MockLocalBilinearFormIntegratorBase> lbfis;
  FormLanguage::List<MockGlobalBilinearFormIntegratorBase> gbfis;
  
  // Create BilinearFormAssemblyInput
  BilinearFormAssemblyInput<MockTrialFES, MockTestFES> input(
    trialFES, testFES, lbfis, gbfis);
  
  // Test getter methods
  const MockTrialFES& retrievedTrialFES = input.getTrialFES();
  EXPECT_EQ(&retrievedTrialFES, &trialFES);
  
  const MockTestFES& retrievedTestFES = input.getTestFES();
  EXPECT_EQ(&retrievedTestFES, &testFES);
  
  auto& retrievedLBFIs = input.getLocalBFIs();
  EXPECT_EQ(&retrievedLBFIs, &lbfis);
  
  auto& retrievedGBFIs = input.getGlobalBFIs();
  EXPECT_EQ(&retrievedGBFIs, &gbfis);
}

// Test BilinearFormAssemblyInput type traits
TEST_F(AssemblyInputTest, BilinearFormAssemblyInputTypes)
{
  using InputType = BilinearFormAssemblyInput<MockTrialFES, MockTestFES>;
  
  // Test type aliases
  static_assert(std::is_same_v<typename InputType::TrialFESType, MockTrialFES>);
  static_assert(std::is_same_v<typename InputType::TestFESType, MockTestFES>);
  static_assert(std::is_same_v<typename InputType::TrialFESScalarType, Real>);
  static_assert(std::is_same_v<typename InputType::TestFESScalarType, Real>);
  static_assert(std::is_same_v<typename InputType::ScalarType, Real>);
  
  // Test that the list types are properly defined
  using LocalListType = typename InputType::LocalBilinearFormIntegratorBaseListType;
  using GlobalListType = typename InputType::GlobalBilinearFormIntegratorBaseListType;
  
  static_assert(std::is_base_of_v<FormLanguage::Base, LocalListType>);
  static_assert(std::is_base_of_v<FormLanguage::Base, GlobalListType>);
}

// Test BilinearFormAssemblyInput const correctness
TEST_F(AssemblyInputTest, BilinearFormAssemblyInputConstCorrectness)
{
  MockTrialFES trialFES;
  MockTestFES testFES;
  FormLanguage::List<MockLocalBilinearFormIntegratorBase> lbfis;
  FormLanguage::List<MockGlobalBilinearFormIntegratorBase> gbfis;
  
  BilinearFormAssemblyInput<MockTrialFES, MockTestFES> input(
    trialFES, testFES, lbfis, gbfis);
  const BilinearFormAssemblyInput<MockTrialFES, MockTestFES>& constInput = input;
  
  // Test that const methods work
  const MockTrialFES& retrievedTrialFES = constInput.getTrialFES();
  EXPECT_EQ(&retrievedTrialFES, &trialFES);
  
  const MockTestFES& retrievedTestFES = constInput.getTestFES();
  EXPECT_EQ(&retrievedTestFES, &testFES);
  
  auto& retrievedLBFIs = constInput.getLocalBFIs();
  EXPECT_EQ(&retrievedLBFIs, &lbfis);
  
  auto& retrievedGBFIs = constInput.getGlobalBFIs();
  EXPECT_EQ(&retrievedGBFIs, &gbfis);
}

// Test copy semantics
TEST_F(AssemblyInputTest, CopySemantics)
{
  MockFES fes;
  FormLanguage::List<MockLinearFormIntegratorBase> lfis;
  
  LinearFormAssemblyInput<MockFES> input1(fes, lfis);
  
  // Test copy constructor
  LinearFormAssemblyInput<MockFES> input2 = input1;
  
  // Both should reference the same underlying objects
  EXPECT_EQ(&input1.getFES(), &input2.getFES());
  EXPECT_EQ(&input1.getLFIs(), &input2.getLFIs());
  
  // Test copy assignment
  MockFES fes2;
  FormLanguage::List<MockLinearFormIntegratorBase> lfis2;
  LinearFormAssemblyInput<MockFES> input3(fes2, lfis2);
  
  input3 = input1;
  EXPECT_EQ(&input3.getFES(), &input1.getFES());
  EXPECT_EQ(&input3.getLFIs(), &input1.getLFIs());
}

// Test BilinearFormTupleAssemblyInput basic functionality
TEST_F(AssemblyInputTest, BilinearFormTupleAssemblyInputBasic)
{
  using TupleType = Tuple<int, double, std::string>;
  using OffsetsType = std::array<Pair<size_t, size_t>, 3>;
  
  OffsetsType offsets = {{{0, 0}, {1, 1}, {2, 2}}};
  TupleType tuple = makeTuple(42, 3.14, std::string("test"));
  
  BilinearFormTupleAssemblyInput<int, double, std::string> input(
    10, 20, offsets, tuple);
  
  // Test basic getters
  EXPECT_EQ(input.getRows(), 10);
  EXPECT_EQ(input.getCols(), 20);
  
  const auto& retrievedOffsets = input.getOffsets();
  EXPECT_EQ(retrievedOffsets.size(), 3);
  EXPECT_EQ(retrievedOffsets[0].first, 0);
  EXPECT_EQ(retrievedOffsets[0].second, 0);
  EXPECT_EQ(retrievedOffsets[1].first, 1);
  EXPECT_EQ(retrievedOffsets[1].second, 1);
  
  const auto& retrievedTuple = input.getInputs();
  EXPECT_EQ(get<0>(retrievedTuple), 42);
  EXPECT_DOUBLE_EQ(get<1>(retrievedTuple), 3.14);
  EXPECT_EQ(get<2>(retrievedTuple), "test");
}

// Test BilinearFormTupleAssemblyInput type traits
TEST_F(AssemblyInputTest, BilinearFormTupleAssemblyInputTypes)
{
  using InputType = BilinearFormTupleAssemblyInput<int, double, std::string>;
  
  // Test Size constant
  static_assert(InputType::Size == 3);
  
  // Test Offsets type
  using OffsetsType = typename InputType::Offsets;
  static_assert(std::is_same_v<OffsetsType, std::array<Pair<size_t, size_t>, 3>>);
}

// Test empty tuple case
TEST_F(AssemblyInputTest, BilinearFormTupleAssemblyInputEmpty)
{
  using EmptyInputType = BilinearFormTupleAssemblyInput<>;
  
  // Test Size for empty tuple
  static_assert(EmptyInputType::Size == 0);
  
  using EmptyOffsetsType = typename EmptyInputType::Offsets;
  static_assert(std::is_same_v<EmptyOffsetsType, std::array<Pair<size_t, size_t>, 0>>);
  
  EmptyOffsetsType emptyOffsets = {};
  Tuple<> emptyTuple = makeTuple();
  
  EmptyInputType emptyInput(0, 0, emptyOffsets, emptyTuple);
  EXPECT_EQ(emptyInput.getRows(), 0);
  EXPECT_EQ(emptyInput.getCols(), 0);
}

// Test different tuple sizes
TEST_F(AssemblyInputTest, BilinearFormTupleAssemblyInputDifferentSizes)
{
  // Test single element tuple
  using SingleInputType = BilinearFormTupleAssemblyInput<int>;
  static_assert(SingleInputType::Size == 1);
  
  // Test large tuple
  using LargeInputType = BilinearFormTupleAssemblyInput<int, double, float, char, bool>;
  static_assert(LargeInputType::Size == 5);
}