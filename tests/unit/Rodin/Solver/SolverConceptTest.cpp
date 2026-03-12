/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

using namespace Rodin;

// Simple concrete test class that would inherit from LinearSolverBase
// We can't test LinearSolverBase directly since it's abstract, but we can test
// the concepts and type system it uses
template <class LinearSystem>
class TestSolver
{
  public:
    using ScalarType = typename FormLanguage::Traits<LinearSystem>::ScalarType;
    using VectorType = typename FormLanguage::Traits<LinearSystem>::VectorType;
    using OperatorType = typename FormLanguage::Traits<LinearSystem>::OperatorType;
    
    TestSolver() = default;
    
    // Mock solve method
    void solve() { m_solved = true; }
    bool isSolved() const { return m_solved; }
    
  private:
    bool m_solved = false;
};

class SolverConceptTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Test that the LinearSystem concept works with our types
TEST_F(SolverConceptTest, LinearSystemTypes)
{
  using TestLinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using TestSolverType = TestSolver<TestLinearSystem>;
  
  // Test type extraction
  static_assert(std::is_same_v<typename TestSolverType::ScalarType, Real>);
  static_assert(std::is_same_v<typename TestSolverType::VectorType, Math::Vector<Real>>);
  static_assert(std::is_same_v<typename TestSolverType::OperatorType, Math::SparseMatrix<Real>>);
}

// Test basic solver construction and interface
TEST_F(SolverConceptTest, BasicSolverInterface)
{
  using TestLinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using TestSolverType = TestSolver<TestLinearSystem>;
  
  // Test construction
  TestSolverType solver;
  EXPECT_FALSE(solver.isSolved());
  
  // Test solve operation
  solver.solve();
  EXPECT_TRUE(solver.isSolved());
}

// Test solver with different scalar types
TEST_F(SolverConceptTest, DifferentScalarTypes)
{
  // Test with Real
  using RealLinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using RealSolver = TestSolver<RealLinearSystem>;
  
  static_assert(std::is_same_v<typename RealSolver::ScalarType, Real>);
  
  // Test with Complex
  using ComplexLinearSystem = Math::LinearSystem<Math::SparseMatrix<Complex>, Math::Vector<Complex>>;
  using ComplexSolver = TestSolver<ComplexLinearSystem>;
  
  static_assert(std::is_same_v<typename ComplexSolver::ScalarType, Complex>);
  
  // Test construction of both types
  RealSolver realSolver;
  ComplexSolver complexSolver;
  
  EXPECT_FALSE(realSolver.isSolved());
  EXPECT_FALSE(complexSolver.isSolved());
}

// Test solver copy semantics
TEST_F(SolverConceptTest, CopySemantics)
{
  using TestLinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using TestSolverType = TestSolver<TestLinearSystem>;
  
  TestSolverType solver1;
  solver1.solve();
  EXPECT_TRUE(solver1.isSolved());
  
  // Test copy constructor
  TestSolverType solver2 = solver1;
  EXPECT_TRUE(solver2.isSolved());
  
  // Test copy assignment
  TestSolverType solver3;
  EXPECT_FALSE(solver3.isSolved());
  solver3 = solver1;
  EXPECT_TRUE(solver3.isSolved());
}

// Test solver move semantics
TEST_F(SolverConceptTest, MoveSemantics)
{
  using TestLinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using TestSolverType = TestSolver<TestLinearSystem>;
  
  TestSolverType solver1;
  solver1.solve();
  EXPECT_TRUE(solver1.isSolved());
  
  // Test move constructor
  TestSolverType solver2 = std::move(solver1);
  EXPECT_TRUE(solver2.isSolved());
  
  // Test move assignment
  TestSolverType solver3;
  TestSolverType solver4;
  solver4.solve();
  solver3 = std::move(solver4);
  EXPECT_TRUE(solver3.isSolved());
}

// Test multiple solver instances independence
TEST_F(SolverConceptTest, MultipleInstances)
{
  using TestLinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using TestSolverType = TestSolver<TestLinearSystem>;
  
  TestSolverType solver1;
  TestSolverType solver2;
  TestSolverType solver3;
  
  // Initially all should be unsolved
  EXPECT_FALSE(solver1.isSolved());
  EXPECT_FALSE(solver2.isSolved());
  EXPECT_FALSE(solver3.isSolved());
  
  // Solve only one
  solver2.solve();
  
  // Check independence
  EXPECT_FALSE(solver1.isSolved());
  EXPECT_TRUE(solver2.isSolved());
  EXPECT_FALSE(solver3.isSolved());
}

// Test solver with different matrix and vector combinations
TEST_F(SolverConceptTest, DifferentMatrixVectorCombinations)
{
  // Real sparse matrix with real vector
  using RealSparseSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using RealSparseSolver = TestSolver<RealSparseSystem>;
  
  // Complex sparse matrix with complex vector  
  using ComplexSparseSystem = Math::LinearSystem<Math::SparseMatrix<Complex>, Math::Vector<Complex>>;
  using ComplexSparseSolver = TestSolver<ComplexSparseSystem>;
  
  // Test that both can be constructed and used
  RealSparseSolver realSolver;
  ComplexSparseSolver complexSolver;
  
  realSolver.solve();
  complexSolver.solve();
  
  EXPECT_TRUE(realSolver.isSolved());
  EXPECT_TRUE(complexSolver.isSolved());
  
  // Verify type consistency
  static_assert(std::is_same_v<typename RealSparseSolver::ScalarType, Real>);
  static_assert(std::is_same_v<typename ComplexSparseSolver::ScalarType, Complex>);
}

// Test solver type traits and metaprogramming
TEST_F(SolverConceptTest, TypeTraitsAndMetaprogramming)
{
  using TestLinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using TestSolverType = TestSolver<TestLinearSystem>;
  
  // Test basic type properties
  static_assert(std::is_default_constructible_v<TestSolverType>);
  static_assert(std::is_copy_constructible_v<TestSolverType>);
  static_assert(std::is_copy_assignable_v<TestSolverType>);
  static_assert(std::is_move_constructible_v<TestSolverType>);
  static_assert(std::is_move_assignable_v<TestSolverType>);
  
  // Test that ScalarType is arithmetic
  static_assert(std::is_arithmetic_v<typename TestSolverType::ScalarType>);
  
  // Test vector type properties
  using VectorType = typename TestSolverType::VectorType;
  static_assert(std::is_default_constructible_v<VectorType>);
  
  // Test operator type properties  
  using OperatorType = typename TestSolverType::OperatorType;
  static_assert(std::is_default_constructible_v<OperatorType>);
}

// Test const correctness
TEST_F(SolverConceptTest, ConstCorrectness)
{
  using TestLinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using TestSolverType = TestSolver<TestLinearSystem>;
  
  TestSolverType solver;
  solver.solve();
  
  // Test const methods
  const TestSolverType& constSolver = solver;
  EXPECT_TRUE(constSolver.isSolved());
  
  // const methods should not modify state
  bool solved1 = constSolver.isSolved();
  bool solved2 = constSolver.isSolved();
  EXPECT_EQ(solved1, solved2);
}

// Test solver lifecycle and RAII
TEST_F(SolverConceptTest, LifecycleAndRAII)
{
  using TestLinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  using TestSolverType = TestSolver<TestLinearSystem>;
  
  // Test scope-based lifetime
  {
    TestSolverType localSolver;
    localSolver.solve();
    EXPECT_TRUE(localSolver.isSolved());
  }  // Should destroy properly
  
  // Test dynamic allocation
  auto* dynamicSolver = new TestSolverType();
  dynamicSolver->solve();
  EXPECT_TRUE(dynamicSolver->isSolved());
  delete dynamicSolver;
  
  // Test smart pointers
  auto smartSolver = std::make_unique<TestSolverType>();
  smartSolver->solve();
  EXPECT_TRUE(smartSolver->isSolved());
  // Should automatically clean up
}