/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Numerical lifecycle oracles for sequential matrix and vector assemblers.
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Assembly/Sequential.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Independently owned assemblers recover an exact two-cell mass matrix.
  template <class Scalar, class Operator>
  void checkAssemblerMatrixLifecycle()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    P1<Scalar> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    using Form = BilinearForm<typename decltype(u)::SolutionType, decltype(space),
      decltype(space), Operator>;
    BilinearForm form(u, v);
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -3);
      else
        return Scalar(2);
    }();
    const auto coefficient = [factor] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return ComplexFunction(factor);
      else
        return RealFunction(factor);
    }();
    form = Integral(coefficient * u, v);
    using Assembler = Assembly::Sequential<Operator, Form>;
    const typename Assembler::InputType input{
      space, space, form.getLocalIntegrators(), form.getGlobalIntegrators()};
    Assembler original;
    auto copied = original;
    auto moved = std::move(copied);
    std::unique_ptr<Assembler> cloned(original.copy());
    Math::Matrix<Scalar> expected(3, 3);
    expected << Scalar(1.0 / 3), Scalar(1.0 / 6), Scalar(0), Scalar(1.0 / 6),
      Scalar(2.0 / 3), Scalar(1.0 / 6), Scalar(0), Scalar(1.0 / 6), Scalar(1.0 / 3);
    expected *= factor;
    auto check = [&](auto& assembler) {
      Operator result;
      for (size_t stage = 0; stage < 2; ++stage)
      {
        assembler.execute(result, input);
        Math::Matrix<Scalar> actual;
        if constexpr (requires { result.rows(); })
          actual = result;
        else
        {
          Math::SparseMatrix<Scalar> sparse(3, 3);
          sparse.setFromTriplets(result.begin(), result.end());
          actual = sparse;
        }
        EXPECT_LT((actual - expected).norm(), 1e-12);
      }
    };
    check(original);
    check(moved);
    check(*cloned);
  }

  /// @brief Vector assembler copies recover exact nodal integrals without accumulation.
  template <class Scalar>
  void checkAssemblerVectorLifecycle()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    P1<Scalar> space(mesh);
    TestFunction v(space);
    LinearForm form(v);
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -3);
      else
        return Scalar(2);
    }();
    const auto coefficient = [factor] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return ComplexFunction(factor);
      else
        return RealFunction(factor);
    }();
    form = Integral(coefficient, v);
    using Assembler = Assembly::Sequential<Math::Vector<Scalar>, decltype(form)>;
    const typename Assembler::InputType input{space, form.getIntegrators()};
    Assembler original;
    auto copied = original;
    auto moved = std::move(copied);
    std::unique_ptr<Assembler> cloned(original.copy());
    Math::Vector<Scalar> expected(3);
    expected << factor * Scalar(0.5), factor, factor * Scalar(0.5);
    auto check = [&](auto& assembler) {
      Math::Vector<Scalar> result;
      for (size_t stage = 0; stage < 2; ++stage)
      {
        assembler.execute(result, input);
        EXPECT_LT((result - expected).norm(), 1e-12);
      }
    };
    check(original);
    check(moved);
    check(*cloned);
  }
}

TEST(Rodin_Assembly_SequentialLifecycle, RealMatricesAndVector)
{
  checkAssemblerMatrixLifecycle<Real, Math::Matrix<Real>>();
  checkAssemblerMatrixLifecycle<Real, Math::SparseMatrix<Real>>();
  checkAssemblerMatrixLifecycle<Real, std::vector<Math::SparseTriplet<Real>>>();
  checkAssemblerVectorLifecycle<Real>();
}

TEST(Rodin_Assembly_SequentialLifecycle, ComplexMatricesAndVector)
{
  checkAssemblerMatrixLifecycle<Complex, Math::Matrix<Complex>>();
  checkAssemblerMatrixLifecycle<Complex, Math::SparseMatrix<Complex>>();
  checkAssemblerMatrixLifecycle<Complex, std::vector<Math::SparseTriplet<Complex>>>();
  checkAssemblerVectorLifecycle<Complex>();
}
