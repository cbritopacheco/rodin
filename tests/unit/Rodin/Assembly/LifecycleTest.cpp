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
#include "Rodin/Variational/H1.h"
#include "Rodin/Assembly/Default.h"
#include "Rodin/Assembly/Sequential.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Within-field affine ties preserve the physical mass equation and admissible energy.
  template <class Operator>
  void checkSingleFieldAffineIdentification()
  {
    for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      auto mesh = LocalMesh::UniformGrid(geometry, {2, 2});
      mesh.getConnectivity().compute(1, 0);
      mesh.getConnectivity().compute(1, 2);
      mesh.getConnectivity().compute(2, 1);
      auto check = [&](auto& space) {
        TrialFunction u(space);
        TestFunction v(space);
        using System = Math::LinearSystem<Operator, Math::Vector<Real>>;
        Problem<System, decltype(u), decltype(v)> problem(u, v);
        Math::SpatialMatrix<Real> transform(2, 2);
        transform(0, 0) = 0;
        transform(0, 1) = 2;
        transform(1, 0) = 0;
        transform(1, 1) = 1;
        problem = Integral(u, v) - Integral(VectorFunction(5.0, 1.0), v) +
          DirichletBC(u, MatrixFunction(transform) * u, VectorFunction(3.0, 0.0));
        GridFunction exact(space), direction(space);
        exact = VectorFunction(5.0, 1.0);
        direction = VectorFunction(size_t(2), [](const Point& point) {
          const Real value = 1 + point.x() + 2 * point.y();
          return Math::SpatialVector<Real>{2 * value, value};
        });
        for (size_t repeat = 0; repeat < 2; ++repeat)
        {
          problem.assemble();
          const auto& system = problem.getLinearSystem();
          const Math::Matrix<Real> matrix = system.getOperator();
          const Math::Vector<Real> rhs = system.getVector();
          EXPECT_LT((matrix * exact.getData() - rhs).norm(), 1e-11);
          const Math::Vector<Real> recovered = matrix.fullPivLu().solve(rhs);
          EXPECT_LT((recovered - exact.getData()).norm(), 1e-10);
          // Integral of |(2f,f)|^2 on the unit square, f = 1 + x + 2y.
          EXPECT_NEAR(direction.getData().dot(matrix * direction.getData()),
            Real(100) / Real(3), 1e-10);
          for (const auto target : {AssemblyTarget::LHS, AssemblyTarget::RHS})
          {
#ifdef RODIN_USE_OPENMP
            EXPECT_NO_THROW(problem.assemble(target));
#else
            EXPECT_THROW(problem.assemble(target), Alert::Exception);
#endif
            EXPECT_LT((Math::Matrix<Real>(system.getOperator()) - matrix).norm(), 1e-14);
            EXPECT_LT((system.getVector() - rhs).norm(), 1e-14);
          }
        }
      };
      P1<Math::SpatialVector<Real>> linear(mesh, 2);
      check(linear);
      H1<2, Math::SpatialVector<Real>> quadratic(
        std::integral_constant<size_t, 2>{}, mesh, 2);
      check(quadratic);
    }
  }

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

TEST(Rodin_Assembly_SequentialLifecycle, DenseSingleFieldAffineIdentification)
{
  checkSingleFieldAffineIdentification<Math::Matrix<Real>>();
}

TEST(Rodin_Assembly_SequentialLifecycle, SparseSingleFieldAffineIdentification)
{
  checkSingleFieldAffineIdentification<Math::SparseMatrix<Real>>();
}
