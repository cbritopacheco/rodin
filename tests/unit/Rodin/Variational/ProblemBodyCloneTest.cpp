/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Partial and untyped formulation clones preserve independent mass systems. */
#include <gtest/gtest.h>
#include <memory>
#include <type_traits>
#include "Rodin/Assembly.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <class Scalar, bool Sparse>
  void checkBodies()
  {
    using Operator =
      std::conditional_t<Sparse, Math::SparseMatrix<Scalar>, Math::Matrix<Scalar>>;
    using Vector = Math::Vector<Scalar>;
    using System = Math::LinearSystem<Operator, Vector>;
    using Function = std::conditional_t<std::is_same_v<Scalar, Complex>,
      ComplexFunction<Complex>, RealFunction<Real>>;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.scale(0.5);
    P1<Scalar> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    const Scalar weight = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, -1);
      else
        return Scalar(2);
    }();
    const Scalar force = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, 2);
      else
        return Scalar(1);
    }();
    const auto mass = Integral(Function(weight) * u, v);
    const auto load = Integral(Function(force), v);
    Math::Matrix<Scalar> expected(3, 3);
    expected << 2, 1, 0, 1, 4, 1, 0, 1, 2;
    expected *= weight / Scalar(12);
    Vector expectedLoad(3);
    expectedLoad << 1, 2, 1;
    expectedLoad *= force / Scalar(4);
    constexpr Real AssemblyTolerance = 1e-12;
    const auto verify = [&](const auto& body) {
      Problem<System, decltype(u), decltype(v)> problem(u, v);
      problem = body;
      problem.assemble();
      const auto& system = problem.getLinearSystem();
      EXPECT_LT((Math::Matrix<Scalar>(system.getOperator()) - expected).norm(),
        AssemblyTolerance);
      EXPECT_LT((system.getVector() - expectedLoad).norm(), AssemblyTolerance);
    };
    auto raw = [&] {
      const auto body = mass - load;
      return std::unique_ptr<std::remove_cvref_t<decltype(body)>>(body.copy());
    }();
    verify(*raw);
    ProblemBodyBase<Scalar> bare(*raw);
    std::unique_ptr<ProblemBodyBase<Scalar>> baseClone(bare.copy());
    bare.getLocalBFIs().clear();
    bare.getLFIs().clear();
    ProblemBody<void, void, Scalar> restored;
    restored.getLocalBFIs().add(baseClone->getLocalBFIs());
    restored.getLFIs().add(baseClone->getLFIs());
    verify(restored);
    {
      using Solution = std::remove_reference_t<decltype(u.getSolution())>;
      BilinearForm<Solution, decltype(space), decltype(space), Operator> assembled(u, v);
      assembled = mass;
      assembled.assemble();
      auto clone = [&] {
        const auto body = assembled - load;
        return std::unique_ptr<std::remove_cvref_t<decltype(body)>>(body.copy());
      }();
      assembled.getOperator().setZero();
      verify(*clone);
    }
    {
      LinearForm<decltype(space), Vector> assembled(v);
      assembled = load;
      assembled.assemble();
      auto clone = [&] {
        const auto body = mass - assembled;
        return std::unique_ptr<std::remove_cvref_t<decltype(body)>>(body.copy());
      }();
      assembled.getVector().setZero();
      verify(*clone);
    }
  }

  /// @brief Real clones preserve raw, base, operator-only and load-only bodies in both storages.
  TEST(ProblemBodyClone, Real)
  {
    checkBodies<Real, false>();
    checkBodies<Real, true>();
  }

  /// @brief Complex clones preserve non-real mass weights and loads independently.
  TEST(ProblemBodyClone, Complex)
  {
    checkBodies<Complex, false>();
    checkBodies<Complex, true>();
  }
}
