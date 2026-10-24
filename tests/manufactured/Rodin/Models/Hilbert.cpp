/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Manufactured scalar and vector Hilbert extension regressions. */
#include <gtest/gtest.h>
#include <cmath>
#include "Rodin/Variational.h"
#include "Rodin/Hilbert/H1a.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Manufactured
{
  /// @brief Exact affine Dirichlet data survive successive applications and parameter changes.
  TEST(Rodin_Manufactured_Hilbert, AffineDirichlet)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {5, 5});
    mesh.scale(Real(0.25));
    mesh.getConnectivity().compute(1, 2);
    P1<Real> fes(mesh);
    TrialFunction trial(fes);
    using Solution = std::remove_cvref_t<decltype(trial.getSolution())>;
    Hilbert::H1a<Solution, decltype(fes)> extension(fes);
    EXPECT_EQ(&extension.getFiniteElementSpace(), &fes);
    EXPECT_EQ(&extension.getProblem(), &std::as_const(extension).getProblem());
    const auto exact =
      RealFunction([](const Point& p) { return 1 + 2 * p(0) - 3 * p(1); });
    const auto load = Integral(exact, extension.getTestFunction());
    using Extension = decltype(extension);
    static_assert(std::is_invocable_v<Extension&, decltype(load)>);
    static_assert(!std::is_invocable_v<const Extension&, decltype(load)>);
    extension += DirichletBC(extension.getTrialFunction(), exact);
    for (const Real alpha : {Real(0.2), Real(0.7)})
    {
      extension.setAlpha(alpha);
      const auto field = extension(Integral(exact, extension.getTestFunction()));
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
        EXPECT_NEAR(field(point), exact(point), 1e-9);
      }
      EXPECT_GT(
        std::as_const(extension).getProblem().getLinearSystem().getOperator().rows(), 0);
    }
  }

  /// @brief A cosine with zero normal flux verifies the scalar regularization coefficient.
  TEST(Rodin_Manufactured_Hilbert, ScalarNeumannCosine)
  {
    constexpr Real Alpha = Real(0.3);
    constexpr Real CoarseErrorTolerance = Real(0.03);
    // P1 second-order error should shrink by one quarter when halving the mesh size.
    constexpr Real MaxRefinementRatio = Real(0.35);
    const Real pi = std::acos(Real(-1));
    Real previous = 0;
    for (const size_t nodes : {size_t(9), size_t(17)})
    {
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
      mesh.scale(Real(1) / (nodes - 1));
      P1<Real> fes(mesh);
      TrialFunction trial(fes);
      using Solution = std::remove_cvref_t<decltype(trial.getSolution())>;
      Hilbert::H1a<Solution, decltype(fes)> extension(fes);
      extension.setAlpha(Alpha);
      const auto load = RealFunction([pi](const Point& p) {
        return (1 + 2 * Alpha * Alpha * pi * pi) * std::cos(pi * p(0)) *
          std::cos(pi * p(1));
      });
      const auto field = extension(Integral(load, extension.getTestFunction()));
      Real error = 0;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
        const Real delta =
          field(point) - std::cos(pi * point(0)) * std::cos(pi * point(1));
        error += cell->getMeasure() * delta * delta;
      }
      error = std::sqrt(error);
      EXPECT_LT(error, CoarseErrorTolerance);
      if (previous > 0)
      {
        EXPECT_LT(error, MaxRefinementRatio * previous);
      }
      previous = error;
    }
  }

  /// @brief Two independent zero-flux cosine components verify the vector Jacobian branch.
  TEST(Rodin_Manufactured_Hilbert, VectorNeumannCosine)
  {
    constexpr Real Alpha = Real(0.4);
    constexpr Real ErrorTolerance = Real(0.015);
    const Real pi = std::acos(Real(-1));
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {17, 17});
    mesh.scale(Real(1) / 16);
    P1<Math::SpatialVector<Real>> fes(mesh, 2);
    TrialFunction trial(fes);
    using Solution = std::remove_cvref_t<decltype(trial.getSolution())>;
    Hilbert::H1a<Solution, decltype(fes)> extension(fes);
    extension.setAlpha(Alpha);
    const auto load = VectorFunction(
      [pi](
        const Point& p) { return (1 + Alpha * Alpha * pi * pi) * std::cos(pi * p(0)); },
      [pi](const Point& p) {
        return 2 * (1 + Alpha * Alpha * pi * pi) * std::cos(pi * p(1));
      });
    const auto field = extension(Integral(load, extension.getTestFunction()));
    Real error = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      const auto value = field(point);
      const Real dx = value(0) - std::cos(pi * point(0));
      const Real dy = value(1) - 2 * std::cos(pi * point(1));
      error += cell->getMeasure() * (dx * dx + dy * dy);
    }
    EXPECT_LT(std::sqrt(error), ErrorTolerance);
  }
}
