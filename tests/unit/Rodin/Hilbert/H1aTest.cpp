/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Hilbert/H1a.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/// @brief Copies and moves solve into their own trial functions.
TEST(H1a, IndependentWorkingProblem)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  P1<Real> space(mesh);
  using Extension =
    Hilbert::H1a<GridFunction<decltype(space), Math::Vector<Real>>, decltype(space)>;
  Extension source(space);
  source.setAlpha(0.2);
  const auto original = source(Integral(RealFunction(2.), source.getTestFunction()));
  EXPECT_LT((original.getData().array() - 2).abs().maxCoeff(), 1e-10);
  Extension copy(source);
  const auto copied = copy(Integral(RealFunction(3.), copy.getTestFunction()));
  EXPECT_LT((copied.getData().array() - 3).abs().maxCoeff(), 1e-10);
  EXPECT_LT(
    (source.getTrialFunction().getSolution().getData().array() - 2).abs().maxCoeff(),
    1e-10);
  auto moved = [&] {
    Extension temporary(copy);
    return Extension(std::move(temporary));
  }();
  const auto result = moved(Integral(RealFunction(4.), moved.getTestFunction()));
  EXPECT_LT((result.getData().array() - 4).abs().maxCoeff(), 1e-10);
  EXPECT_LT(
    (copy.getTrialFunction().getSolution().getData().array() - 3).abs().maxCoeff(),
    1e-10);
}

/// @brief Boundary conditions survive owner changes without borrowing the old problem.
TEST(H1a, CopyRetainsBoundaryConditions)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  mesh.getConnectivity().compute(1, 2);
  P1<Real> space(mesh);
  using Extension =
    Hilbert::H1a<GridFunction<decltype(space), Math::Vector<Real>>, decltype(space)>;
  auto extension = [&] {
    Extension source(space);
    source += DirichletBC(source.getTrialFunction(), RealFunction(2.));
    return Extension(source);
  }();
  const auto solution =
    extension(Integral(RealFunction(2.), extension.getTestFunction()));
  EXPECT_LT((solution.getData().array() - 2).abs().maxCoeff(), 1e-10);
}
