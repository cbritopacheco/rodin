/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local-context P1/P2 Poisson rate certification. */

#include <cassert>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <string>

#include <gtest/gtest.h>
#include <petsc.h>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/PETSc.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::PETScPoisson
{
  template <size_t K, class Exact, class Source, class Gradient>
  ErrorNorms solve(const LocalMesh& mesh, const Exact& exact,
    const Source& source, const Gradient& gradient)
  {
    H1 space(std::integral_constant<size_t, K>{}, mesh);
    PETSc::Variational::TrialFunction u(space);
    PETSc::Variational::TestFunction v(space);
    auto stiffness = Integral(Grad(u), Grad(v));
    auto load = Integral(source, v);
    stiffness.setOrder(12);
    load.setOrder(12);
    Problem problem(u, v);
    problem = stiffness - load + DirichletBC(u, exact);
    PETSc::Solver::CG solver(problem);
    solver.setTolerances(1e-12, 1e-14, 1e5, 20000);
    solver.solve();
    EXPECT_TRUE(std::isfinite(solver.getError()));
    return ErrorNorm::compute(mesh, u.getSolution(), exact, gradient, 12);
  }

  class PETScPoissonTest : public ::testing::TestWithParam<Polytope::Type> {};

  template <size_t K>
  void checkRates(Polytope::Type geometry)
  {
    UniformGridHierarchy hierarchy(geometry,
      K == 1 ? std::initializer_list<size_t>{5, 9, 17}
             : std::initializer_list<size_t>{3, 5, 9});
    const size_t dim = hierarchy.getDimension();
    const Real pi = Math::Constants::pi();
    const RealFunction exact([dim, pi](const Point& p)
    {
      Real value = 1;
      for (size_t i = 0; i < dim; ++i)
        value *= std::sin(pi * p(i));
      return value;
    });
    const RealFunction source([dim, pi](const Point& p)
    {
      Real value = Real(dim) * pi * pi;
      for (size_t i = 0; i < dim; ++i)
        value *= std::sin(pi * p(i));
      return value;
    });
    const VectorFunction gradient(dim, [dim, pi](const Point& p)
    {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t i = 0; i < dim; ++i)
      {
        value(i) = pi * std::cos(pi * p(i));
        for (size_t j = 0; j < dim; ++j)
          if (j != i)
            value(i) *= std::sin(pi * p(j));
      }
      return value;
    });
    ErrorHistory history;
    for (size_t level : hierarchy.getLevels())
    {
      const auto mesh = hierarchy.makeMesh(level);
      history.append(hierarchy.getMeshSize(level),
        solve<K>(mesh, exact, source, gradient));
    }
    ASSERT_EQ(history.getSize(), 3);
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      const auto rate = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message() << "L2 " << coarse.getL2()
        << " -> " << fine.getL2() << ", H1 " << coarse.getH1Seminorm()
        << " -> " << fine.getH1Seminorm() << ", rates " << rate.getL2()
        << ", " << rate.getH1Seminorm());
      EXPECT_GT(rate.getL2(), Real(K) + 0.55);
      EXPECT_LT(rate.getL2(), Real(K) + 1.45);
      EXPECT_GT(rate.getH1Seminorm(), Real(K) - 0.3);
      EXPECT_LT(rate.getH1Seminorm(), Real(K) + 0.45);
    }
  }

  TEST_P(PETScPoissonTest, P1OptimalRates) { checkRates<1>(GetParam()); }
  TEST_P(PETScPoissonTest, P2OptimalRates) { checkRates<2>(GetParam()); }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScPoissonTest,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info)
    {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}

int main(int argc, char** argv)
{
  [[maybe_unused]] const PetscErrorCode ierr =
    PetscInitialize(&argc, &argv, nullptr, nullptr);
  assert(ierr == PETSC_SUCCESS);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  [[maybe_unused]] const PetscErrorCode finalizeError = PetscFinalize();
  assert(finalizeError == PETSC_SUCCESS);
  return result;
}
