/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Semilinear Poisson Newton and h-rate validation. */

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <string>

#include <gtest/gtest.h>

#include "Convergence.h"
#include "../../NonlinearPoisson.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::NonlinearPoisson
{

  template <size_t K>
  void checkRates(Polytope::Type geometry)
  {
    UniformGridHierarchy hierarchy(geometry,
      K == 1 ? std::initializer_list<size_t>{5, 9, 17}
             : std::initializer_list<size_t>{3, 5, 9});
    ErrorHistory history;
    for (size_t level : hierarchy.getLevels())
    {
      const auto mesh = hierarchy.makeMesh(level);
      history.append(
        hierarchy.getMeshSize(level), NonlinearPoissonProblem(mesh, 12).solve<K>());
    }
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      const auto rate = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message()
        << "L2 " << coarse.getL2() << " -> " << fine.getL2() << ", H1 "
        << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm() << ", rates "
        << rate.getL2() << ", " << rate.getH1Seminorm());
      EXPECT_GT(rate.getL2(), Real(K) + 0.5);
      EXPECT_LT(rate.getL2(), Real(K) + 1.5);
      EXPECT_GT(rate.getH1Seminorm(), Real(K) - 0.3);
      EXPECT_LT(rate.getH1Seminorm(), Real(K) + 0.5);
    }
  }

  class NonlinearPoissonTest : public ::testing::TestWithParam<Polytope::Type>
  {};
  TEST_P(NonlinearPoissonTest, P1OptimalRates)
  {
    checkRates<1>(GetParam());
  }
  TEST_P(NonlinearPoissonTest, P2OptimalRates)
  {
    checkRates<2>(GetParam());
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, NonlinearPoissonTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
