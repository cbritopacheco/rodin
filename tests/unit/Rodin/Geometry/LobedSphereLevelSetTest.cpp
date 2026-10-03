/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <LobedSphereLevelSet.h>

#include <array>
#include <cmath>

using namespace Rodin;

namespace Rodin::Tests::Unit
{
  using Vec3 = Math::SpatialVector<Real>;

  TEST(Rodin_Geometry_LobedSphereLevelSet, ZeroLobesIsSphere)
  {
    Examples::LobedSphereLevelSet target;
    target.R0 = Real(0.24);
    target.amp = Real(0.08);
    target.lobes = 0;
    target.phase = Real(0.71);
    for (const Vec3& direction :
      {Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}, Vec3{1, 2, 3}})
    {
      const Vec3 n = direction / direction.norm();
      const Vec3 point = target.c + (target.R0 + target.amp) * n;
      EXPECT_NEAR(target.radius(point), target.R0 + target.amp, 1e-14);
      EXPECT_NEAR(target.phi(point), 0, 1e-14);
      EXPECT_NEAR((target.grad(point) - n).norm(), 0, 1e-14);
    }
  }

  TEST(Rodin_Geometry_LobedSphereLevelSet, BoundedAndAxisBalanced)
  {
    Examples::LobedSphereLevelSet target;
    target.R0 = Real(0.24);
    target.amp = Real(0.08);
    for (int lobes = 1; lobes <= 10; ++lobes)
    {
      target.lobes = lobes;
      const Vec3 px = target.c + Vec3{1, 0, 0};
      const Vec3 py = target.c + Vec3{0, 1, 0};
      const Vec3 pz = target.c + Vec3{0, 0, 1};
      EXPECT_NEAR(target.radius(px), target.radius(py), 1e-14);
      EXPECT_NEAR(target.radius(px), target.radius(pz), 1e-14);
      for (const Vec3& direction :
        {Vec3{1, 1, 1}, Vec3{1, 2, 3}, Vec3{-2, 1, 3}, Vec3{3, -1, -2}})
      {
        const Vec3 point = target.c + direction / direction.norm();
        EXPECT_GE(target.radius(point), target.R0 - target.amp - 1e-14);
        EXPECT_LE(target.radius(point), target.R0 + target.amp + 1e-14);
      }
    }
  }

  TEST(Rodin_Geometry_LobedSphereLevelSet, PhaseRotatesTarget)
  {
    Examples::LobedSphereLevelSet target;
    target.lobes = 5;
    const Vec3 point = target.c + Vec3{Real(0.23), Real(0.12), Real(0.09)};
    const Real baseRadius = target.radius(point);
    const Real phase = Real(0.37);
    const Real cosine = std::cos(phase);
    const Real sine = std::sin(phase);
    const Vec3 x = point - target.c;
    const Vec3 rotated =
      target.c + Vec3{cosine * x(0) - sine * x(1), sine * x(0) + cosine * x(1), x(2)};
    target.phase = phase;
    EXPECT_NEAR(target.radius(rotated), baseRadius, 1e-14);
  }

  TEST(Rodin_Geometry_LobedSphereLevelSet, GradientMatchesFiniteDifferences)
  {
    Examples::LobedSphereLevelSet target;
    target.R0 = Real(0.24);
    target.amp = Real(0.08);
    target.phase = Real(0.37);
    const Real delta = Real(1e-6);
    for (const int lobes : {0, 1, 2, 5, 10})
    {
      target.lobes = lobes;
      for (const Vec3& direction :
        {Vec3{1, 0, 0}, Vec3{0, 0, 1}, Vec3{1, 2, 3}, Vec3{-2, 1, 3}})
      {
        const Vec3 point = target.c + Real(0.26) * direction / direction.norm();
        const Vec3 analytic = target.grad(point);
        for (int i = 0; i < 3; ++i)
        {
          Vec3 plus = point;
          Vec3 minus = point;
          plus(i) += delta;
          minus(i) -= delta;
          const Real numerical =
            (target.phi(plus) - target.phi(minus)) / (Real(2) * delta);
          EXPECT_NEAR(analytic(i), numerical, 2e-8)
            << "lobes=" << lobes << " component=" << i;
        }
      }
    }
  }
}
