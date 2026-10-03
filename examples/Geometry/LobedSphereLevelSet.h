/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXAMPLES_GEOMETRY_LOBEDSPHERELEVELSET_H
#define RODIN_EXAMPLES_GEOMETRY_LOBEDSPHERELEVELSET_H

#include <Rodin/Math.h>

#include <cmath>

namespace Rodin::Examples
{
  struct LobedSphereLevelSet
  {
      using Vec3 = Math::SpatialVector<Real>;

      Vec3 c = {Real(0.5), Real(0.5), Real(0.5)};
      Real R0 = Real(0.25);
      Real amp = Real(0.05);
      Real lobes = Real(6);
      Real phase = Real(0);

      Real radius(const Vec3& p) const
      {
        const Vec3 x = p - c;
        const Real r = x.norm();
        if (r <= Real(1e-14))
          return R0 + amp;
        const Vec3 direction = rotateZ(x / r, -phase);
        const Real frequency = lobes;
        return R0 +
          amp / Real(3) *
          (std::cos(frequency * direction(0)) + std::cos(frequency * direction(1)) +
            std::cos(frequency * direction(2)));
      }

      Real phi(const Vec3& p) const
      {
        return (p - c).norm() - radius(p);
      }

      Vec3 grad(const Vec3& p) const
      {
        const Vec3 x = p - c;
        const Real r = x.norm();
        if (r <= Real(1e-14))
          return Vec3{0, 0, 0};
        const Vec3 n = x / r;
        const Vec3 direction = rotateZ(n, -phase);
        const Real frequency = lobes;
        Vec3 angularGradient(3);
        for (int i = 0; i < 3; ++i)
          angularGradient(i) =
            -amp * frequency / Real(3) * std::sin(frequency * direction(i));
        angularGradient = rotateZ(angularGradient, phase);
        return n - (angularGradient - n * n.dot(angularGradient)) / r;
      }

    private:
      static Vec3 rotateZ(const Vec3& x, Real angle)
      {
        const Real cosine = std::cos(angle);
        const Real sine = std::sin(angle);
        Vec3 out(3);
        out(0) = cosine * x(0) - sine * x(1);
        out(1) = sine * x(0) + cosine * x(1);
        out(2) = x(2);
        return out;
      }
  };
}

#endif
