/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/Flow.h"
#include "Rodin/Advection/Lagrangian.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Advection;

namespace Rodin::Tests::Manufactured::AdvectionLagrangian3D
{
  template <size_t NX, size_t NY, size_t NZ>
  class ManufacturedAdvection3DTest : public ::testing::TestWithParam<Polytope::Type>
  {
  protected:
    Mesh<Context::Local> getMesh()
    {
      const auto geom = GetParam();
      Mesh<Context::Local> mesh = Mesh<Context::Local>::UniformGrid(geom, { NX, NY, NZ });
      mesh.scale(1.0 / (NX - 1));
      mesh.getConnectivity().compute(2, 3);
      mesh.getConnectivity().compute(3, 2);
      return mesh;
    }

    template <class GF>
    static void checkL2CentroidError(
        const Mesh<Context::Local>& mesh,
        const GF& uh,
        const Math::SpatialVector<Real>& vel,
        Real dt,
        Real atol,
        Real rtol)
    {
      const Real pi = Math::Constants::pi();
      auto inside01 = [](Real z) { return z >= 0.0 && z <= 1.0; };

      const size_t cd = mesh.getDimension();
      ASSERT_EQ(cd, 3u);

      Real err2 = 0.0;
      Real ref2 = 0.0;
      size_t n = 0;

      Math::SpatialVector<Real> rc{{0.5, 0.5, 0.5}};

      for (Index c = 0; c < mesh.getPolytopeCount(cd); ++c)
      {
        auto itc = mesh.getPolytope(cd, c);
        const auto& cell = *itc;

        Geometry::Point p(cell, rc);
        const auto pc = p.getPhysicalCoordinates();
        const Real x = pc[0], y = pc[1], z = pc[2];

        const Real x0 = x - vel[0] * dt;
        const Real y0 = y - vel[1] * dt;
        const Real z0 = z - vel[2] * dt;

        if (!(inside01(x0) && inside01(y0) && inside01(z0)))
          continue;

        const Real uex =
          std::sin(pi * x0) * std::sin(pi * y0) * std::sin(pi * z0);

        const Real uhv = uh(p);

        const Real diff = uhv - uex;
        err2 += diff * diff;
        ref2 += uex * uex;
        ++n;
      }

      ASSERT_GT(n, 0u);

      const Real err = std::sqrt(err2);
      const Real ref = std::sqrt(ref2);
      EXPECT_TRUE(err <= atol || (ref > 0 ? err / ref <= rtol : err <= atol))
          << "err=" << err << " ref=" << ref << " n=" << n;
    }
  };

  using ManufacturedAdvection3DTest_10 = ManufacturedAdvection3DTest<15, 15, 15>;

  TEST_P(ManufacturedAdvection3DTest_10, ConstantVelocity_OneStep_L2Interior)
  {
    auto mesh = this->getMesh();
    P1 vh(mesh);

    Math::SpatialVector<Real> vel{{0.15, -0.25, 0.20}};
    const Real dt = 0.9;

    auto velocity = VectorFunction{
      RealFunction([&](const Point&) { return vel[0]; }),
      RealFunction([&](const Point&) { return vel[1]; }),
      RealFunction([&](const Point&) { return vel[2]; })
    };

    auto u0 = sin(Math::Constants::pi() * F::x)
            * sin(Math::Constants::pi() * F::y)
            * sin(Math::Constants::pi() * F::z);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Advection::Lagrangian lagrangian(u, v, u0, velocity);
    lagrangian.step(dt);

    const auto& uh = u.getSolution();

    const Real atol = 5e-3;
    const Real rtol = 5e-2;
    this->checkL2CentroidError(mesh, uh, vel, dt, atol, rtol);
  }

  TEST_P(ManufacturedAdvection3DTest_10, ConstantVelocity_TwoHalfSteps_vs_OneFullStep)
  {
    auto mesh = this->getMesh();
    P1 vh(mesh);

    Math::SpatialVector<Real> vel{{-0.1, 0.05, 0.08}};
    const Real dt = 0.08;

    auto velocity = VectorFunction{
      RealFunction([&](const Point&) { return vel[0]; }),
      RealFunction([&](const Point&) { return vel[1]; }),
      RealFunction([&](const Point&) { return vel[2]; })
    };

    auto u0 = sin(Math::Constants::pi() * F::x)
            * sin(Math::Constants::pi() * F::y)
            * sin(Math::Constants::pi() * F::z);

    // two half steps
    TrialFunction u2(vh);
    TestFunction  v2(vh);
    {
      Advection::Lagrangian L(u2, v2, u0, velocity);
      L.step(0.5 * dt);
      L.step(0.5 * dt);
    }

    // one full step
    TrialFunction u1(vh);
    TestFunction  v1(vh);
    {
      Advection::Lagrangian L(u1, v1, u0, velocity);
      L.step(dt);
    }

    const auto& uh1 = u1.getSolution();
    const auto& uh2 = u2.getSolution();

    const size_t cd = mesh.getDimension();
    Math::SpatialVector<Real> rc{{0.5, 0.5, 0.5}};

    Real num2 = 0.0, den2 = 0.0;
    size_t n = 0;

    auto inside01 = [](Real z) { return z >= 0.0 && z <= 1.0; };

    for (Index c = 0; c < mesh.getPolytopeCount(cd); ++c)
    {
      auto itc = mesh.getPolytope(cd, c);
      const auto& cell = *itc;

      Geometry::Point p(cell, rc);
      const auto pc = p.getPhysicalCoordinates();
      const Real x = pc[0], y = pc[1], z = pc[2];

      const Real x0 = x - vel[0] * dt;
      const Real y0 = y - vel[1] * dt;
      const Real z0 = z - vel[2] * dt;
      if (!(inside01(x0) && inside01(y0) && inside01(z0)))
        continue;

      const Real a = uh1(p);
      const Real b = uh2(p);
      const Real d = a - b;
      num2 += d * d;
      den2 += 0.5 * (a * a + b * b);
      ++n;
    }

    ASSERT_GT(n, 0u);
    const Real rel = std::sqrt(num2) / (den2 > 0 ? std::sqrt(den2) : 1.0);
    EXPECT_LE(rel, 5e-2) << "two half-steps diverge from one full step too much";
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    ManufacturedAdvection3DTest_10,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge)
  );
}
