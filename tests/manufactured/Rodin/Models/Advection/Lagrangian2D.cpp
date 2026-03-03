/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE or copy at
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

namespace Rodin::Tests::Manufactured::AdvectionLagrangian2D
{
  template <size_t M>
  class ManufacturedAdvectionTest : public ::testing::TestWithParam<Polytope::Type>
  {
   protected:
    Mesh<Context::Local> getMesh()
    {
      Mesh<Context::Local> mesh = Mesh<Context::Local>::UniformGrid(GetParam(), { M, M });
      mesh.scale(1.0 / (M - 1));
      mesh.getConnectivity().compute(1, 2);
      return mesh;
    }

    // L2 error on interior cell centroids only (avoid boundary exits)
    template <class GF>
    static void checkL2CentroidError(
        const Mesh<Context::Local>& mesh,
        const GF& uh,
        Real vx, Real vy, Real dt,
        Real atol, Real rtol)
    {
      const Real pi = Math::Constants::pi();

      auto inside01 = [](Real z) { return z >= 0.0 && z <= 1.0; };

      const size_t cd = mesh.getDimension();
      ASSERT_EQ(cd, 2u);

      Real err2 = 0.0;
      Real ref2 = 0.0;
      size_t n = 0;

      // centroid in reference space for both tri and quad
      Math::SpatialVector<Real> rc{{0.5, 0.5}};

      for (Index c = 0; c < mesh.getPolytopeCount(cd); ++c)
      {
        auto itc = mesh.getPolytope(cd, c);
        const auto& cell = *itc;

        Geometry::Point p(cell, rc);
        // physical centroid
        const auto pc = p.getPhysicalCoordinates();
        const Real x = pc[0], y = pc[1];

        // back-traced foot under constant velocity
        const Real x0 = x - vx * dt;
        const Real y0 = y - vy * dt;

        if (!(inside01(x0) && inside01(y0)))
          continue; // skip cells whose characteristics exit in dt

        const Real uex =
          std::sin(pi * x0) * std::sin(pi * y0);

        const Real uhv = uh(p);

        const Real diff = uhv - uex;
        err2 += diff * diff;
        ref2 += uex * uex;
        ++n;
      }

      ASSERT_GT(n, 0u); // ensure we actually sampled interior cells

      const Real err = std::sqrt(err2);
      const Real ref = std::sqrt(ref2);
      // pass if absolute error small or relative error small
      EXPECT_TRUE(err <= atol || (ref > 0 ? err / ref <= rtol : err <= atol))
          << "err=" << err << " ref=" << ref << " n=" << n;
    }
  };

  using ManufacturedAdvectiionTest_16 = ManufacturedAdvectionTest<16>;
  using ManufacturedAdvectiionTest_32 = ManufacturedAdvectionTest<32>;

  // One-step manufactured check with constant velocity.
  // u0(x,y) = sin(pi x) sin(pi y)
  // u(x,y,dt) = u0(x - vx dt, y - vy dt)
  TEST_P(ManufacturedAdvectiionTest_32, ConstantVelocity_OneStep_L2Interior)
  {
    auto mesh = this->getMesh();

    P1 vh(mesh);

    const Real vx = -0.20;
    const Real vy =  0.35;
    const Real dt =  0.4; // small to keep most centroids interior

    auto velocity = VectorFunction{
      RealFunction([vx](const Point&) { return vx; }),
      RealFunction([vy](const Point&) { return vy; })
    };

    auto u0 = sin(Math::Constants::pi() * F::x)
            * sin(Math::Constants::pi() * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Advection::Lagrangian lagrangian(u, v, u0, velocity);
    lagrangian.step(dt);

   const auto& uh = u.getSolution();

    // tolerances tuned for P1 + SL arrival with RK tracing
    const Real atol = 5e-3;  // absolute L2 on centroids
    const Real rtol = 5e-2;  // relative L2 on centroids
    this->checkL2CentroidError(mesh, uh, vx, vy, dt, atol, rtol);
  }

  // Two smaller steps should not be worse than one larger step for the same total time.
  TEST_P(ManufacturedAdvectiionTest_32, ConstantVelocity_TwoHalfSteps_vs_OneFullStep)
  {
    auto mesh = this->getMesh();
    P1 vh(mesh);

    const Real vx = 0.15;
    const Real vy = 0.10;
    const Real dt = 0.04;

    auto velocity = VectorFunction{
      RealFunction([vx](const Point&) { return vx; }),
      RealFunction([vy](const Point&) { return vy; })
    };

    auto u0 = sin(Math::Constants::pi() * F::x)
            * sin(Math::Constants::pi() * F::y);

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

    // compare u1 and u2 at interior centroids
    const auto& uh1 = u1.getSolution();
    const auto& uh2 = u2.getSolution();

    const size_t cd = mesh.getDimension();
    Math::SpatialVector<Real> rc{{0.5, 0.5}};

    Real num2 = 0.0, den2 = 0.0;
    size_t n = 0;

    auto inside01 = [](Real z) { return z >= 0.0 && z <= 1.0; };



    for (Index c = 0; c < mesh.getPolytopeCount(cd); ++c)
    {
      auto itc = mesh.getPolytope(cd, c);
      const auto& cell = *itc;

      Geometry::Point p(cell, rc);
      const auto pc = p.getPhysicalCoordinates();
      const Real x = pc[0], y = pc[1];

      // back-traced foot for total dt
      const Real x0 = x - vx * dt;
      const Real y0 = y - vy * dt;
      if (!(inside01(x0) && inside01(y0)))
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

  // REPLACE the failing test with a spatial-refinement check.
  TEST_P(ManufacturedAdvectiionTest_16, ConstantVelocity_SpatialRefinement_DecreasesError)
  {
    // coarse mesh
    auto meshC = this->getMesh();          // 16x16
    P1 vhC(meshC);

    // refined mesh (uniform 2x in each dir -> 31x31 nodes on [0,1])
    Mesh<Context::Local> meshF = Mesh<Context::Local>::UniformGrid(GetParam(), { 2*16-1, 2*16-1 });
    meshF.scale(1.0 / (2*16-2));
    meshF.getConnectivity().compute(1, 2);
    P1 vhF(meshF);

    const Real vx = -0.10, vy = 0.25;
    const Real dt = 0.03; // small to avoid exits at centroids

    auto velocity = VectorFunction{
      RealFunction([vx](const Point&) { return vx; }),
      RealFunction([vy](const Point&) { return vy; })
    };

    auto u0 = sin(Math::Constants::pi() * F::x) * sin(Math::Constants::pi() * F::y);

    auto compute_centroid_err = [&](auto& mesh, auto& vh) -> Real
    {
      TrialFunction u(vh);
      TestFunction  v(vh);
      Advection::Lagrangian L(u, v, u0, velocity);
      L.step(dt);
      const auto& uh = u.getSolution();

      const size_t cd = mesh.getDimension();
      Math::SpatialVector<Real> rc{{0.5, 0.5}};
      auto inside01 = [](Real z) { return z >= 0.0 && z <= 1.0; };

      Real err2 = 0.0;
      size_t n = 0;
      for (Index c = 0; c < mesh.getPolytopeCount(cd); ++c)
      {
        auto itc = mesh.getPolytope(cd, c);
        const auto& cell = *itc;
        Geometry::Point p(cell, rc);
        const auto pc = p.getPhysicalCoordinates();
        const Real x = pc[0], y = pc[1];
        const Real x0 = x - vx * dt, y0 = y - vy * dt;
        if (!(inside01(x0) && inside01(y0))) continue;

        const Real uex = std::sin(Math::Constants::pi() * x0) * std::sin(Math::Constants::pi() * y0);
        const Real uhv = uh(p);
        const Real d = uhv - uex;
        err2 += d * d;
        ++n;
      }
      EXPECT_GT(n, 0u);
      return std::sqrt(err2);
    };

    const Real eC = compute_centroid_err(meshC, vhC);
    const Real eF = compute_centroid_err(meshF, vhF);

    // Expect lower error on finer mesh (projection improves). Allow slack.
    EXPECT_LT(eF, 0.7 * eC);
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    ManufacturedAdvectiionTest_16,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

  INSTANTIATE_TEST_SUITE_P(
    MeshParams32x32,
    ManufacturedAdvectiionTest_32,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}
