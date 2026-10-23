/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief TEMPORARY. Convergence of the stabilized P1--P1 Stokes discretization.
 *
 * The Kelvin-ball example discretizes Stokes with equal-order P1 velocity and
 * pressure, stabilized by a pressure-gradient term
 * @f$ \tau \int \nabla p \cdot \nabla q @f$. The existing manufactured tests
 * cover the unstabilized H1--P1 pair against an exact affine solution, so they
 * say nothing about the rates of this pair, nor about which mesh size belongs
 * in @f$ \tau @f$.
 *
 * These tests measure the rates against a manufactured solution, for the two
 * candidate stabilizations:
 * - a single size taken from the background mesh, @f$ \tau = \alpha h^2/\mu @f$;
 * - the size of each cell, @f$ \tau_K = \alpha h_K^2/\mu @f$.
 * Each is measured on a uniform mesh, where the two agree up to a constant,
 * and on a graded mesh, where they do not.
 *
 * The manufactured solution is divergence-free by construction, being the curl
 * of @f$ (\sin \pi y \sin \pi z, \sin \pi z \sin \pi x, \sin \pi x \sin \pi y) @f$,
 * and satisfies @f$ -\Delta u = 2\pi^2 u @f$, so that
 * @f$ f = 2 \pi^2 \mu u + \nabla p @f$.
 */

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Manufactured::StokesStabilized
{
  static constexpr Real Mu = 1;
  static constexpr Real StabilizationFactor = 0.05;

  /// @brief Which mesh size enters the stabilization.
  enum class Stabilization
  {
    Background,
    Cell
  };

  /// @brief Edge of the regular tetrahedron of the same measure as @p cell.
  static Real cellSize(const Polytope& cell)
  {
    static const Real scale = std::cbrt(6 * std::sqrt(Real(2)));
    return scale * std::cbrt(cell.getMeasure());
  }

  static Math::SpatialVector<Real> velocity(const Geometry::Point& point)
  {
    const auto x = point.getPhysicalCoordinates();
    Math::SpatialVector<Real> value(3);
    value(0) =
      M_PI * std::sin(M_PI * x(0)) * (std::cos(M_PI * x(1)) - std::cos(M_PI * x(2)));
    value(1) =
      M_PI * std::sin(M_PI * x(1)) * (std::cos(M_PI * x(2)) - std::cos(M_PI * x(0)));
    value(2) =
      M_PI * std::sin(M_PI * x(2)) * (std::cos(M_PI * x(0)) - std::cos(M_PI * x(1)));
    return value;
  }

  static Real pressure(const Geometry::Point& point)
  {
    const auto x = point.getPhysicalCoordinates();
    return std::cos(M_PI * x(0)) * std::sin(M_PI * x(1)) * std::sin(M_PI * x(2));
  }

  static Math::SpatialVector<Real> force(const Geometry::Point& point)
  {
    const auto x = point.getPhysicalCoordinates();
    Math::SpatialVector<Real> value = 2 * M_PI * M_PI * Mu * velocity(point);
    value(0) +=
      -M_PI * std::sin(M_PI * x(0)) * std::sin(M_PI * x(1)) * std::sin(M_PI * x(2));
    value(1) +=
      M_PI * std::cos(M_PI * x(0)) * std::cos(M_PI * x(1)) * std::sin(M_PI * x(2));
    value(2) +=
      M_PI * std::cos(M_PI * x(0)) * std::sin(M_PI * x(1)) * std::cos(M_PI * x(2));
    return value;
  }

  /// @brief Unit cube of @p n points per edge, optionally graded towards x = 0.
  static Mesh<Context::Local> makeMesh(size_t n, bool graded)
  {
    Mesh<Context::Local> mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(Real(1) / static_cast<Real>(n - 1));
    if (graded)
    {
      // A smooth monotone map keeps the mesh valid while making the cell size
      // vary by the refinement factor across the cube.
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      {
        auto coordinates = mesh.getVertexCoordinates(vertex);
        coordinates(0) = coordinates(0) * coordinates(0);
        mesh.setVertexCoordinates(vertex, coordinates);
      }
    }
    mesh.getConnectivity().compute(2, 3);
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    return mesh;
  }

  struct Errors
  {
      Real velocity = 0;
      Real pressure = 0;
  };

  /// @brief Solves the stabilized system and returns the two @f$ L^2 @f$ errors.
  static Errors solve(size_t n, bool graded, Stabilization stabilization)
  {
    Mesh mesh = makeMesh(n, graded);
    const Real h = Real(1) / static_cast<Real>(n - 1);

    P1 vh(mesh, mesh.getSpaceDimension());
    P1 ph(mesh);
    P0g gauge(mesh);

    VectorFunction exactVelocity(mesh.getSpaceDimension(), velocity);
    RealFunction exactPressure(pressure);
    VectorFunction exactForce(mesh.getSpaceDimension(), force);

    const Real background = StabilizationFactor * h * h / Mu;
    RealFunction tau([stabilization, background](const Geometry::Point& point) {
      if (stabilization == Stabilization::Background)
        return background;
      const Real size = cellSize(point.getPolytope());
      return StabilizationFactor * size * size / Mu;
    });

    TrialFunction u(vh);
    TrialFunction p(ph);
    TrialFunction lambda(gauge);
    TestFunction v(vh);
    TestFunction q(ph);
    TestFunction eta(gauge);

    const auto Du = Real(0.5) * (Jacobian(u) + Jacobian(u).T());
    const auto Dv = Real(0.5) * (Jacobian(v) + Jacobian(v).T());

    Problem stokes(u, p, lambda, v, q, eta);
    stokes = Integral(Real(2) * Mu * Du, Dv) - Integral(p, Div(v)) - Integral(Div(u), q) -
      Integral(tau * Grad(p), Grad(q)) + Integral(lambda, q) + Integral(p, eta) -
      Integral(exactForce, v) + DirichletBC(u, exactVelocity);
    stokes.assemble();

    Solver::SparseLU solver(stokes);
    solver.solve();
    EXPECT_TRUE(solver.success());

    P1 sh(mesh);
    GridFunction velocityError(sh);
    velocityError = Pow(Frobenius(u.getSolution() - exactVelocity), 2);

    GridFunction mean(sh);
    mean = p.getSolution() - exactPressure;
    const Real shift = Integral(mean).compute() / mesh.getMeasure(mesh.getDimension());
    GridFunction pressureError(sh);
    pressureError = Pow((p.getSolution() - shift) - exactPressure, 2);

    return {std::sqrt(Integral(velocityError).compute()),
      std::sqrt(Integral(pressureError).compute())};
  }

  static Real rate(Real coarse, Real fine)
  {
    return std::log(coarse / fine) / std::log(Real(2));
  }

  static void expectConvergence(bool graded, Stabilization stabilization)
  {
    const Errors coarse = solve(5, graded, stabilization);
    const Errors medium = solve(9, graded, stabilization);
    const Errors fine = solve(17, graded, stabilization);

    const Real velocityRate = rate(medium.velocity, fine.velocity);
    const Real pressureRate = rate(medium.pressure, fine.pressure);
    std::cout << "[RATES] graded=" << graded << " stabilization="
              << (stabilization == Stabilization::Cell ? "cell" : "background")
              << " eU=" << coarse.velocity << ',' << medium.velocity << ','
              << fine.velocity << " rateU=" << velocityRate << " eP=" << coarse.pressure
              << ',' << medium.pressure << ',' << fine.pressure
              << " rateP=" << pressureRate << std::endl;

    // The pair is stabilized, not consistent, so the pressure carries the
    // O(tau) perturbation: first order is what it can deliver.
    EXPECT_GT(velocityRate, Real(1.5));
    EXPECT_GT(pressureRate, Real(0.8));
    EXPECT_LT(fine.velocity, medium.velocity);
    EXPECT_LT(fine.pressure, medium.pressure);
  }

  /// @brief Rates on a uniform mesh with a single background size in the stabilization.
  TEST(Manufactured_StokesStabilized, UniformMeshWithBackgroundSize)
  {
    expectConvergence(false, Stabilization::Background);
  }

  /// @brief Rates on a uniform mesh with the size of each cell in the stabilization.
  TEST(Manufactured_StokesStabilized, UniformMeshWithCellSize)
  {
    expectConvergence(false, Stabilization::Cell);
  }

  /// @brief Rates on a graded mesh with a single background size in the stabilization.
  TEST(Manufactured_StokesStabilized, GradedMeshWithBackgroundSize)
  {
    expectConvergence(true, Stabilization::Background);
  }

  /// @brief Rates on a graded mesh with the size of each cell in the stabilization.
  TEST(Manufactured_StokesStabilized, GradedMeshWithCellSize)
  {
    expectConvergence(true, Stabilization::Cell);
  }
}
