/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using namespace Rodin::External;

static constexpr Scalar hmax = 0.01; // Maximal size of a triangle's edge
static constexpr Scalar hmin = 0.1 * hmax;

static constexpr Attribute dQ = 2; // Attribute of box boundary
static constexpr Attribute dCurrent = 4; // Attribute of box boundary
static constexpr Attribute dGround = 5; // Attribute of box boundary

static const Math::Vector x0{{0.5, 0.5}}; // Center of domain
static constexpr Scalar epsilon = 0.05; // Radius of B_e(x_0)
static_assert(epsilon > 0);

static constexpr Scalar pi = Math::Constants::pi();

// static constexpr Scalar gamma_ek =
//   (sectorArea * gammaA1 + (inhomogeinityArea - sectorArea) * gammaA2) / inhomogeinityArea;

const constexpr Scalar m = 50;
static constexpr Scalar gamma_ek = 1e-12;

static constexpr Scalar R0 = 0.2; // Radius of B_R(x_0)
static constexpr Scalar R1 = R0 + 10 * hmax; // Radius of B_R(x_0)

static Solver::SparseLU solver;

int main(int, char**)
{
  Alert::Info() << "Epsilon: " << epsilon
                << Alert::NewLine
                << "m: " << m
                << Alert::NewLine
                << "Inhomogeinity constant: " << gamma_ek
                << Alert::Raise;

  // Define mesh
  MMG::Mesh mesh;
  mesh.load("Q1.medit.mesh", IO::FileFormat::MEDIT);
  mesh.save("Q.mesh");

  // // Refine mesh
  // {
  //   P1 fes(mesh);
  //   MMG::ScalarGridFunction sizeMap(fes);
  //   sizeMap =
  //     [&](const Geometry::Point& p)
  //     {
  //       const Scalar r = (p.getCoordinates() - x0).norm();
  //       return r;
  //     };
  //   sizeMap *= (hmax - hmin);
  //   sizeMap += hmin;
  //   MMG::Adapt().setAngleDetection(false).setHMax(hmax).setHMin(hmin).adapt(mesh, sizeMap);
  // }
  // mesh.save("Q1.medit.mesh", IO::FileFormat::MEDIT);
  // std::exit(1);

  // Define finite element spaces
  P1 vh(mesh);
  P1 gh(mesh, mesh.getSpaceDimension());

  // Define conductivity
  ScalarFunction gamma =
    [&](const Point& p)
    {
      return 2 + sin(pi * m * p.x()) * cos(pi * m * p.y());
    };

  ScalarFunction gamma_e =
    [&](const Point& p)
    {
      const Scalar r = (p.getCoordinates() - x0).norm();
      if (r > epsilon)
        return gamma(p);
      else
        return gamma_ek;
    };

  GridFunction conductivity(vh);
  conductivity = gamma;
  conductivity.save("Conductivity.gf");

  conductivity = gamma_e;
  conductivity.save("Conductivity_E.gf");

  // Define variational problems
  TrialFunction u(vh);
  TestFunction  v(vh);
  ScalarFunction phi = 1;

  const auto f =
    [&](const Scalar& x)
    {
      if (x > 0)
        return std::exp(-1.0 / x);
      else
        return 0.0;
    };

  // Define source
  const auto g =
    [&](const Scalar& x)
    {
      return f(x) / (f(x) + f(1 - x));
    };

  ScalarFunction h = 1;
  //   [&](const Geometry::Point& p)
  //   {
  //     const Scalar r = (p.getCoordinates() - x0).norm();
  //     return g((r - R0) / (R1 - R0));
  //   };

  Problem poisson(u, v);
  poisson = Integral(gamma * Grad(u), Grad(v))
          + DirichletBC(u, phi).on(dCurrent)
          + DirichletBC(u, ScalarFunction(0)).on(dGround);

  Problem perturbed(u, v);
  perturbed = Integral(gamma_e * Grad(u), Grad(v))
            + DirichletBC(u, phi).on(dCurrent)
            + DirichletBC(u, ScalarFunction(0)).on(dGround);

  // Solve the background problem
  Alert::Info() << "Solving background equation." << Alert::Raise;
  poisson.assemble().solve(solver);
  const auto u0 = u.getSolution();

  // u0.save("Background.gf");

  Alert::Info() << "Computing its gradient." << Alert::Raise;
  GridFunction g0(gh);
  g0 = Grad(u0);
  g0.save("BackgroundGradient.gf");

  // Solve the perturbed problem
  Alert::Info() << "Solving perturbed equation." << Alert::Raise;

  perturbed.assemble().solve(solver);
  const auto u_e = u.getSolution();
  GridFunction g_e(gh);
  u_e.save("Perturbed.gf");

  Alert::Info() << "Computing its gradient." << Alert::Raise;
  g_e = Grad(u_e);
  g_e.save("PerturbedGradient.gf");

  GridFunction diff(vh);
  GridFunction gdiff(gh);
  //diff = u0 - u_e;
  diff = Frobenius(Grad(u0) - Grad(u_e));
  gdiff.save("GDiff.gf");
  diff.setWeights();
  diff.save("Diff.gf");

  return 0;
}


