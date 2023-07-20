/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static constexpr Attribute A1 = 1; // Attribute of triangle
static constexpr Attribute A2 = 2; // Attribute of complement
static constexpr Attribute dQ = 2; // Attribute of box boundary

static const Math::Vector x0{{0.5, 0.5}}; // Center of domain
static constexpr Scalar epsilon = 0.01; // Radius of B_e(x_0)
static constexpr Scalar R0 = 0.1; // Radius of B_R(x_0)
static constexpr Scalar R1 = 0.2; // Radius of B_R(x_0)
static_assert(R1 > R0);
static_assert(R0 > epsilon);
static_assert(epsilon > 0);

static constexpr Scalar pi = Math::Constants::pi();
static constexpr Scalar sectorArea = pi * epsilon * epsilon / 4;
static constexpr Scalar inhomogeinityArea = pi * epsilon * epsilon;

static constexpr Scalar gammaA1 = 10;
static constexpr Scalar gammaA2 = 1;
static constexpr Scalar inhomogeinityConstant =
  (sectorArea * gammaA1 + (inhomogeinityArea - sectorArea) * gammaA2) / inhomogeinityArea;

static Solver::SparseLU solver;

int main(int, char**)
{
  Alert::Info() << "Epsilon: " << epsilon
                << Alert::NewLine
                << "(R0, R1): (" << R0 << ", " << R1 << ")"
                << Alert::NewLine
                << "Area of sector: " << sectorArea
                << Alert::NewLine
                << "Area of inhomogeinity: " << inhomogeinityArea
                << Alert::NewLine
                << "A1 Conductivity: " << gammaA1
                << Alert::NewLine
                << "A2 Conductivity: " << gammaA2
                << Alert::NewLine
                << "Inhomogeinity constant: " << inhomogeinityConstant
                << Alert::Raise;

  // Define mesh
  Mesh mesh;
  mesh.load("Q.medit.mesh", IO::FileFormat::MEDIT);
  mesh.save("Q.mesh");

  // Define finite element spaces
  P1 vh(mesh);
  P1 gh(mesh, mesh.getSpaceDimension());

  // Define source
  const auto f =
    [&](const Scalar& x)
    {
      if (x > 0)
        return std::exp(-1.0 / x);
      else
        return 0.0;
    };

  const auto g =
    [&](const Scalar& x)
    {
      return f(x) / (f(x) + f(1 - x));
    };

  ScalarFunction h =
    [&](const Geometry::Point& p)
    {
      const Scalar r = (p.getCoordinates() - x0).norm();
      return g((r - epsilon) / (R1 - R0));
    };

  GridFunction source(vh);
  source = h;
  source.save("Source.gf");

  // Define conductivity
  // ScalarFunction gamma =
  //   [](const Point& p)
  //   {
  //     return 2 + sin(pi * m * p.x()) * cos(pi * m * p.y());
  //   };

  ScalarFunction gamma =
    [](const Geometry::Point& p)
    {
      if (p.getPolytope().getAttribute() == A1)
        return gammaA1;
      else
        return gammaA2;
    };

  ScalarFunction gamma_e =
    [&](const Geometry::Point& p)
    {
      const Scalar r = (p.getCoordinates() - x0).norm();
      if (r > epsilon)
        return gamma(p);
      else
        return inhomogeinityConstant;
    };

  // Define variational problems
  TrialFunction u(vh);
  TestFunction  v(vh);
  ScalarFunction phi = 0;

  Problem poisson(u, v);
  poisson = Integral(gamma * Grad(u), Grad(v))
          - Integral(h, v)
          + DirichletBC(u, phi).on(dQ);

  Problem perturbed(u, v);
  perturbed = Integral(gamma_e * Grad(u), Grad(v))
            - Integral(h, v)
            + DirichletBC(u, phi).on(dQ);

  // Define gradient grid function
  GridFunction gradient(gh);

  // Solve the background problem
  Alert::Info() << "Solving background equation." << Alert::Raise;
  poisson.solve(solver);
  const auto u0 = u.getSolution();
  u0.save("Background.gf");

  Alert::Info() << "Computing its gradient." << Alert::Raise;
  gradient = Grad(u.getSolution());
  gradient.save("BackgroundGradient.gf");

  // Solve the perturbed problem
  Alert::Info() << "Solving perturbed equation." << Alert::Raise;

  Solver::SparseLU solver;
  perturbed.solve(solver);
  const auto u_e = u.getSolution();
  u_e.save("Perturbed.gf");

  Alert::Info() << "Computing its gradient." << Alert::Raise;
  gradient = Grad(u.getSolution());
  gradient.save("PerturbedGradient.gf");

  GridFunction diff(vh);
  diff = Pow(u_e - u0, 2);
  diff.setWeights();
  const Scalar err = sqrt(Integral(diff).compute());

  Alert::Info() << "Error: " << err << Alert::Raise;

  return 0;
}


