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

static constexpr Real hmax = 0.01; // Maximal size of a triangle's edge
static constexpr Real hmin = 0.1 * hmax;

static constexpr Attribute dQ = 2; // Attribute of box boundary
static constexpr Attribute dCurrent = 4; // Attribute of box boundary
static constexpr Attribute dGround = 5; // Attribute of box boundary

static const Math::Vector<Real> x0{{0.5, 0.5}}; // Center of domain
static constexpr Real epsilon = 0.1; // Radius of B_e(x_0)
static_assert(epsilon > 0);

static constexpr Real pi = Math::Constants::pi();

const constexpr Real angle = M_PI / 4;
const constexpr Real waveNumber = 40;
const constexpr Real m = 40;
static constexpr Real gamma_ek = 2;

static constexpr Real R0 = 0.2; // Radius of B_R(x_0)
static constexpr Real R1 = R0 + 10 * hmax; // Radius of B_R(x_0)

int main(int, char**)
{
  Alert::Info() << "Epsilon: " << epsilon
                << Alert::NewLine
                << "m: " << m
                << Alert::NewLine
                << "waveNumber: " << waveNumber
                << Alert::NewLine
                << "Inhomogeinity constant: " << gamma_ek
                << Alert::Raise;

  // Define mesh
  Mesh mesh;
  mesh.load("Q1.medit.mesh", IO::FileFormat::MEDIT);
  mesh.save("Q.mesh");
  mesh.getConnectivity().compute(1, 2);

  // Define finite element spaces
  P1 vh(mesh);
  P1 gh(mesh, mesh.getSpaceDimension());

  // Define oscillatory screen
  RealFunction h =
    [&](const Point& p)
    {
      return 2 + sin(2 * pi * m * p.x()) * sin(2 * pi * m * p.y());
    };

  // Define conductivity
  RealFunction gamma = 1;

  RealFunction gamma_e =
    [&](const Point& p)
    {
      const Real r = (p.getCoordinates() - x0).norm();
      if (r > epsilon)
        return gamma(p);
      else
        return gamma_ek;
    };

  GridFunction screen(vh);
  screen = h;
  screen.save("Screen.gf");

  GridFunction conductivity(vh);
  conductivity = gamma;
  conductivity.save("Conductivity.gf");

  conductivity = gamma_e;
  conductivity.save("Conductivity_E.gf");

  VectorFunction xi = { std::cos(angle), std::sqrt(angle) };
  RealFunction phi =
    [&](const Point& p)
    { return cos(waveNumber * p.getCoordinates().dot(xi(p))); };

  GridFunction miaow(vh);
  miaow = phi;
  miaow.save("PlaneWave.gf");

  // Define variational problems
  TrialFunction u(vh);
  TestFunction  v(vh);

  Problem helmholtz(u, v);
  helmholtz = Integral(gamma * Grad(u), Grad(v))
            - Integral(waveNumber * waveNumber * h * u, v)
            + DirichletBC(u, phi);

  Problem perturbed(u, v);
  perturbed = Integral(gamma_e * Grad(u), Grad(v))
            - Integral(waveNumber * waveNumber * h * u, v)
            + DirichletBC(u, phi);

  // Solve the background problem
  Alert::Info() << "Solving background equation." << Alert::Raise;
  Solver::SparseLU(helmholtz).solve();
  const auto u0 = u.getSolution();

  u0.save("Background.gf");

  Alert::Info() << "Computing its gradient." << Alert::Raise;
  GridFunction g0(gh);
  g0 = Grad(u0);
  g0.save("BackgroundGradient.gf");

  // Solve the perturbed problem
  Alert::Info() << "Solving perturbed equation." << Alert::Raise;

  Solver::SparseLU(perturbed).solve();
  const auto u_e = u.getSolution();
  GridFunction g_e(gh);
  u_e.save("Perturbed.gf");

  Alert::Info() << "Computing its gradient." << Alert::Raise;
  g_e = Grad(u_e);
  g_e.save("PerturbedGradient.gf");

  RealFunction chi =
    [&](const Point& p)
    { return Real((p.getCoordinates() - x0).norm() > 0.25); };

  GridFunction diff(vh);
  diff = chi * Pow(u0 - u_e, 2);
  diff.setWeights();

  GridFunction indicator(vh);
  indicator = chi;
  indicator.save("Indicator.gf");

  const Real error = sqrt(Integral(diff).compute());

  Alert::Info() << "L2 Error: " << error << Alert::Raise;

  return 0;
}




