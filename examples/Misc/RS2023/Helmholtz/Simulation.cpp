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

static constexpr Scalar hmax = 0.01; // Maximal size of a triangle's edge
static constexpr Scalar hmin = 0.1 * hmax;

static constexpr Attribute dQ = 2; // Attribute of box boundary
static constexpr Attribute dCurrent = 4; // Attribute of box boundary
static constexpr Attribute dGround = 5; // Attribute of box boundary

static const Math::Vector x0{{0.5, 0.5}}; // Center of domain

static constexpr Scalar pi = Math::Constants::pi();

static constexpr Scalar R0 = 0.2; // Radius of B_R(x_0)
static constexpr Scalar R1 = R0 + 10 * hmax; // Radius of B_R(x_0)

static Solver::SparseLU solver;

struct Data
{
  const Scalar m;
  const Scalar epsilon;
  const Scalar waveNumber;
  const Scalar conductivity;
};

int main(int, char**)
{
  // Load mesh
  Mesh mesh;
  mesh.load("Q1.medit.mesh", IO::FileFormat::MEDIT);

  // Define evaluation grid
  constexpr size_t N = 20;
  Math::Vector m_r = Math::Vector::LinSpaced(N, 0, 0.5 * 1. / hmax);
  Math::Vector epsilon_r = Math::Vector::LinSpaced(N, hmax, 0.2);
  Math::Vector waveNumber_r = Math::Vector::LinSpaced(25, 0.1, 25);
  Math::Vector conductivity_r{{ 1e-12, 1.0, 1e12 }};

  std::vector<Data> grid;
  for (const Scalar m : m_r)
    for (const Scalar epsilon : epsilon_r)
      for (const Scalar waveNumber : waveNumber_r)
        for (const Scalar g : conductivity_r)
          grid.push_back({ m, epsilon, waveNumber, g });

  std::stringstream filename;
  filename << "L2_Grid_HMax" << std::setw(4) << hmax << ".live.csv";
  std::ofstream out(filename.str());
  out << "m,epsilon,waveNumber,conductivity,error";
  for (const auto& data : grid)
  {
    Alert::Info() << "Epsilon: " << data.epsilon
                  << Alert::NewLine
                  << "m: " << data.m
                  << Alert::NewLine
                  << "waveNumber: " << data.waveNumber
                  << Alert::NewLine
                  << "Inhomogeinity constant: " << data.conductivity
                  << Alert::Raise;

    // Define finite element spaces
    P1 vh(mesh);

    // Define oscillatory screen
    ScalarFunction h =
      [&](const Point& p)
      {
        return 2 + sin(2 * pi * data.m * p.x()) * sin(2 * pi * data.m * p.y());
      };

    // Define conductivity
    ScalarFunction gamma = 1;

    ScalarFunction gamma_e =
      [&](const Point& p)
      {
        const Scalar r = (p.getCoordinates() - x0).norm();
        if (r > data.epsilon)
          return gamma(p);
        else
          return data.conductivity;
      };

    // Define boundary data
    ScalarFunction phi = 1;

    // Define variational problems
    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem helmholtz(u, v);
    helmholtz = Integral(gamma * Grad(u), Grad(v))
              - Integral(data.waveNumber * data.waveNumber * h * u, v)
              + DirichletBC(u, phi).on(dCurrent)
              + DirichletBC(u, ScalarFunction(0)).on(dGround);

    Problem perturbed(u, v);
    perturbed = Integral(gamma_e * Grad(u), Grad(v))
              - Integral(data.waveNumber * data.waveNumber * h * u, v)
              + DirichletBC(u, phi).on(dCurrent)
              + DirichletBC(u, ScalarFunction(0)).on(dGround);

    // Solve the background problem
    Alert::Info() << "Solving background equation." << Alert::Raise;
    helmholtz.solve(solver);
    const auto u0 = std::move(u.getSolution());
    Alert::Success() << "Done." << Alert::Raise;

    Alert::Info() << "Solving perturbed equation." << Alert::Raise;
    perturbed.solve(solver);
    const auto ue = std::move(u.getSolution());
    Alert::Success() << "Done." << Alert::Raise;

    mesh.save("Q.mesh");
    u0.save("u0.gf");
    ue.save("ue.gf");

    GridFunction diff(vh);
    diff = Pow(u0 - ue, 2);
    diff.setWeights();
    const Scalar error = Integral(diff);

    Alert::Success() << "L2 Error: " << error << Alert::Raise;

    out.flush();
    std::exit(1);
  }

  return 0;
}




