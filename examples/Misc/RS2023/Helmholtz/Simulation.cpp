/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <chrono>
#include <fstream>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static constexpr Scalar hmax = 0.01; // Maximal size of a triangle's edge
static constexpr Scalar hmin = 0.01 * hmax;

static constexpr Attribute dQ = 2; // Attribute of box boundary
static constexpr Attribute dCurrent = 4; // Attribute of box boundary
static constexpr Attribute dGround = 5; // Attribute of box boundary

static const Math::Vector x0{{0.5, 0.5}}; // Center of domain

static constexpr Scalar pi = Math::Constants::pi();
static constexpr Scalar gamma_ek = 1.93559;

// static constexpr Scalar gamma_ek =
//   (sectorArea * gammaA1 + (inhomogeinityArea - sectorArea) * gammaA2) / inhomogeinityArea;


static Solver::SparseLU solver;

struct Data
{
  Scalar m;
  Scalar epsilon;
  Scalar error;
};

int main(int, char**)
{
  const size_t N = 100;
  Math::Vector ms = Math::Vector::LinSpaced(N, 0, 0.5 * 1. / hmax);
  Math::Vector es = Math::Vector::LinSpaced(N, hmax, 0.2);
  es.reverseInPlace();
  ms.reverseInPlace();
  std::vector<Data> grid;
  grid.reserve(ms.size() * es.size());

  // Define mesh
  Mesh mesh;
  mesh.load("Q1.medit.mesh", IO::FileFormat::MEDIT);

  // Define finite element spaces
  P1 vh(mesh);
  P1 gh(mesh, mesh.getSpaceDimension());

  std::stringstream filename;
  filename << "Helmholtz_L2_Grid_"
           << "hmax" << std::setw(4) << hmax << "_"
           << "gammaek" << gamma_ek
           << ".live.csv";

  std::ofstream out(filename.str());
  for (auto m : ms)
  {
    for (auto epsilon : es)
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      Alert::Info() << "---------------------"
                    << Alert::NewLine
                    << "Epsilon: " << epsilon
                    << Alert::NewLine
                    << "m: " << m
                    << Alert::NewLine
                    << "Inhomogeinity constant: " << gamma_ek
                    << Alert::Raise;

      // Define conductivity
      ScalarFunction gamma =
        [&](const Point& p)
        {
          return 2 + sin(2 * pi * m * p.x()) * cos(2 * pi * m * p.y());
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

      Problem poisson(u, v);
      poisson = Integral(gamma * Grad(u), Grad(v))
              + DirichletBC(u, phi).on(dCurrent)
              + DirichletBC(u, ScalarFunction(0)).on(dGround);

      Problem perturbed(u, v);
      perturbed = Integral(gamma_e * Grad(u), Grad(v))
                + DirichletBC(u, phi).on(dCurrent)
                + DirichletBC(u, ScalarFunction(0)).on(dGround);

      // Define gradient grid function

      // Solve the background problem
      Alert::Info() << "Solving background equation." << Alert::Raise;
      poisson.solve(solver);
      const auto u0 = u.getSolution();

      // Solve the perturbed problem
      Alert::Info() << "Solving perturbed equation." << Alert::Raise;

      perturbed.solve(solver);
      const auto u_e = u.getSolution();

      GridFunction diff(vh);
      diff = Pow(u0 - u_e, 2);// - Pow(Frobenius(g0 - g_e), 2);
      diff.setWeights();
      const Scalar err = sqrt(Integral(diff).compute());

      Alert::Info() << "Error: " << err << Alert::Raise;
      grid.push_back({ m, epsilon, err });

      out << m << "," << epsilon << "," << err << '\n';
      out.flush();

      auto t1 = std::chrono::high_resolution_clock::now();

      auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
      Alert::Info() << "Time: " << delta.count() << "ms" << Alert::Raise;
    }
  }

  out.close();
  return 0;
}




