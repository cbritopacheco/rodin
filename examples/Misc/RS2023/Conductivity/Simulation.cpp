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

static constexpr Real hmax = 0.005; // Maximal size of a triangle's edge
static constexpr Real hmin = 0.01 * hmax;

static constexpr Attribute dQ = 2; // Attribute of box boundary
static constexpr Attribute dCurrent = 4; // Attribute of box boundary
static constexpr Attribute dGround = 5; // Attribute of box boundary

static const Math::Vector<Real> x0{{0.5, 0.5}}; // Center of domain

static constexpr Real pi = Math::Constants::pi();
static constexpr Real gamma_ek = 1.93559;

struct Data
{
  Real m;
  Real epsilon;
  Real error;
};

int main(int, char**)
{
  const size_t N = 100;
  Math::Vector<Real> ms = Math::Vector<Real>::LinSpaced(N, 0, 0.5 * 1. / hmax);
  Math::Vector<Real> es = Math::Vector<Real>::LinSpaced(N, hmax, 0.2);
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
  filename << "hmax" << std::setw(4) << hmax << "_"
           << "gammaek" << gamma_ek;
  mesh.save("Q_" + filename.str() + ".mesh");
  std::ofstream out("L2_Grid_" + filename.str() + ".live.csv");

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
      RealFunction gamma =
        [&](const Point& p)
        {
          return 2 + sin(2 * pi * m * p.x()) * cos(2 * pi * m * p.y());
        };

      RealFunction gamma_e =
        [&](const Point& p)
        {
          const Real r = (p.getCoordinates() - x0).norm();
          if (r > epsilon)
            return gamma(p);
          else
            return gamma_ek;
        };

      // GridFunction conductivity(vh);
      // conductivity = gamma;
      // conductivity.save("Conductivity.gf");

      // conductivity = gamma_e;
      // conductivity.save("Conductivity_E.gf");

      // Define variational problems
      TrialFunction u(vh);
      TestFunction  v(vh);
      RealFunction phi = 1;

      const auto f =
        [&](const Real& x)
        {
          if (x > 0)
            return std::exp(-1.0 / x);
          else
            return 0.0;
        };

      // Define source
      const auto g =
        [&](const Real& x)
        {
          return f(x) / (f(x) + f(1 - x));
        };

      RealFunction h = 1;
      //   [&](const Geometry::Point& p)
      //   {
      //     const Real r = (p.getCoordinates() - x0).norm();
      //     return g((r - R0) / (R1 - R0));
      //   };

      Problem poisson(u, v);
      poisson = Integral(gamma * Grad(u), Grad(v))
              + DirichletBC(u, phi).on(dCurrent)
              + DirichletBC(u, RealFunction(0)).on(dGround);

      Problem perturbed(u, v);
      perturbed = Integral(gamma_e * Grad(u), Grad(v))
                + DirichletBC(u, phi).on(dCurrent)
                + DirichletBC(u, RealFunction(0)).on(dGround);

      // Define gradient grid function

      // Solve the background problem
      Alert::Info() << "Solving background equation." << Alert::Raise;
      Solver::SparseLU(poisson).solve();
      const auto u0 = u.getSolution();

      // u0.save("Background.gf");

      // Alert::Info() << "Computing its gradient." << Alert::Raise;
      // GridFunction g0(gh);
      // g0 = Grad(u0);
      // // g0.save("BackgroundGradient.gf");

      // Solve the perturbed problem
      Alert::Info() << "Solving perturbed equation." << Alert::Raise;

      Solver::SparseLU(perturbed).solve();
      const auto u_e = u.getSolution();
      // GridFunction g_e(gh);
      // // u_e.save("Perturbed.gf");

      // Alert::Info() << "Computing its gradient." << Alert::Raise;
      // g_e = Grad(u_e);
      // gradient.setWeights();
      // gradient.save("PerturbedGradient.gf");

      GridFunction diff(vh);
      diff = Pow(u0 - u_e, 2);// - Pow(Frobenius(g0 - g_e), 2);
      diff.setWeights();
      const Real err = sqrt(Integral(diff).compute());

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



