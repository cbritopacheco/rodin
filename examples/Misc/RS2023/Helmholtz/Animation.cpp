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
static constexpr Scalar hmin = 0.1 * hmax;

static constexpr Attribute dQ = 2; // Attribute of box boundary
static constexpr Attribute dCurrent = 4; // Attribute of box boundary
static constexpr Attribute dGround = 5; // Attribute of box boundary

static const Math::Vector x0{{0.5, 0.5}}; // Center of domain

static constexpr Scalar pi = Math::Constants::pi();

static constexpr Scalar R0 = 0.2; // Radius of B_R(x_0)
static constexpr Scalar R1 = R0 + 10 * hmax; // Radius of B_R(x_0)

struct Data
{
  const Scalar m;
  const Scalar epsilon;
  const Scalar waveNumber;
  const Scalar conductivity;
};

template <class T>
std::vector<std::vector<T>> splitVector(const std::vector<T>& vec, int n) {
    std::vector<std::vector<T>> result;
    int size = vec.size();
    int partSize = size / n;
    int remainder = size % n;
    int startIndex = 0;

    for (int i = 0; i < n; i++) {
        int currentPartSize = partSize + (i < remainder ? 1 : 0);
        std::vector<T> currentPart(vec.begin() + startIndex, vec.begin() + startIndex + currentPartSize);
        result.push_back(currentPart);
        startIndex += currentPartSize;
    }

    return result;
}

void run(int i, const std::vector<Data>& grid);

int main(int, char**)
{
  // Define evaluation grid
  Math::Vector m_r = Math::Vector::LinSpaced(50, 0, 0.5 * 1. / hmax);
  Math::Vector epsilon_r = Math::Vector::LinSpaced(1. / hmax, hmax, 0.2);
  Math::Vector waveNumber_r = Math::Vector::LinSpaced(1.0 / hmax, 1, 1.0 / hmax);
  Math::Vector conductivity_r{{ 1e-12, 0.5, 1.0, 2.0, 1e12 }};

  std::vector<Data> grid;
  grid.reserve(m_r.size() * epsilon_r.size() * waveNumber_r.size() * conductivity_r.size());
  for (const Scalar m : m_r)
    for (const Scalar epsilon : epsilon_r)
      for (const Scalar waveNumber : waveNumber_r)
        for (const Scalar g : conductivity_r)
          grid.push_back({ m, epsilon, waveNumber, g });

}

void run(int id, const std::vector<Data>& grid)
{
  // Load mesh
  Mesh mesh;
  mesh.load("Q1.medit.mesh", IO::FileFormat::MEDIT);
  mesh.save("out/Q.mesh");

  Solver::SparseLU solver;

  Alert::Info() << "Grid size: " << grid.size() << Alert::Raise;

  size_t i = 0;
  auto t0 = std::chrono::high_resolution_clock::now();
  for (const auto& data : grid)
  {
      auto t1 = std::chrono::high_resolution_clock::now();
      const auto delta = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0);
      Alert::Info() << i << " / " << grid.size() << " ---- "
                    << (100 * Scalar(i) / Scalar(grid.size())) << "%"
                    << Alert::NewLine
                    << "Total elapsed time: " << delta.count() << "s"
                    << Alert::Raise;

    // Alert::Info() << ">>>>>>>>>>>>>>>>>>>>"
    //               << Alert::NewLine
    //               << i << " / " << grid.size()
    //               << Alert::NewLine
    //               << "Epsilon: " << data.epsilon
    //               << Alert::NewLine
    //               << "m: " << data.m
    //               << Alert::NewLine
    //               << "waveNumber: " << data.waveNumber
    //               << Alert::NewLine
    //               << "Inhomogeinity constant: " << data.conductivity
    //               << Alert::Raise;

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

    ScalarFunction chi_e =
      [&](const Point& p)
      { return Scalar((p.getCoordinates() - x0).norm() > 2 * data.epsilon); };

    GridFunction diff(vh);
    diff = chi_e * Pow(u0 - ue, 2);
    diff.setWeights();
    const Scalar error = Integral(diff);
    Alert::Success() << "L2 Error: " << error
                     << Alert::NewLine
                     << "<<<<<<<<<<<<<<<<<<<<"
                     << Alert::Raise;

    std::stringstream ss;
    ss << std::setfill('0') << std::setw(4) << i;
    u0.save("out/Background_" + ss.str() + ".gf");
    ue.save("out/Perturbed_" + ss.str() + ".gf");

    i++;
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  const auto delta = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0);

  Alert::Success() << "Finished " << id << " in " << delta.count() << "s" << Alert::Raise;
}

