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
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::External;

inline
Real K(const Point& x, const Point& y)
{
  return 1. / (4 * M_PI * ((x - y).norm()));
}

inline
Real exact(const Point& x)
{
  return 4. / (M_PI * std::sqrt(abs(1 - x.squaredNorm())));
}

int main(int, char**)
{
  Eigen::initParallel();
  Eigen::setNbThreads(8);
  Threads::getGlobalThreadPool().reset(8);
  MMG::Mesh mesh;
  mesh.load("D1.mesh");

  mesh.getConnectivity().compute(1, 0);

  Real max = -1;
  Real length = 0;
  for (auto it = mesh.getFace(); it; ++it)
  {
    length += it->getMeasure();
    if (it->getMeasure() > max)
      max = it->getMeasure();
  }
  Alert::Info() << "Length: " << length / mesh.getFaceCount() << Alert::NewLine << Alert::Raise;
  std::cout << max << std::endl;



  std::exit(1);


  std::ofstream data("data.txt");
  data << "eps,hmax,avg,err\n";
  Real eps0 = 1e-5, eps1 = 0.1;
  Real hmax0 = 0.1, hmax1 = 1;
  size_t n = 50;
  for (double eps = eps0; eps < eps1; eps += (eps1 - eps0) / n)
  {
    for (double hmax = hmax0; hmax < hmax1; hmax += (hmax1 - hmax0) / n)
    {
      Alert::Info() << "eps: " << eps << Alert::NewLine
                    << "hmax: " << hmax << Alert::Raise;

      Alert::Info() << "Optimizing." << Alert::Raise;

      double hmin = hmax / 10;
      MMG::Optimizer().setHMax(hmax)
                      .setHMin(hmin)
                      .setHausdorff(hmax / 20)
                      .optimize(mesh);

      Alert::Info() << "Connectivity." << Alert::Raise;
      mesh.getConnectivity().compute(1, 2);

      P1 fes(mesh);
      TrialFunction u(fes);
      TestFunction  v(fes);

      Problem eq(u, v);
      eq = Integral(0.0001 * Grad(u), Grad(v))
         + Integral(Potential(K, u), v)
         - Integral(v)
         ;
      Alert::Info() << "Assembling." << Alert::Raise;
      eq.assemble();

      Alert::Info() << "Solving." << Alert::Raise;
      Solver::CG(eq).solve();

      mesh.save("u.mesh");
      u.getSolution().save("u.gf");


      GridFunction trunc(fes);
      trunc = [&](const Point& p)
      {
        if (p.norm() < 0.8)
          return u.getSolution()(p);
        else
          return 0.0;
      };

      double max = trunc.max();

      trunc = [&](const Point& p)
      {
        if (p.norm() < 0.8)
          return u.getSolution()(p);
        else
          return max;
      };

      trunc.save("trunc.gf");
      mesh.save("trunc.mesh");

      GridFunction uEx(fes);
      uEx = [&](const Point& p)
      {
        if (p.norm() < 0.8)
          return exact(p);
        else
          return 0.0;
      };

      uEx = trunc.max();

      uEx = [&](const Point& p)
      {
        if (p.norm() < 0.8)
          return exact(p);
        else
          return max;
      };
      uEx.save("uEx.gf");
      mesh.save("uEx.mesh");
      std::exit(1);

      u.getSolution().setWeights();

      Alert::Info() << "Getting data." << Alert::Raise;

      GridFunction one(fes);
      one = RealFunction(1);
      one.setWeights();

      GridFunction phi(fes);
      phi = Potential(K, u.getSolution());
      phi.setWeights();

      GridFunction resError(fes);
      resError = Pow(one - phi, 2);
      resError.setWeights();


      GridFunction err(fes);
      err = [&](const Point& p)
      {
        if (p.norm() < 0.9)
          return pow(u.getSolution()(p) - exact(p), 2);
        else
          return 0.0;
      };
      err.setWeights();

      const Real Aerr = pow(Integral(u.getSolution()) - 8, 2);
      const Real Rerr = Integral(resError);
      const Real Eerr = Integral(err);

      data << eps << "," << hmax << "," << Aerr << "," << Rerr << "," << Eerr << '\n';
      data.flush();
    }
  }

  return 0;
}
