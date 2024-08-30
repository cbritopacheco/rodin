/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Math.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::External;

const Real E = 100;
const Real nu = 0.48;
const Real lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2 * nu));
const Real mu = E / (1.0 + nu);

inline
void K(Math::Matrix<Real>& res, const Point& x, const Point& y)
{
  const auto norm = (x - y).norm();
  res.resize(3, 3);
  const Real L00 = (1 - nu) / (2 * M_PI * mu * norm) + nu * (x(0) - y(0)) * (x(0) - y(0)) / (2 * M_PI * mu * norm * norm * norm);
  const Real L01 = nu * (x(0) - y(0)) * (x(1) - y(1)) / (2 * M_PI * mu * norm * norm * norm);

  const Real L10 = nu * (x(1) - y(1)) * (x(0) - y(0)) / (2 * M_PI * mu * norm * norm * norm);
  const Real L11 = (1 - nu) / (2 * M_PI * mu * norm) + nu * (x(1) - y(1)) * (x(1) - y(1)) / (2 * M_PI * mu * norm * norm * norm);

  const Real L20 = -(1 - 2 * nu) * (x(0) - y(0)) / (4 * M_PI * mu * norm * norm);
  const Real L21 = -(1 - 2 * nu) * (x(1) - y(1)) / (4 * M_PI * mu * norm * norm);
  const Real L22 = (1 - nu) / (2 * M_PI * mu * norm);

  const Real L02 = -L20;
  const Real L12 = -L21;

  res << L00, L01, L02,
         L10, L11, L12,
         L20, L21, L22;
}

int main(int, char**)
{
  Eigen::initParallel();
  Eigen::setNbThreads(8);
  Threads::getGlobalThreadPool().reset(8);

  std::cout << "lambda: " << lambda << std::endl;
  std::cout << "mu: " << mu << std::endl;

  MMG::Mesh mesh;
  mesh.load("D1.o.mesh", IO::FileFormat::MEDIT);
  mesh.save("D1.mfem.mesh");

  VectorFunction e1{1, 0, 0};
  VectorFunction e2{0, 1, 0};
  VectorFunction e3{0, 0, 1};

  std::ofstream data("e1.txt");
  data << "eps,hmax,err\n";
  Real eps0 = 1e-5, eps1 = 0.1;
  Real hmax0 = 0.5, hmax1 = 1;
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

      P1 vfes(mesh, 3);
      P1 sfes(mesh);
      TrialFunction u(vfes);
      TestFunction  v(vfes);

      DenseProblem eq(u, v);
      eq = eps * LinearElasticityIntegral(u, v)(lambda, mu)
         + Integral(Potential(K, u), v)
         - Integral(e1, v);

      Alert::Info() << "Assembling." << Alert::Raise;
      eq.assemble();

      Alert::Info() << "Solving." << Alert::Raise;
      Solver::CG(eq).solve();

      mesh.save("u.mesh");
      u.getSolution().save("u.gf");

      GridFunction frob(sfes);
      frob = Frobenius(u.getSolution());
      u.getSolution() /= 10 * frob.max();

      mesh.save("scaled.mesh");
      u.getSolution().save("scaled.gf");
      std::exit(1);

      u.getSolution().setWeights();

      Alert::Info() << "Getting data." << Alert::Raise;
      GridFunction phi(sfes);
      phi = Pow(Frobenius(e1 - Potential(K, u.getSolution())), 2);
      phi.setWeights();

      Real err = Integral(phi).compute();

      data << eps << "," << hmax << "," << err << '\n';
      data.flush();
    }
  }

  return 0;
}

