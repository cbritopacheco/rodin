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
#include <Rodin/Threads/ThreadPool.h>
#include <Rodin/Solver.h>

#include "Tools.h"

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using namespace Rodin::Examples::BoundaryOptimization;

// Parameters
const Attribute PML = 6;
const Attribute SoundHard = 7;
const Attribute SoundSoft = 8;
const Attribute dBox = 16;
const Attribute dSoundSoft = 11;
const size_t maxIt = 2000;
const Real epsilon = 1e-6;
const Real ell = 1e-12;
const Real tgv = 1e+12;
const Real alpha = 4;
const Real waveLength = 50;
const Real waveNumber = 2 * Math::Constants::pi() / waveLength;

const VectorFunction xi = { 0, 0, 1 };
const Real impedance = 1.0 / 10;
const Real R = 50;

using RealFES = P1<Real>;
using ComplexFES = P1<Complex>;
using VectorFES = P1<Math::Vector<Real>>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

template <class ProblemType>
void solve(ProblemType& pb)
{
  Solver::UMFPack solver(pb);
  solver.solve();
  solver.printStatus();
  // Solver::GMRES solver(pb);
  // solver.setMaxIterations(2000);
  // solver.setTolerance(1e-12);
  // solver.solve();
}

int main(int, char**)
{
  Eigen::initParallel();
  Eigen::setNbThreads(8);
  Threads::getGlobalThreadPool().reset(8);
  std::cout << Eigen::nbThreads() << std::endl;

  Alert::Info() << "Loading mesh..." << Alert::Raise;
  MMG::Mesh mesh;
  mesh.load("AcousticCloaking.medit.mesh", IO::FileFormat::MEDIT);

  Real hmax = waveLength / 2;
  Real hmin = waveLength / 32;
  Real hausd = waveLength / 64;
  Real hgrad = 1.2;

  Alert::Info() << "Initializing sound-soft region..." << Alert::Raise;
  {
    P1 vh(mesh);

    GridFunction dist(vh);
    dist = [&](const Point& p)
    {
      Math::SpatialVector<Real> c(3);
      c << -45.189708, 0, -11.330562;
      return (p - c).norm() - 2;
    };

    mesh.save("Dist.mesh", IO::FileFormat::MEDIT);
    dist.save("Dist.sol", IO::FileFormat::MEDIT);

    mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                      .split(SoundSoft, { SoundSoft, SoundHard })
                                      .split(SoundHard, { SoundSoft, SoundHard })
                                      .noSplit(PML)
                                      .noSplit(dBox)
                                      .setHMax(hmax)
                                      .setHMin(hmin)
                                      .setGradation(hgrad)
                                      .setHausdorff(hausd)
                                      .surface()
                                      .discretize(dist);
  }

  mesh.save("Omega0.mfem.mesh", IO::FileFormat::MFEM);
  mesh.save("Omega0.mesh", IO::FileFormat::MEDIT);

  std::ofstream fObj("obj.txt");
  size_t i = 0;
  size_t regionCount;
  while (i < maxIt)
  {
    //bool topologicalStep = i < 5 || (i < 200  && i % 20 == 0);
    bool topologicalStep = true;
    bool geometricStep = !topologicalStep;

    const Real k = 0.5 * (hmax + hmin);
    const Real dt = 2 * hmin;
    const Real radius = 2 * hmin;
    Alert::Info() << "Iteration: " << i                         << Alert::NewLine
                  << "HMax:      " << Alert::Notation(hmax)     << Alert::NewLine
                  << "HMin:      " << Alert::Notation(hmin)     << Alert::NewLine
                  << "Hausdorff: " << Alert::Notation(hausd)    << Alert::NewLine
                  << "HGrad:     " << Alert::Notation(hgrad)    << Alert::NewLine
                  << "ell:     " << Alert::Notation(ell)    << Alert::NewLine
                  << "dt:        " << Alert::Notation(dt)       << Alert::Raise;

    Alert::Info() << "Optimizing the domain..." << Alert::Raise;
    MMG::Optimizer().setHMax(hmax)
                    .setHMin(hmin)
                    .setHausdorff(hausd)
                    .setGradation(hgrad)
                    .setAngleDetection(false)
                    .optimize(mesh);

    Alert::Info() << "Computing required connectivity..." << Alert::Raise;
    mesh.getConnectivity().compute(2, 3); // Computes boundary
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(1, 0);

    Alert::Info() << "RMC..." << Alert::Raise;
    // regionCount = rmc(mesh, { SoundSoft }, SoundHard);

    Alert::Info() << "Skinning mesh..." << Alert::Raise;
    auto dOmega = mesh.skin();

    Alert::Info() << "Tracing mesh..." << Alert::Raise;
    dOmega.trace({ {{ SoundHard, SoundSoft }, dSoundSoft } });
    dOmega.save("dK.mesh", IO::FileFormat::MEDIT);

    ComplexFunction wave =
      [&](const Point& p)
      {
        const Real d = waveNumber * p.getCoordinates().dot(xi(p));
        return Complex(cos(d), sin(d));
      };

    BoundaryNormal normal(mesh);
    ComplexFunction dnWave = waveNumber * Complex(0, 1) * Dot(xi, normal) * wave;

    Alert::Info() << "Building finite element spaces..." << Alert::Raise;

    ComplexFES cfes(mesh);
    RealFES rfes(mesh);
    RealFES drfes(dOmega);
    VectorFES drvfes(dOmega, mesh.getSpaceDimension());

    GridFunction phi(rfes);
    phi.getFiniteElementSpace().getMesh().save("Phi.mesh");

    phi = Re(wave);
    phi.save("PhiRe.gf");

    phi = Im(wave);
    phi.save("PhiIm.gf");

    Alert::Info() << "Distancing absorbing domain..." << Alert::Raise;
    auto dist = MMG::Distancer(drfes).setInteriorDomain(SoundSoft)
                                     .distance(dOmega);

    dist.save("Dist.gf");
    dist.getFiniteElementSpace().getMesh().save("Dist.mesh");

    Alert::Info() << "Assembling state equation..." << Alert::Raise;
    TrialFunction u(cfes);
    TestFunction  v(cfes);
    Problem state(u, v);
    state = Integral(Grad(u), Grad(v))
          - waveNumber * waveNumber * Integral(u, v)
          + BoundaryIntegral(dnWave, v).over(SoundSoft, SoundHard)
          + Complex(0, waveNumber / impedance) * BoundaryIntegral(u, v).over(SoundSoft)
          + Complex(0, waveNumber / impedance) * BoundaryIntegral(wave, v).over(SoundSoft)
          + Complex(1 / R, waveNumber) * FaceIntegral(u, v).over(PML)
          ;

    state.assemble();

    Alert::Info() << "Solving..." << Alert::Raise;
    solve(state);

    GridFunction total(rfes);
    total = [&](const Point& x) { return std::norm(wave(x) + u.getSolution()(x)); };
    total.save("Total.gf");
    rfes.getMesh().save("Total.mesh");

    GridFunction pressure(rfes);
    Math::SpatialVector<Real> c({{-50, 0, 0}});
    pressure = [&](const Point& p)
    {
      if ((p - c).norm() < 75)
        return std::norm(u.getSolution()(p));
      else return 0.0;
    };
    pressure.save("Scattered.gf");
    rfes.getMesh().save("Scattered.mesh");

    Alert::Info() << "Computing objective..." << Alert::Raise;
    pressure.setWeights();
    const Real J = 0.5 * Integral(pressure).compute();
    const Real area = mesh.getArea(SoundSoft);
    const Real objective = J + ell * area;

    Alert::Info() << "Objective: " << Alert::Notation(objective) << Alert::NewLine
                  << "J: " << Alert::Notation(J) << Alert::NewLine
                  << "Area: " << Alert::Notation(area) << Alert::NewLine
                  << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

    Alert::Info() << "Assembling adjoint equation..." << Alert::Raise;
    TrialFunction p(cfes);
    TestFunction  q(cfes);
    Problem adjoint(p, q);
    adjoint = Integral(Grad(p), Grad(q))
            - waveNumber * waveNumber * Integral(p, q).over(SoundSoft)
            + Integral(u.getSolution(), q)
            + Complex(0, -waveNumber / impedance) * BoundaryIntegral(p, q).over(SoundSoft)
            + Complex(1 / R, -waveNumber) * FaceIntegral(p, q).over(PML)
            ;

    adjoint.assemble();

    Alert::Info() << "Solving..." << Alert::Raise;
    solve(adjoint);

    pressure = [&](const Point& x) { return std::norm(p.getSolution()(x)); };
    pressure.save("Adjoint.gf");
    p.getSolution().getFiniteElementSpace().getMesh().save("Adjoint.mesh");

    if (topologicalStep)
    {
      Alert::Info() << "Topological optimization..." << Alert::Raise;
      Alert::Info() << "Inserting sound-soft region..." << Alert::Raise;
      TrialFunction s(drfes);
      TestFunction  t(drfes);
      Problem topo(s, t);
      topo = alpha * alpha * Integral(Grad(s), Grad(t))
           + Integral(s, t)
           - Integral(
               Im(waveNumber / impedance * Conjugate(u.getSolution() + wave) * p.getSolution()), t).over(SoundHard)
           + tgv * Integral(s, t).over(dBox, SoundSoft);

      topo.assemble();

      solve(topo);

      s.getSolution().save("Topo.gf");
      drfes.getMesh().save("Topo.mesh");

      Alert::Info() << "Computing nucleation locations..." << Alert::Raise;
      auto cs = locations(s.getSolution());

      Alert::Info() << "Nucleating " << Alert::Notation(cs.size()) << " holes..."
                    << Alert::Raise;
      holes(radius, dist, cs);
    }
    else if (geometricStep)
    {
      Alert::Info() << "Geometrical optimization..." << Alert::Raise;

      Alert::Info() << "Computing conormal..." << Alert::Raise;
      GridFunction conormal(drvfes);
      conormal.project(Grad(dist), { SoundSoft, SoundHard });
      conormal.stableNormalize();

      conormal.getFiniteElementSpace().getMesh().save("Conormal.mesh");
      conormal.save("Conormal.gf");

      Alert::Info() << "Computing shape gradient..." << Alert::Raise;

      TrialFunction theta(drvfes);
      TestFunction  w(drvfes);
      Problem hilbert(theta, w);
      hilbert = alpha * alpha * Integral(Jacobian(theta), Jacobian(w))
              + Integral(theta, w)
              - FaceIntegral(
                  Im(waveNumber / impedance * Conjugate(u.getSolution() + wave) * p.getSolution()),
                  Dot(conormal, w)).over(dSoundSoft)
              ;

      hilbert.assemble();

      solve(hilbert);

      GridFunction norm(drfes);
      norm = Frobenius(theta.getSolution());
      theta.getSolution() /= norm.max();

      theta.getFiniteElementSpace().getMesh().save("Theta.mesh");
      theta.getSolution().save("Theta.gf");

      Alert::Info() << "Advecting support..." << Alert::Raise;
      MMG::Advect(dist, theta.getSolution()).step(dt);
    }

    dist.save("Dist.gf");
    dist.getFiniteElementSpace().getMesh().save("Dist.mesh");

    Alert::Info() << "Meshing the support..." << Alert::Raise;
    RealFES workaroundfes(mesh);
    GridFunction workaround(workaroundfes);
    workaround.projectOnBoundary(dist);
    mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                      .split(SoundSoft, { SoundSoft, SoundHard })
                                      .split(SoundHard, { SoundSoft, SoundHard })
                                      .noSplit(PML)
                                      .noSplit(dBox)
                                      .setHMax(hmax)
                                      .setHMin(hmin)
                                      .setGradation(hgrad)
                                      .setHausdorff(hausd)
                                      .surface()
                                      .discretize(workaround);

    Alert::Info() << "Saving files..." << Alert::Raise;
    mesh.save("Omega.mesh", IO::FileFormat::MEDIT);
    dOmega.save("out/dOmega." + std::to_string(i) +  ".mesh", IO::FileFormat::MEDIT);
    dOmega.save("out/dOmega.mfem." + std::to_string(i) +  ".mesh", IO::FileFormat::MFEM);

    Alert::Success() << "Completed Iteration: " << i << '\n' << Alert::Raise;
    i++;
  }
  return 0;
}

