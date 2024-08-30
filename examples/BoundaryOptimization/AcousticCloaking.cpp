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
static constexpr Attribute Aircraft = 3;
static constexpr Attribute Air = 4;
static constexpr Attribute ROI = 5;

static constexpr Attribute PML = 6;
static constexpr Attribute SoundHard = 7;
static constexpr Attribute SoundSoft = 8;
static constexpr Attribute Window = 9;

static constexpr Attribute dBox = 10;
static constexpr Attribute dSoundSoft = 11;

static constexpr size_t maxIt = 2000;

static constexpr Real epsilon = 1e-6;
static constexpr Real ell = 0;
static constexpr Real tgv = std::numeric_limits<float>::max();
static constexpr Real alpha = 2;
static constexpr Real angle = 0;
static constexpr Real waveLength = 50;
static constexpr Real waveNumber = 2 * Math::Constants::pi() / waveLength;

using RealFES = P1<Real>;
using ComplexFES = P1<Complex>;
using VectorFES = P1<Math::Vector<Real>>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;


int main(int, char**)
{
  Eigen::initParallel();
  Eigen::setNbThreads(8);
  Threads::getGlobalThreadPool().reset(8);
  std::cout << Eigen::nbThreads() << std::endl;

  MMG::Mesh mesh;
  mesh.load("../resources/mmg/AircraftPML.medit.mesh", IO::FileFormat::MEDIT);

  Real hmax = waveLength / 2;
  Real hmin = waveLength / 8;
  Real hausd = 0.5 * hmin;
  Real hgrad = 1.2;

  std::ofstream fObj("obj.txt");
  size_t i = 0;
  size_t regionCount;
  while (i < maxIt)
  {
    // bool topologicalStep = (i < 10);// || (i < 200  && i % 20 == 0);
    bool topologicalStep = true;
    bool geometricStep = !topologicalStep;

    const Real k = 0.5 * (hmax + hmin);
    const Real dt = alpha * hmin;
    const Real radius = 0.5 * k;
    Alert::Info() << "Iteration: " << i                         << Alert::NewLine
                  << "HMax:      " << Alert::Notation(hmax)     << Alert::NewLine
                  << "HMin:      " << Alert::Notation(hmin)     << Alert::NewLine
                  << "Hausdorff: " << Alert::Notation(hausd)    << Alert::NewLine
                  << "HGrad:     " << Alert::Notation(hgrad)    << Alert::NewLine
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
    regionCount = rmc(mesh, { SoundSoft }, SoundHard);

    Alert::Info() << "Found " << Alert::Notation(regionCount) << " regions."
      << Alert::Raise;

    Alert::Info() << "Trimming mesh..." << Alert::Raise;
    auto perturbed = mesh.trim(Aircraft);
    perturbed.save("Perturbed.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "Skinning mesh..." << Alert::Raise;
    auto dPerturbed = perturbed.skin();

    Alert::Info() << "Tracing mesh..." << Alert::Raise;
    dPerturbed.trace({ {{ SoundHard, SoundSoft }, dSoundSoft } });
    dPerturbed.save("dK.mesh", IO::FileFormat::MEDIT);

    // Parameters
    VectorFunction xi = { 0, 0, 1 };
    Real impedance = 1;
    Real R = 50;

    ComplexFunction wave =
      [&](const Point& p)
      {
        const Real d = waveNumber * p.getCoordinates().dot(xi(p));
        return Complex(cos(d), -sin(d));
      };

    BoundaryNormal normal(perturbed);

    ComplexFunction dnWave = waveNumber * Complex(0, 1) * Dot(xi, normal) * wave;

    // Bump function
    auto h = [](Real r)
    {
      if (r < -1.0)
        return 1.0;
      else if (r > 1.0)
        return 0.0;
      else
        return 1.0 - 1.0 / (1.0 + std::exp(4 * r / (r * r - 1.0)));
    };

    Alert::Info() << "Building finite element spaces..." << Alert::Raise;

    ComplexFES cfes(perturbed);

    RealFES rfes(perturbed);
    RealFES drfes(dPerturbed);

    GridFunction phi(rfes);
    phi.getFiniteElementSpace().getMesh().save("Phi.mesh");
    phi = Re(wave);
    phi.save("PhiRe.gf");
    phi = Im(wave);
    phi.save("PhiIm.gf");

    GridFunction dnRe(rfes);
    dnRe.projectOnBoundary(Re(dnWave));
    dnRe.save("dnWaveRe.gf");

    GridFunction dnIm(rfes);
    dnIm.projectOnBoundary(Im(dnWave));
    dnIm.save("dnWaveIm.gf");

    Alert::Info() << "Distancing absorbing domain..." << Alert::Raise;
    auto dist = MMG::Distancer(drfes).setInteriorDomain(SoundSoft)
                                     .distance(dPerturbed);

    dist.save("Dist.gf");
    dist.getFiniteElementSpace().getMesh().save("Dist.mesh");

    RealFunction he =
      [&](const Geometry::Point& p) { return h(dist(p) / epsilon) / epsilon; };

    Alert::Info() << "Assembling state equation..." << Alert::Raise;
    TrialFunction u(cfes);
    TestFunction  v(cfes);
    Problem state(u, v);
    state = Integral(Grad(u), Grad(v))
          - waveNumber * waveNumber * Integral(u, v)
          - Complex(-1 / R, waveNumber) * FaceIntegral(u, v).over(PML)
          - BoundaryIntegral(dnWave, v).over(SoundSoft, SoundHard, Window)
          - Complex(0, waveNumber / impedance) * BoundaryIntegral(u, v).over(SoundSoft)
          - Complex(0, waveNumber / impedance) * BoundaryIntegral(wave, v).over(SoundSoft);

    state.assemble();

    Alert::Info() << "Solving..." << Alert::Raise;
    Solver::UMFPack solver(state);
    solver.solve();
    solver.printStatus();

    GridFunction pressure(rfes);
    pressure.getFiniteElementSpace().getMesh().save("Pressure.mesh");
    pressure = [&](const Point& p) { return std::abs(u.getSolution()(p)); };
    pressure.save("Scattered.gf");

    GridFunction uRe(rfes);
    uRe = Re(u.getSolution());
    uRe.save("uRe.gf");

    GridFunction uIm(rfes);
    uIm = Im(u.getSolution());
    uIm.save("uIm.gf");

    GridFunction total(cfes);
    total = wave + u.getSolution();

    pressure = [&](const Point& p) { return std::abs(total(p)); };
    pressure.save("Total.gf");

    std::exit(1);

    // Alert::Info() << "Solving state equation..." << Alert::Raise;
    // TrialFunction u(spfes);
    // TestFunction  v(spfes);
    // Problem state(u, v);
    // state = Integral(Grad(u), Grad(v))
    //       - Integral(waveNumber * waveNumber * gamma * u, v)
    //       + FaceIntegral(he * u, v).over({ Absorbing, Reflecting, Window })
    //       + DirichletBC(u, phi).on(PlaneWave, Dirichlet);
    //       //+ DirichletBC(u, Zero()).on(Dirichlet);
    // Solver::CG(state).solve();

    // u.getSolution().save("State.gf");
    // u.getSolution().getFiniteElementSpace().getMesh().save("State.mesh");

    // GridFunction diff(spfes);
    // diff = 0.5 * Pow(u.getSolution() - u0.getSolution(), 2);
    // diff.save("Diff.gf");
    // diff.getFiniteElementSpace().getMesh().save("Diff.mesh");

    // Alert::Info() << "Solving adjoint equation..." << Alert::Raise;
    // TrialFunction p(spfes);
    // TestFunction  q(spfes);

    // Problem adjoint(p, q);
    // adjoint = Integral(Grad(p), Grad(q))
    //         - Integral(waveNumber * waveNumber * gamma * p, q)
    //         + Integral(u.getSolution() - u0.getSolution(), q)
    //         + FaceIntegral(he * p, q).over({ Absorbing, Reflecting, Window })
    //         + DirichletBC(u, Zero()).on({ Dirichlet, PlaneWave });
    //         ;
    // Solver::CG(adjoint).solve();

    // p.getSolution().save("Adjoint.gf");
    // p.getSolution().getFiniteElementSpace().getMesh().save("Adjoint.mesh");

    // Alert::Info() << "Computing objective..." << Alert::Raise;
    // diff.setWeights();
    // if (i > 0)
    // {
    //   oldAugmented = augmented;
    //   oldObjective = objective;
    //   oldConstraint = constraint;
    // }

    // const Real J = Integral(diff).compute();
    // const Real area = mesh.getArea(Absorbing);
    // const Real perimeter = dPerturbed.getMeasure(1, dAbsorbing);
    // objective = J;
    // constraint = (area / targetArea - 1);
    // augmented = objective;// + ellA * constraint + 0.5 * bA * constraint * constraint;

    // Alert::Info() << "Objective: " << Alert::Notation(objective) << Alert::NewLine
    //               << "Augmented: " << Alert::Notation(augmented) << Alert::NewLine
    //               << "Absorbing Area: " << Alert::Notation(area) << Alert::NewLine
    //               << "Constraint: " << Alert::Notation(constraint) << Alert::NewLine
    //               << Alert::Raise;
    // fObj << augmented << "\n";
    // fObj.flush();

    // if (topologicalStep)
    // {
    //   Alert::Info() << "Topological optimization..." << Alert::Raise;
    //   Alert::Info() << "Inserting supporting region..." << Alert::Raise;
    //   TrialFunction s(dsfes);
    //   TestFunction  t(dsfes);
    //   Problem topo(s, t);
    //   topo = alpha * alpha * dt * dt * Integral(Grad(s), Grad(t))
    //        + Integral(s, t)
    //        + Integral(u.getSolution() * p.getSolution(), t)
    //        + tgv * Integral(s, t).over(Absorbing, Window, Dirichlet, PlaneWave);
    //   Solver::CG(topo).solve();

    //   s.getSolution().save("Topo.gf");
    //   dsfes.getMesh().save("Topo.mesh");

    //   GridFunction raw(dsfes);
    //   raw = -u.getSolution() * p.getSolution();
    //   raw.save("Raw.gf");
    //   dsfes.getMesh().save("Raw.mesh");

    //   Alert::Info() << "Computing nucleation locations..." << Alert::Raise;
    //   const Real tc = s.getSolution().max();
    //   std::vector<Point> cs;
    //   for (auto it = dPerturbed.getVertex(); !it.end(); ++it)
    //   {
    //     const Point p(*it, it->getTransformation(),
    //         Polytope::getVertex(0, Polytope::Type::Point), it->getCoordinates());
    //     const Real tp = s.getSolution()(p);
    //     if (tp > 1e-6 && abs(tp / tc) > (1 - 1e-12))
    //     {
    //       cs.emplace_back(std::move(p));
    //       break;
    //     }
    //   }

    //   if (cs.size() > 0)
    //   {

    //     Alert::Info() << "Nucleating " << Alert::Notation(cs.size()) << " holes..."
    //                   << Alert::Raise;
    //     holes(radius, dist, cs);
    //   }
    //   else
    //   {
    //     Alert::Info() << "No holes to nucleate !" << Alert::Raise;
    //   }
    // }
    // else if (geometricStep)
    // {
    //   Alert::Info() << "Geometrical optimization..." << Alert::Raise;

    //   Alert::Info() << "Computing conormal..." << Alert::Raise;

    //   GridFunction conormal(dvfes);
    //   conormal.projectOnFaces(Average(Grad(dist)));
    //   conormal.stableNormalize();

    //   conormal.getFiniteElementSpace().getMesh().save("Conormal.mesh");
    //   conormal.save("Conormal.gf");

    //   Alert::Info() << "Computing shape gradient..." << Alert::Raise;

    //   TrialFunction theta(dvfes);
    //   TestFunction  w(dvfes);
    //   Problem hilbert(theta, w);
    //   hilbert = alpha * alpha * dt * dt * Integral(Jacobian(theta), Jacobian(w))
    //           + Integral(theta, w)
    //           + 1.0 / epsilon * FaceIntegral(
    //               Dot(u.getSolution(), p.getSolution()), Dot(conormal, w)).over(dAbsorbing)
    //           + constraint * ellA * FaceIntegral(conormal, w).over(dAbsorbing)
    //           + 0.5 * 2 * bA * constraint * FaceIntegral(conormal, w).over(dAbsorbing)
    //           + tgv * Integral(theta, w).over(Dirichlet, PlaneWave, Window)
    //           ;
    //   Solver::CG(hilbert).solve();

    //   GridFunction norm(dsfes);
    //   norm = Frobenius(theta.getSolution());
    //   theta.getSolution() /= norm.max();

    //   // dPerturbed.save("Theta.mesh");
    //   // theta.getSolution().save("Theta.gf");

    //   dPerturbed.save("Theta.mesh", IO::FileFormat::MEDIT);
    //   theta.getSolution().save("Theta.sol", IO::FileFormat::MEDIT);

    //   Alert::Info() << "Advecting support..." << Alert::Raise;
    //   MMG::Advect(dist, theta.getSolution()).step(dt);
    // }

    // dist.save("NewDist.gf");
    // dist.getFiniteElementSpace().getMesh().save("NewDist.mesh");

    // Alert::Info() << "Meshing the support..." << Alert::Raise;
    // try
    // {
    //   GridFunction workaround(sfes);
    //   workaround.projectOnFaces(dist, { Reflecting, Window, Absorbing, Dirichlet, PlaneWave });
    //   workaround.save("Workaround.gf");
    //   workaround.getFiniteElementSpace().getMesh().save("Workaround.mesh");
    //   mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
    //                                     .split(Absorbing, { Absorbing, Reflecting })
    //                                     .split(Reflecting, { Absorbing, Reflecting })
    //                                     .noSplit(Window)
    //                                     .noSplit(Dirichlet)
    //                                     .noSplit(PlaneWave)
    //                                     .setHMax(hmax)
    //                                     .setHMin(hmin)
    //                                     .setGradation(hgrad)
    //                                     .setHausdorff(hausd)
    //                                     .surface()
    //                                     .discretize(workaround);
    //   hmax = 1.1 * hmax > resolution * waveLength ? resolution * waveLength : hmax * 1.1;
    // }
    // catch (Alert::Exception& e)
    // {
    //   Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
    //   hmax = 0.9 * hmax < 0.01 ? 0.01 : hmax * 0.9;
    //   continue;
    // }

    // Alert::Info() << "Saving files..." << Alert::Raise;
    // mesh.save("Omega.mesh", IO::FileFormat::MEDIT);
    // dPerturbed.save("out/dOmega." + std::to_string(i) +  ".mesh", IO::FileFormat::MEDIT);
    // dPerturbed.save("out/dOmega.mfem." + std::to_string(i) +  ".mesh", IO::FileFormat::MFEM);

    // Alert::Success() << "Completed Iteration: " << i << '\n' << Alert::Raise;

    // ellA += bA * constraint;
    // if (i > 0)
    //   bA *= 1 + alpha * alpha * abs(oldAugmented - augmented);
    // bA = fmin(bTarget, bA);

    // i++;
  }
  return 0;
}

