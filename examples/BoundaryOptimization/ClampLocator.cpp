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
#include <Rodin/Threads/ThreadPool.h>

#include "Tools.h"

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using namespace Rodin::Examples::BoundaryOptimization;

// Parameters
static constexpr Geometry::Attribute Gamma = 1; // dn u = 0
static constexpr Geometry::Attribute Clamp = 2; // u = 1
static constexpr Geometry::Attribute Locator = 3; // u = f
static constexpr Geometry::Attribute GammaT = 4; // dn u = g
static constexpr Geometry::Attribute Fixed = 5; // dn u = g

static constexpr Geometry::Attribute dClamp = 2;
static constexpr Geometry::Attribute dLocator = 3;

static constexpr size_t maxIt = 1000;

static constexpr Real epsilon = 1e-6;
static constexpr Real ellClamp = 1e-4;
static constexpr Real ellLocator = 1e-4;
static constexpr Real radius = 0.2;
static constexpr Real tgv = std::numeric_limits<float>::max();
static const Real alpha = 4;

using RealFES = P1<Real>;
using VectorFES = P1<Math::Vector<Real>>;
using RealGridFunction = GridFunction<RealFES>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

int main(int, char**)
{
  Eigen::initParallel();
  Eigen::setNbThreads(8);
  MMG::Mesh mesh;
  // mesh.load("../resources/mmg/MechanicalPiece_2.mesh", IO::FileFormat::MEDIT);
  //mesh.scale(1.0 / 10);
  mesh.load("out/dOmega.245.mesh", IO::FileFormat::MEDIT);
  mesh.save("miaow.mesh");
  std::exit(1);
  for (auto it = mesh.getFace(); it; ++it)
  {
    if (it->getAttribute() == 29 || it->getAttribute() == 24 || it->getAttribute() == 25)
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, Gamma);
    else
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, Fixed);
  }
  mesh.save("Omega0.mesh", IO::FileFormat::MEDIT);

  Real hmax0 = 0.2;
  Real hmax = hmax0;
  Real hmin = hmax0 / 10.0;
  Real hausd = 0.5 * hmin;
  Real hgrad = 1.2;

  Alert::Info() << "Initializing clamp region..." << Alert::Raise;
  {
    P1 vh(mesh);

    GridFunction dist(vh);
    dist = [&](const Point& p)
    {
      Math::SpatialVector<Real> c(3);
      c << 9, 2, -5;
      return (p - c).norm() - 1;
    };

    mesh = MMG::ImplicitDomainMesher().setAngleDetection(true)
                                      .split(Gamma, { Clamp, Gamma })
                                      .noSplit(Fixed)
                                      .noSplit(Clamp)
                                      .noSplit(Locator)
                                      .noSplit(GammaT)
                                      .setHMax(hmax)
                                      .setHMin(hmin)
                                      .setGradation(hgrad)
                                      .setHausdorff(hausd)
                                      .surface()
                                      .discretize(dist);
  }

  Alert::Info() << "Initializing tool region..." << Alert::Raise;
  {
    P1 vh(mesh);

    GridFunction dist(vh);
    dist = [&](const Point& p)
    {
      Math::SpatialVector<Real> c(3);
      c << -5.5, 0, -5;
      return (p - c).norm() - 0.5;
    };

    mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                      .split(Fixed, { GammaT, Fixed })
                                      .noSplit(Gamma)
                                      .noSplit(Clamp)
                                      .noSplit(Locator)
                                      .noSplit(GammaT)
                                      .setHMax(hmax)
                                      .setHMin(hmin)
                                      .setHausdorff(hausd)
                                      .surface()
                                      .discretize(dist);
  }

  mesh.save("Omega0.mesh", IO::FileFormat::MEDIT);

  std::ofstream fObj("obj.txt");
  size_t i = 0;
  size_t regionCount;
  Real objective = 0;
  while (i < maxIt)
  {
    const Real hmin = hmax / 10.0;
    const Real hausd = 0.5 * hmin;
    const Real hgrad = 1.2;
    const Real k = 0.5 * (hmax + hmin);
    const Real dt = 2 * k;
    const Real radius = k;

    bool topologicalStepLocator = (i == 0) || (i >= 10 && i < 100 && i % 10 == 0);
    bool topologicalStepClamp = (i == 1) || (i >= 11 && i < 101 && (i - 1) % 10 == 0);
    bool geometricStepLocator = !topologicalStepLocator && (i % 2 == 0) && i > 0;
    bool geometricStepClamp = !topologicalStepClamp && (i - 1) % 2 == 0 && i > 1;

    Alert::Info() << "Iteration: " << i                         << Alert::NewLine
                  << "HMax:      " << Alert::Notation(hmax)     << Alert::NewLine
                  << "HMin:      " << Alert::Notation(hmin)     << Alert::NewLine
                  << "Hausdorff: " << Alert::Notation(hausd)    << Alert::NewLine
                  << "HGrad:     " << Alert::Notation(hgrad)    << Alert::NewLine
                  << "dt:        " << Alert::Notation(dt)       << Alert::Raise;

    try
    {
      Alert::Info() << "Optimizing the domain..." << Alert::Raise;
      MMG::Optimizer().setHMax(hmax)
                      .setHMin(hmin)
                      .setHausdorff(hausd)
                      .setGradation(hgrad)
                      .setAngleDetection(false)
                      .optimize(mesh);
    }
    catch (Alert::Exception& e)
    {
      Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
      hmax = 0.5 * hmax < hmax0 ? hmax0 : hmax * 0.9;
      continue;
    }

    Alert::Info() << "Computing required connectivity..." << Alert::Raise;
    mesh.getConnectivity().compute(2, 3); // Computes boundary
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(1, 2);

    Alert::Info() << "RMC..." << Alert::Raise;
    regionCount = rmc(mesh, { Clamp, Locator }, Gamma);

    Alert::Info() << "Found " << Alert::Notation(regionCount) << " regions."
                  << Alert::Raise;

    Alert::Info() << "Skinning mesh..." << Alert::Raise;
    auto dOmega = mesh.skin();
    dOmega.trace({
        {{ Clamp, Gamma }, dClamp },
        {{ Locator, Gamma }, dLocator }});

    dOmega.save("dOmega.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "Building finite element spaces..." << Alert::Raise;
    RealFES sfes(mesh);
    VectorFES vfes(mesh, mesh.getSpaceDimension());
    RealFES dsfes(dOmega);
    VectorFES dvfes(dOmega, dOmega.getSpaceDimension());

    Alert::Info() << "Distancing clamp and locator domains..." << Alert::Raise;
    auto distClamp =
      MMG::Distancer(dsfes).setInteriorDomain(Clamp)
                           .distance(dOmega);
    auto distLocator =
      MMG::Distancer(dsfes).setInteriorDomain(Locator)
                           .distance(dOmega);

    // Parameters
    const Real lambda = 0.5769, mu = 0.3846;

    auto f = -0.1 * BoundaryNormal(mesh); // Locator
    VectorFunction g{1, 0, 0}; // GammaT

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

    RealFunction heClamp =
      [&](const Geometry::Point& p) { return h(distClamp(p) / epsilon) / epsilon; };

    RealFunction heLocator =
      [&](const Geometry::Point& p) { return h(distLocator(p) / epsilon) / epsilon; };

    TrialFunction u(vfes);
    TestFunction  v(vfes);
    Problem state(u, v);
    state = LinearElasticityIntegral(u, v)(lambda, mu)
          - BoundaryIntegral(heClamp * u, v).over(Gamma, Clamp)
          - BoundaryIntegral(f, v).over(Locator)
          - BoundaryIntegral(g, v).over(GammaT)
          ;
    state.assemble();
    Alert::Info() << "Solving state equation..." << Alert::Raise;
    Solver::CG(state).solve();

    u.getSolution().save("u.gf");
    mesh.save("u.mesh");

    Alert::Info() << "Solving adjoint equation..." << Alert::Raise;
    TrialFunction p(vfes);
    TestFunction  q(vfes);
    Problem adjoint(p, q);
    adjoint = LinearElasticityIntegral(p, q)(lambda, mu)
            + Integral(u.getSolution() / mesh.getVolume(), q)
            - BoundaryIntegral(heClamp * p, q).over(Gamma, Clamp);
    Solver::CG(adjoint).solve();

    p.getSolution().save("p.gf");
    mesh.save("p.mesh");

    Alert::Info() << "Computing objective..." << Alert::Raise;
    RealGridFunction j(sfes);
    j = 0.5 * Pow(Frobenius(u.getSolution()), 2)/ mesh.getVolume();
    j.setWeights();

    const Real J = Integral(j).compute();
    const Real pLocator = ellLocator * mesh.getPerimeter(Locator);
    const Real pClamp = ellClamp * mesh.getPerimeter(Clamp);
    objective = J + pClamp + pLocator;

    Alert::Info() << "Objective: " << Alert::Notation(objective) << Alert::NewLine
                  << "J: " << J << Alert::NewLine
                  << "pClamp: " << pClamp << Alert::NewLine
                  << "pLocator: " << pLocator
                  << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

    if (topologicalStepLocator)
    {
      Alert::Info() << "Computing topological sensitivity for the locator..." << Alert::Raise;
      TrialFunction s(dsfes);
      TestFunction  t(dsfes);
      Problem topo(s, t);
      topo = alpha * alpha * Integral(Grad(s), Grad(t))
           + Integral(s, t)
           - Integral(Dot(f, p.getSolution()), t)
           + tgv * Integral(s, t).over(Fixed, GammaT, Clamp, Locator);
      Solver::CG(topo).solve();

      s.getSolution().save("s.gf");
      dsfes.getMesh().save("s.mesh");

      Alert::Info() << "Computing nucleation locations..." << Alert::Raise;
      auto cs = locations(s.getSolution());

      Alert::Info() << "Nucleating " << Alert::Notation(cs.size()) << " holes." << Alert::Raise;
      holes(radius, distLocator, cs);
    }
    else if (topologicalStepClamp)
    {
      Alert::Info() << "Computing topological sensitivity for the clamp..." << Alert::Raise;
      TrialFunction s(dsfes);
      TestFunction  t(dsfes);
      Problem topo(s, t);
      topo = alpha * alpha * Integral(Grad(s), Grad(t))
           + Integral(s, t)
           + Integral(Dot(u.getSolution(), p.getSolution()), t)
           + tgv * Integral(s, t).over(Fixed, GammaT, Clamp, Locator);
      Solver::CG(topo).solve();

      s.getSolution().save("s.gf");
      dsfes.getMesh().save("s.mesh");

      Alert::Info() << "Computing nucleation locations..." << Alert::Raise;
      auto cs = locations(s.getSolution());

      Alert::Info() << "Nucleating " << Alert::Notation(cs.size()) << " holes." << Alert::Raise;
      holes(radius, distClamp, cs);
    }
    else if (geometricStepLocator || geometricStepClamp)
    {
      Alert::Info() << "Computing conormal..." << Alert::Raise;
      GridFunction conormalClamp(dvfes);
      conormalClamp = Grad(distClamp);
      conormalClamp.stableNormalize();

      GridFunction conormalLocator(dvfes);
      conormalLocator = Grad(distLocator);
      conormalLocator.stableNormalize();

      Alert::Info() << "Computing shape gradient..." << Alert::Raise;
      TrialFunction theta(dvfes);
      TestFunction  w(dvfes);
      Problem hilbert(theta, w);
      hilbert = Integral(alpha * alpha * Jacobian(theta), Jacobian(w))
              + Integral(theta, w)
              + 1.0 / epsilon * FaceIntegral(
                  Dot(u.getSolution(), p.getSolution()), Dot(conormalClamp, w)).over(dClamp)
              - FaceIntegral(
                  Dot(f.traceOf(Locator), p.getSolution()), Dot(conormalLocator, w)).over(dLocator)
              + ellLocator * FaceIntegral(conormalLocator, w).over(dLocator)
              + ellClamp * FaceIntegral(conormalClamp, w).over(dClamp);
      Solver::CG(hilbert).solve();

      GridFunction norm(dsfes);
      norm = Frobenius(theta.getSolution());
      theta.getSolution() /= norm.max();

      dOmega.save("Theta.mesh");
      theta.getSolution().save("Theta.gf");

      if (geometricStepLocator)
      {
        Alert::Info() << "Advecting locator..." << Alert::Raise;
        MMG::Advect(distLocator, theta.getSolution()).step(dt);
      }
      else if (geometricStepClamp)
      {
        Alert::Info() << "Advecting clamp..." << Alert::Raise;
        MMG::Advect(distClamp, theta.getSolution()).step(dt);
      }
      else
      {
        Alert::Exception() << "Unhandled case." << Alert::Raise;
      }
    }
    else
    {
      Alert::Exception() << "Unhandled case." << Alert::Raise;
    }

    if (geometricStepLocator || topologicalStepLocator)
    {
      Alert::Info() << "Meshing the anode..." << Alert::Raise;
      try
      {
        GridFunction workaround(sfes);
        workaround.projectOnBoundary(distLocator);
        mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                          .split(Gamma, { Locator, Gamma })
                                          .split(Locator, { Locator, Gamma })
                                          .noSplit(Clamp)
                                          .noSplit(Fixed)
                                          .noSplit(GammaT)
                                          .setHMax(hmax)
                                          .setHMin(hmin)
                                          .setGradation(hgrad)
                                          .setHausdorff(hausd)
                                          .surface()
                                          .discretize(workaround);
        hmax = 1.1 * hmax > hmax0 ? hmax0 : hmax * 1.1;
      }
      catch (Alert::Exception& e)
      {
        Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
        hmax = 0.9 * hmax < 0.1 ? 0.1 : hmax * 0.9;
        continue;
      }
    }
    else if (geometricStepClamp || topologicalStepClamp)
    {
      Alert::Info() << "Meshing the cathode..." << Alert::Raise;
      try
      {
        GridFunction workaround(sfes);
        workaround.projectOnBoundary(distClamp);
        mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                          .split(Gamma, { Clamp, Gamma })
                                          .split(Clamp, { Clamp, Gamma })
                                          .noSplit(Locator)
                                          .noSplit(Fixed)
                                          .noSplit(GammaT)
                                          .setHMax(hmax)
                                          .setHMin(hmin)
                                          .setGradation(hgrad)
                                          .setHausdorff(hausd)
                                          .surface()
                                          .discretize(workaround);
        hmax = 1.1 * hmax > hmax0 ? hmax0 : hmax * 1.1;
      }
      catch (Alert::Exception& e)
      {
        Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
        hmax = 0.9 * hmax < 0.1 ? 0.1 : hmax * 0.9;
        continue;
      }
    }
    else
    {
      Alert::Exception() << "Unhandled case." << Alert::Raise;
    }

    Alert::Info() << "Saving files..." << Alert::Raise;
    mesh.save("Omega.mesh", IO::FileFormat::MEDIT);
    dOmega.save("out/dOmega." + std::to_string(i) +  ".mesh", IO::FileFormat::MEDIT);
    dOmega.save("out/dOmega.mfem." + std::to_string(i) +  ".mesh", IO::FileFormat::MFEM);

    Alert::Success() << "Completed Iteration: " << i << '\n' << Alert::Raise;
    i++;
  }
  return 0;
}
