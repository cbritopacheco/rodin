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

#include <chrono>

#include <Rodin/Models/Hilbert/H1a.h>

#include "Tools.h"

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using namespace Rodin::Examples::BoundaryOptimization;

// Parameters
static constexpr Geometry::Attribute Anode = 2; // u = 1
static constexpr Geometry::Attribute Cathode = 3; // u = 0
static constexpr Geometry::Attribute Gamma = 4; // dn u = 0
static constexpr Geometry::Attribute Inlet = 6; // dn u = 0
static constexpr Geometry::Attribute FixedAnode = 8;
static constexpr Geometry::Attribute FixedCathode = 9;

static constexpr Geometry::Attribute dAnode = 3;
static constexpr Geometry::Attribute dCathode = 2;

static constexpr size_t maxIt = 10000;

static constexpr Real epsilon = 1e-3;
static constexpr Real ellAnode = 1e-3;
static constexpr Real ellCathode = 1e-3;
static constexpr Real lambdaAnode = 1e-3;
static constexpr Real lambdaCathode = 1e-3;
static constexpr Real tgv = 1e+12;
static const Real alpha = 4;

using RealFES = P1<Real>;
using VectorFES = P1<Math::Vector<Real>>;
using RealGridFunction = GridFunction<RealFES>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

int main(int, char**)
{
  MMG::Mesh mesh;
  mesh.load("../resources/examples/BoundaryOptimization/EMM.medit.mesh", IO::FileFormat::MEDIT);
  // mesh.load("Omega.mesh", IO::FileFormat::MEDIT);
  // Threads::getGlobalThreadPool().reset(8);

  Real hmax0 = 0.5;
  Real hmax = hmax0;

  std::ofstream fObj("obj.txt");
  size_t i = 0;
  size_t regionCount;
  while (i < maxIt)
  {
    auto t0 = std::chrono::system_clock::now();

    // Convert the current time to time since epoch
    const Real hmin = hmax / 10.0;
    const Real hausd = 0.5 * hmin;
    const Real hgrad = 1.2;
    const Real k = 0.5 * (hmax + hmin);
    const Real dt = 2 * k;
    const Real radius = k;

    bool topologicalStepAnode = (i == 0) || (i >= 10 && i < 100 && i % 10 == 0);
    bool topologicalStepCathode = (i == 1) || (i >= 11 && i < 101 && (i - 1) % 10 == 0);
    bool geometricStepAnode = !topologicalStepAnode && (i % 2 == 0) && i > 0;
    bool geometricStepCathode = !topologicalStepCathode && (i - 1) % 2 == 0 && i > 1;

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
      hmax = 0.9 * hmax < 0.1 ? 0.1 : hmax * 0.9;
      continue;
    }

    Alert::Info() << "Computing required connectivity..." << Alert::Raise;
    mesh.getConnectivity().compute(2, 3); // Computes boundary
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(1, 2);

    mesh.save("Omega.mfem.mesh");

    Alert::Info() << "RMC..." << Alert::Raise;
    regionCount = rmc(mesh, { Cathode, Anode }, Gamma);

    Alert::Info() << "Found " << Alert::Notation(regionCount) << " regions."
      << Alert::Raise;

    Alert::Info() << "Skinning mesh..." << Alert::Raise;
    auto dOmega = mesh.skin();
    dOmega.trace({
        {{Anode, Gamma}, dAnode}, {{Cathode, Gamma}, dCathode}});
    dOmega.save("dOmega.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "Building finite element spaces..." << Alert::Raise;
    RealFES sfes(mesh);
    VectorFES vfes(mesh, mesh.getSpaceDimension());
    RealFES dsfes(dOmega);
    VectorFES dvfes(dOmega, dOmega.getSpaceDimension());

    Alert::Info() << "Distancing cathode and anode domains..." << Alert::Raise;
    auto distAnode =
      MMG::Distancer(dsfes).setInteriorDomain({ Anode, FixedAnode })
                           .distance(dOmega);
    auto distCathode =
      MMG::Distancer(dsfes).setInteriorDomain({ Cathode, FixedCathode })
                           .distance(dOmega);

    // dOmega.save("out/Cathode." + std::to_string(i) + ".mesh");
    // distCathode.save("out/Cathode." + std::to_string(i) + ".gf");

    // dOmega.save("out/Anode." + std::to_string(i) + ".mesh");
    // distAnode.save("out/Anode." + std::to_string(i) + ".gf");

    // Parameters
    RealFunction uIn = 1.0;
    RealFunction gamma = 1.0;

    auto h = [](Real r)
    {
      if (r < -1.0)
        return 1.0;
      else if (r > 1.0)
        return 0.0;
      else
        return 1.0 - 1.0 / (1.0 + std::exp(4 * r / (r * r - 1.0)));
    };

    RealFunction heAnode =
      [&](const Geometry::Point& p) { return h(distAnode(p) / epsilon) / epsilon; };
    RealFunction heCathode =
      [&](const Geometry::Point& p) { return h(distCathode(p) / epsilon) / epsilon; };

    TrialFunction u(sfes);
    TestFunction  v(sfes);
    Problem state(u, v);
    state = Integral(gamma * Grad(u), Grad(v))
          + BoundaryIntegral((heAnode + heCathode) * u, v)
          - BoundaryIntegral(heAnode * uIn, v)
          ;

    Alert::Info() << "Solving state equation..." << Alert::Raise;
    Solver::CG(state).solve();
    u.getSolution().save("U.gf");

    TrialFunction gtr(vfes);
    TestFunction gte(vfes);
    Problem potential(gtr, gte);
    potential = Integral(alpha * alpha * Jacobian(gtr), Jacobian(gte))
              + Integral(gtr, gte)
              - Integral(Grad(u.getSolution()), gte)
              ;
    Solver::CG(potential).solve();

    gtr.getSolution().save("miaow.gf");

    GridFunction miaow(vfes);
    miaow = VectorFunction{0, 0, 0};

    TrialFunction p(sfes);
    TestFunction  q(sfes);
    Problem adjoint(p, q);
    adjoint = Integral(gamma * Grad(p), Grad(q))
            + BoundaryIntegral((heAnode + heCathode) * p, q)
            + Integral(0.5 * miaow, Grad(q));

    Alert::Info() << "Solving adjoint equation..." << Alert::Raise;
    Solver::CG(adjoint).solve();
    p.getSolution().save("P.gf");

    Alert::Info() << "Computing objective..." << Alert::Raise;
    GridFunction j(sfes);
    j = 0.5 * Pow(Frobenius(Grad(u.getSolution())), 2);
    j.setWeights();

    const Real J = Integral(j).compute();
    const Real objective =
      J - ellAnode * mesh.getPerimeter(Anode) - ellCathode * mesh.getPerimeter(Cathode)
      - lambdaAnode * mesh.getMeasure(1, Anode) - lambdaCathode * mesh.getMeasure(1, Cathode)
      ;
    Alert::Info() << "Objective: " << Alert::Notation(objective) << Alert::NewLine
                  << "J: " << Alert::Notation(J)
                  << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

    if (topologicalStepAnode)
    {
      Alert::Info() << "Computing topological sensitivity for the anode..." << Alert::Raise;
      TrialFunction s(dsfes);
      TestFunction  t(dsfes);
      Problem topo(s, t);
      topo = Integral(alpha * alpha * Grad(s), Grad(t))
           + Integral(s, t)
           + Integral((uIn - u.getSolution()) * p.getSolution(), t)
           + tgv * Integral(s, t).over(Inlet, Cathode, Anode, FixedAnode, FixedCathode);
      Solver::CG(topo).solve();

      dOmega.save("sA.mesh");
      s.getSolution().save("sA.gf");

      Alert::Info() << "Computing nucleation locations..." << Alert::Raise;
      auto cs = locations(s.getSolution());

      Alert::Info() << "Nucleating " << Alert::Notation(cs.size()) << " holes." << Alert::Raise;
      holes(radius, distAnode, cs);
    }
    else if (topologicalStepCathode)
    {
      Alert::Info() << "Computing topological sensitivity for the cathode..." << Alert::Raise;
      TrialFunction s(dsfes);
      TestFunction  t(dsfes);
      Problem topo(s, t);
      topo = Integral(alpha * alpha * Grad(s), Grad(t))
           + Integral(s, t)
           - Integral((u.getSolution()) * p.getSolution(), t)
           + tgv * Integral(s, t).over(Inlet, Cathode, Anode, FixedAnode, FixedCathode);
      Solver::CG(topo).solve();

      dOmega.save("sC.mesh");
      s.getSolution().save("sC.gf");

      Alert::Info() << "Computing nucleation locations..." << Alert::Raise;
      auto cs = locations(s.getSolution());

      Alert::Info() << "Nucleating " << Alert::Notation(cs.size()) << " holes." << Alert::Raise;
      holes(radius, distCathode, cs);
    }
    else if (geometricStepAnode || geometricStepCathode)
    {
      Alert::Info() << "Computing conormal..." << Alert::Raise;
      GridFunction conormalCathode(dvfes);
      conormalCathode = Grad(distCathode);
      conormalCathode.stableNormalize();

      GridFunction conormalAnode(dvfes);
      conormalAnode = Grad(distAnode);
      conormalAnode.stableNormalize();

      Alert::Info() << "Computing shape gradient..." << Alert::Raise;
      TrialFunction theta(dvfes);
      TestFunction  w(dvfes);
      Problem hilbert(theta, w);
      hilbert = Integral(alpha * alpha * Jacobian(theta), Jacobian(w))
              + Integral(theta, w)
              + 1.0 / epsilon * FaceIntegral(
                  (uIn - u.getSolution()) * p.getSolution(), Dot(conormalAnode, w)).over(dAnode)
              - 1.0 / epsilon * FaceIntegral(
                  u.getSolution() * p.getSolution(), Dot(conormalCathode, w)).over(dCathode)
              + ellAnode * FaceIntegral(conormalAnode, w).over(dAnode)
              + ellCathode * FaceIntegral(conormalCathode, w).over(dCathode)
              + lambdaAnode * FaceIntegral(
                  Div(conormalAnode).traceOf(Anode) * conormalAnode, w).over(dAnode)
              + lambdaCathode * FaceIntegral(
                  Div(conormalCathode).traceOf(Cathode) * conormalCathode, w).over(dCathode)
              + tgv * Integral(theta, w).over(Inlet, FixedAnode, FixedCathode);
      Solver::CG(hilbert).solve();

      GridFunction norm(dsfes);
      norm = Frobenius(theta.getSolution());
      theta.getSolution() /= norm.max();

      dOmega.save("Theta.mesh");
      theta.getSolution().save("Theta.gf");

      if (geometricStepAnode)
      {
        Alert::Info() << "Advecting anode..." << Alert::Raise;
        MMG::Advect(distAnode, theta.getSolution()).step(dt);
      }
      else if (geometricStepCathode)
      {
        Alert::Info() << "Advecting cathode..." << Alert::Raise;
        MMG::Advect(distCathode, theta.getSolution()).step(dt);
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

    if (geometricStepAnode || topologicalStepAnode)
    {
      Alert::Info() << "Meshing the anode..." << Alert::Raise;
      try
      {
        GridFunction workaround(sfes);
        workaround.projectOnBoundary(distAnode);
        mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                          .split(Anode, { Anode, Gamma })
                                          .split(Gamma, { Anode, Gamma })
                                          .noSplit(Cathode)
                                          .noSplit(FixedCathode)
                                          .noSplit(FixedAnode)
                                          .noSplit(Inlet)
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
    else if (geometricStepCathode || topologicalStepCathode)
    {
      Alert::Info() << "Meshing the cathode..." << Alert::Raise;
      try
      {
        GridFunction workaround(sfes);
        workaround.projectOnBoundary(distCathode);
        mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                          .split(Cathode, { Cathode, Gamma })
                                          .split(Gamma, { Cathode, Gamma })
                                          .noSplit(Anode)
                                          .noSplit(FixedCathode)
                                          .noSplit(FixedAnode)
                                          .noSplit(Inlet)
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

    auto t1 = std::chrono::system_clock::now();

    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t1 -
        t0).count() << std::endl;
  }

  return 0;
}

