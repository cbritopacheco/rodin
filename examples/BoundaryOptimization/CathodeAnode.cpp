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

static constexpr Scalar epsilon = 0.001;
static constexpr Scalar ellAnode = 1e-4;
static constexpr Scalar ellCathode = 1e-4;
static constexpr Scalar lambdaAnode = 1e-5;
static constexpr Scalar lambdaCathode = 1e-5;
static constexpr Scalar tgv = std::numeric_limits<float>::max();
static const Scalar alpha = 4;

using ScalarFES = P1<Scalar>;
using VectorFES = P1<Math::Vector<Scalar>>;
using ScalarGridFunction = GridFunction<ScalarFES>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

int main(int, char**)
{
  MMG::Mesh mesh;
  // mesh.load("Unnamed1-Fusion001.o.mesh", IO::FileFormat::MEDIT);
  // mesh.save("Omega0.mesh");
  mesh.load("Omega.mesh", IO::FileFormat::MEDIT);

  Scalar hmax = 0.5;

  // Alert::Info() << "Initializing anode region..." << Alert::Raise;
  // {
  //   P1 vh(mesh);

  //   GridFunction distAnode(vh);
  //   distAnode = [&](const Point& p)
  //   {
  //     Math::SpatialVector c(3);
  //     c << 2.5, 35, 7;
  //     return (p - c).norm() - 1;
  //   };

  //   mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
  //                                     .split(Anode, { Anode, Gamma })
  //                                     .split(Gamma, { Anode, Gamma })
  //                                     .noSplit(Cathode)
  //                                     .noSplit(Inlet)
  //                                     .setHMax(hmax)
  //                                     .surface()
  //                                     .discretize(distAnode);
  // }

  // Alert::Info() << "Initializing cathode region..." << Alert::Raise;
  // {
  //   P1 vh(mesh);

  //   GridFunction distCathode(vh);
  //   distCathode = [&](const Point& p)
  //   {
  //     Math::SpatialVector c(3);
  //     c << 2.5, 15, -2.5;
  //     return (p - c).norm() - 1;
  //   };

  //   mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
  //                                     .split(Cathode, { Cathode, Gamma })
  //                                     .split(Gamma, { Cathode, Gamma })
  //                                     .noSplit(Anode)
  //                                     .noSplit(Inlet)
  //                                     .setHMax(hmax)
  //                                     .surface()
  //                                     .discretize(distCathode);
  // }

  std::ofstream fObj("obj.txt");
  size_t i = 208;
  size_t regionCount;
  while (i < maxIt)
  {
    const Scalar hmin = hmax / 10.0;
    const Scalar hausd = 0.5 * hmin;
    const Scalar hgrad = 1.2;
    const Scalar k = 0.5 * (hmax + hmin);
    const Scalar dt = 2 * k;
    const Scalar radius = k;

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

    Alert::Info() << "RMC..." << Alert::Raise;
    regionCount = rmc(mesh, { Cathode, Anode }, Gamma);

    Alert::Info() << "Found " << Alert::Notation(regionCount) << " regions."
      << Alert::Raise;

    Alert::Info() << "Skinning mesh..." << Alert::Raise;
    auto dOmega = mesh.skin();
    dOmega.trace({{{Anode, Gamma}, dAnode}, {{Cathode, Gamma}, dCathode}});
    dOmega.save("dOmega.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "Building finite element spaces..." << Alert::Raise;
    ScalarFES sfes(mesh);
    VectorFES vfes(mesh, mesh.getSpaceDimension());
    ScalarFES dsfes(dOmega);
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

    Solver::CG cg;

    // Parameters
    ScalarFunction uIn = 1.0;
    ScalarFunction gamma = 1.0;

    auto h = [](Scalar r)
    {
      if (r < -1.0)
        return 1.0;
      else if (r > 1.0)
        return 0.0;
      else
        return 1.0 - 1.0 / (1.0 + std::exp(4 * r / (r * r - 1.0)));
    };

    ScalarFunction heAnode =
      [&](const Geometry::Point& p) { return h(distAnode(p) / epsilon) / epsilon; };
    ScalarFunction heCathode =
      [&](const Geometry::Point& p) { return h(distCathode(p) / epsilon) / epsilon; };

    TrialFunction u(sfes);
    TestFunction  v(sfes);
    Problem state(u, v);
    state = Integral(gamma * Grad(u), Grad(v))
          + BoundaryIntegral((heAnode + heCathode) * u, v)
          - BoundaryIntegral(heAnode * uIn, v)
          ;
    Alert::Info() << "Solving state equation..." << Alert::Raise;
    state.solve(cg);
    u.getSolution().save("U.gf");
    mesh.save("Omega.mfem.mesh");

    TrialFunction gtr(vfes);
    TestFunction gte(vfes);
    Problem potential(gtr, gte);
    potential = Integral(alpha * alpha * Jacobian(gtr), Jacobian(gte))
              + Integral(gtr, gte)
              - Integral(0.5 * Grad(u.getSolution()) / mesh.getVolume(), gte)
              ;
    potential.solve(cg);

    Alert::Info() << "Solving adjoint equation..." << Alert::Raise;
    TrialFunction p(sfes);
    TestFunction  q(sfes);
    Problem adjoint(p, q);
    adjoint = Integral(gamma * Grad(p), Grad(q))
            + BoundaryIntegral((heAnode + heCathode) * p, q)
            - Integral(gtr.getSolution(), Grad(q));
    adjoint.solve(cg);
    p.getSolution().save("P.gf");

    GridFunction sC(sfes);
    sC = u.getSolution() * p.getSolution();
    sC.save("sC.gf");

    GridFunction sA(sfes);
    sA = (uIn - u.getSolution()) * p.getSolution();
    sA.save("sA.gf");

    Alert::Info() << "Computing objective..." << Alert::Raise;
    GridFunction j(sfes);
    j = 0.5 * Frobenius(Grad(u.getSolution()));
    j.setWeights();
    mesh.save("j.mesh");
    j.save("j.gf");
    const Scalar Ij = Integral(j).compute();
    const Scalar objective =
      Ij - ellAnode * mesh.getPerimeter(Anode) - ellCathode * mesh.getPerimeter(Cathode);
    Alert::Info() << "Objective: " << Alert::Notation(objective) << Alert::NewLine
                  << "Norm: " << Alert::Notation(Ij)
                  << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

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
            - 1.0 / epsilon * FaceIntegral(
                u.getSolution() * p.getSolution(), Dot(conormalAnode, w)).over(dAnode)
            + 1.0 / epsilon * FaceIntegral(
                (uIn - u.getSolution()) * p.getSolution(), Dot(conormalCathode, w)).over(dCathode)
            + ellAnode * FaceIntegral(conormalAnode, w).over(dAnode)
            + ellCathode * FaceIntegral(conormalCathode, w).over(dCathode)
            + lambdaAnode * FaceIntegral(
                Div(conormalAnode).traceOf(Anode) * conormalAnode, w).over(dAnode)
            + lambdaCathode * FaceIntegral(
                Div(conormalCathode).traceOf(Cathode) * conormalCathode, w).over(dCathode)
            + tgv * Integral(theta, w).over(Inlet, FixedAnode, FixedCathode);
    hilbert.solve(cg);

    GridFunction norm(dsfes);
    norm = Frobenius(theta.getSolution());
    theta.getSolution() /= norm.max();

    dOmega.save("Theta.mesh");
    theta.getSolution().save("Theta.gf");

    // dOmega.save("out/Theta." + std::to_string(i) + ".mesh");
    // theta.getSolution().save("out/Theta." + std::to_string(i) + ".gf");

    // if (true)
    if (((i % 10 == 0) || ((i - 1) % 10 == 0)) && i < 100)
    {
      Alert::Info() << "Computing topological sensitivity..." << Alert::Raise;
      TrialFunction s(dsfes);
      TestFunction  t(dsfes);
      Problem topo(s, t);
      if (i % 10 == 0)
      {
        topo = Integral(alpha * alpha * Grad(s), Grad(t))
             + Integral(s, t)
             - Integral((uIn - u.getSolution()) * p.getSolution(), t)
             + tgv * Integral(s, t).over(Inlet, Cathode, Anode, FixedAnode, FixedCathode);
      }
      else
      {
        topo = Integral(alpha * alpha * Grad(s), Grad(t))
             + Integral(s, t)
             + Integral((u.getSolution()) * p.getSolution(), t)
             + tgv * Integral(s, t).over(Inlet, Cathode, Anode, FixedAnode, FixedCathode);
      }
      topo.solve(cg);

      s.getSolution().save("Topo.sol", IO::FileFormat::MEDIT);
      dsfes.getMesh().save("Topo.mesh", IO::FileFormat::MEDIT);

      // s.getSolution().save("out/Topo." + std::to_string(i) + ".gf");
      // dsfes.getMesh().save("out/Topo." + std::to_string(i) + ".mesh");

      Alert::Info() << "Computing nucleation locations..." << Alert::Raise;
      const Scalar tc = s.getSolution().max();
      std::vector<Point> cs;
      for (auto it = dOmega.getVertex(); !it.end(); ++it)
      {
        const Point p(*it, it->getTransformation(),
            Polytope::getVertices(Polytope::Type::Point).col(0), it->getCoordinates());
        const Scalar tp = s.getSolution()(p);
        if (tp > 0 && (tp / tc) > (1 - 1e-12))
        {
          cs.emplace_back(std::move(p));
          break;
        }
      }

      Alert::Info() << "Nucleating " << Alert::Notation(cs.size()) << " holes..."
                    << Alert::Raise;
      if (cs.size())
      {
        if (i % 10 == 0)
        {
          auto holes =
            [&](const Point& v)
            {
              Scalar d = distAnode(v);
              for (const auto& c : cs)
              {
                const Scalar dd = (v - c).norm() - radius;
                d = std::min(d, dd);
              }
              return d;
            };
          distAnode = holes;
        }
        else
        {
          auto holes =
            [&](const Point& v)
            {
              Scalar d = distCathode(v);
              for (const auto& c : cs)
              {
                const Scalar dd = (v - c).norm() - radius;
                d = std::min(d, dd);
              }
              return d;
            };
          distCathode = holes;
        }
      }
    }

    if (i % 2 == 0)
    {
      Alert::Info() << "Advecting anode..." << Alert::Raise;

      if (i > 1)
        MMG::Advect(distAnode, theta.getSolution()).step(dt);

      // dOmega.save("out/AdvectedAnode." + std::to_string(i) + ".mesh");
      // distAnode.save("out/AdvectedAnode." + std::to_string(i) + ".gf");

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
        hmax = 1.1 * hmax > 0.5 ? 0.5 : hmax * 1.1;
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
      Alert::Info() << "Advecting cathode..." << Alert::Raise;
      if (i > 1)
        MMG::Advect(distCathode, theta.getSolution()).step(dt);

      // dOmega.save("out/AdvectedCathode." + std::to_string(i) + ".mesh");
      // distCathode.save("out/AdvectedCathode." + std::to_string(i) + ".gf");

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
        hmax = 1.1 * hmax > 0.5 ? 0.5 : hmax * 1.1;
      }
      catch (Alert::Exception& e)
      {
        Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
        hmax = 0.9 * hmax < 0.1 ? 0.1 : hmax * 0.9;
        continue;
      }
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

