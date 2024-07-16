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
#include <Rodin/Solver/AppleAccelerate.h>

#include "Tools.h"

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using namespace Rodin::Examples::BoundaryOptimization;

// Parameters
static constexpr Geometry::Attribute Gamma = 1;
static constexpr Geometry::Attribute Unsupported = 2;
static constexpr Geometry::Attribute Support = 11;
static constexpr Geometry::Attribute dSupport = 111;

static constexpr size_t maxIt = 2000;

static constexpr Real epsilon = 1e-6;
static constexpr Real ellP = 0;
static constexpr Real tgv = std::numeric_limits<float>::max();
static constexpr Real alpha = 2;

static Real bA = epsilon;
static Real bTarget = 1.0 / epsilon;
static Real ellA = 0;
static Real targetArea = NAN;

using RealFES = P1<Real>;
using VectorFES = P1<Math::Vector<Real>>;
using RealGridFunction = GridFunction<RealFES>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

int main(int, char**)
{
  Eigen::initParallel();
  Eigen::setNbThreads(8);
  std::cout << Eigen::nbThreads() << std::endl;
  MMG::Mesh mesh;
  // mesh.load("Omega0.mesh", IO::FileFormat::MEDIT);
  // mesh.load("Mechanical.mesh", IO::FileFormat::MEDIT);
  mesh.load("WaterTankVolume.o.mesh", IO::FileFormat::MEDIT);
  // mesh.load("Omega.mesh", IO::FileFormat::MEDIT);

  // mesh.load("MechanicalNoMat.mesh", IO::FileFormat::MEDIT);

  MMG::Mesh D1;
  D1.load("D1.medit.o.mesh", IO::FileFormat::MEDIT);

  P1 h1d1(D1, 3);
  GridFunction phi1(h1d1);
  GridFunction phi2(h1d1);
  GridFunction phi3(h1d1);
  phi1.load("Phi_1.gf");
  phi2.load("Phi_2.gf");
  phi3.load("Phi_3.gf");

  Math::Matrix<Real> aniso(3, 3);
  {
    P1 h1d1s(D1);
    GridFunction integ(h1d1s);
    integ = phi1.x();
    integ.setWeights();
    aniso(0, 0) = Integral(integ).compute();

    integ = phi1.y();
    integ.save("integ.gf");
    integ.setWeights();
    aniso(0, 1) = Integral(integ).compute();


    integ = phi1.z();
    integ.setWeights();
    aniso(0, 2) = Integral(integ).compute();

    integ = phi2.x();
    integ.setWeights();
    aniso(1, 0) = Integral(integ).compute();

    integ = phi2.y();
    integ.setWeights();
    aniso(1, 1) = Integral(integ).compute();

    integ = phi2.z();
    integ.setWeights();
    aniso(1, 2) = Integral(integ).compute();

    integ = phi3.x();
    integ.setWeights();
    aniso(2, 0) = Integral(integ).compute();

    integ = phi3.y();
    integ.setWeights();
    aniso(2, 1) = Integral(integ).compute();

    integ = phi3.z();
    integ.setWeights();
    aniso(2, 2) = Integral(integ).compute();
  }

  std::cout << std::endl << aniso << std::endl;

  D1.save("miaow.mesh", IO::FileFormat::MEDIT);
  phi1.save("miaow.sol", IO::FileFormat::MEDIT);

  // Threads::getGlobalThreadPool().reset(6);

  Real hmax = 0.4;
  Real hmin = hmax / 10.0;
  Real hausd = 0.5 * hmin;
  Real hgrad = 1.2;

  Alert::Info() << "Initializing unsupported region..." << Alert::Raise;
  {
    P1 vh(mesh);

    GridFunction dist(vh);
    dist = [&](const Point& p)
    {
      Math::SpatialVector<Real> c(3);
      c << 0, 0, 5;
      // return (p - c).norm() - (6 + 0.622876);
      return (p - c).norm() - 6.5;
    };

    mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                      .split(Gamma, { Unsupported, Gamma })
                                      .split(Unsupported, { Unsupported, Gamma })
                                      .setHMax(hmax)
                                      .setHMin(hmin)
                                      .setHausdorff(hausd)
                                      .surface()
                                      .discretize(dist);
  }

  Alert::Info() << "Initializing support region..." << Alert::Raise;
  {
    P1 vh(mesh);

    // std::vector<Math::SpatialVector> cs;
    // for (double x = -8; x < 8; x += 2)
    // {
    //   for (double y = -8; y < 8; y += 2)
    //   {
    //     for (double z = -6; z < 4; z += 2)
    //     {
    //       Math::SpatialVector c(3);
    //       c << x, y, z;
    //       cs.push_back(std::move(c));
    //     }
    //   }
    // }

    GridFunction dist(vh);
    dist = [&](const Point& p)
    {
      // double d = (p - cs[0]).norm() - sqrt(2) / 2.0;
      // for (size_t i = 1; i < cs.size(); i++)
      // {
      //   d = std::min(d, (p - cs[i]).norm() - sqrt(2) / 2.0);
      // }
      // return d;
      Math::SpatialVector<Real> c(3);
      c << 0, 0, -5.86;
      return (p - c).norm() - 0.5 * alpha * (hmax + hmin);
    };

    // dist *= -1.0;

    mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                      .split(Gamma, { Support, Gamma })
                                      .split(Support, { Support, Gamma })
                                      .noSplit(Unsupported)
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
  Real augmented = 0, oldAugmented = 1e+5;
  Real objective = 0, oldObjective = 1e+5;
  Real constraint = 0, oldConstraint = 1e+5;
  while (i < maxIt)
  {
    bool topologicalStep = (i < 20) || (i < 200  && i % 10 == 0);
    bool geometricStep = !topologicalStep;

    hmin = hmax / 10.0;
    hausd = hmin;
    hgrad = 1.2;

    const Real k = 0.5 * (hmax + hmin);
    const Real dt = k;
    const Real radius = k;

    Alert::Info() << "Iteration: " << i                         << Alert::NewLine
                  << "HMax:      " << Alert::Notation(hmax)     << Alert::NewLine
                  << "HMin:      " << Alert::Notation(hmin)     << Alert::NewLine
                  << "Hausdorff: " << Alert::Notation(hausd)    << Alert::NewLine
                  << "HGrad:     " << Alert::Notation(hgrad)    << Alert::NewLine
                  << "bA:        " << Alert::Notation(bA)    << Alert::NewLine
                  << "ellA:      " << Alert::Notation(ellA)    << Alert::NewLine
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
      hmax = 0.5 * hmax < 0.1 ? 0.1 : hmax * 0.9;
      continue;
    }

    Alert::Info() << "Computing required connectivity..." << Alert::Raise;
    mesh.getConnectivity().compute(2, 3); // Computes boundary
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(1, 0);

    targetArea = 0.2 * mesh.getPerimeter();

    Alert::Info() << "Target area: " << Alert::Notation(targetArea) << Alert::Raise;

    Alert::Info() << "RMC..." << Alert::Raise;
    regionCount = rmc(mesh, { Support }, Gamma);

    Alert::Info() << "Found " << Alert::Notation(regionCount) << " regions."
      << Alert::Raise;

    Alert::Info() << "Skinning mesh..." << Alert::Raise;
    auto dOmega = mesh.skin();
    dOmega.trace({
        {{ Support, Gamma }, dSupport },
        {{ Unsupported, Support }, dSupport }
        });

    dOmega.save("dOmega.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "Building finite element spaces..." << Alert::Raise;
    RealFES sfes(mesh);
    VectorFES vfes(mesh, mesh.getSpaceDimension());
    RealFES dsfes(dOmega);
    VectorFES dvfes(dOmega, dOmega.getSpaceDimension());


    Alert::Info() << "Distancing support domain..." << Alert::Raise;
    auto dist = MMG::Distancer(dsfes).setInteriorDomain(Support)
                                     .distance(dOmega);

    // Parameters
    const Real lambda = 0.5769, mu = 0.3846;

    VectorFunction g{0, 0, -0.01};

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

    RealFunction he =
      [&](const Geometry::Point& p) { return h(dist(p) / epsilon) / epsilon; };

    Alert::Info() << "Solving state equation..." << Alert::Raise;
    TrialFunction u(vfes);
    TestFunction  v(vfes);
    Problem state(u, v);
    state = LinearElasticityIntegral(u, v)(lambda, mu)
          - Integral(g, v)
          + FaceIntegral(he * u, v).over(Support);
          ;
    Solver::CG(state).solve();

    u.getSolution().save("u.gf");
    mesh.save("u.mesh");

    Alert::Info() << "Solving adjoint equation..." << Alert::Raise;
    TrialFunction p(vfes);
    TestFunction  q(vfes);
    Problem adjoint(p, q);
    adjoint = LinearElasticityIntegral(p, q)(lambda, mu)
            + Integral(u.getSolution() / mesh.getVolume(), q)
            + FaceIntegral(he * p, q).over(Support);
    Solver::CG(adjoint).solve();

    p.getSolution().save("p.gf");
    mesh.save("p.mesh");

    Alert::Info() << "Computing objective..." << Alert::Raise;
    RealGridFunction j(sfes);
    j = Frobenius(u.getSolution()) / mesh.getVolume();
    j.setWeights();
    if (i > 0)
    {
      oldAugmented = augmented;
      oldObjective = objective;
      oldConstraint = constraint;
    }

    const Real J = Integral(j).compute();
    const Real area = mesh.getPerimeter(Support);
    const Real perimeter = dOmega.getMeasure(1, dSupport);
    objective = J;
    constraint = (area / targetArea - 1);
    augmented =
      objective + ellA * constraint + 0.5 * bA * constraint * constraint;

    Alert::Info() << "Objective: " << Alert::Notation(objective) << Alert::NewLine
                  << "Augmented: " << Alert::Notation(augmented) << Alert::NewLine
                  << "Support Area: " << Alert::Notation(area) << Alert::NewLine
                  << "Constraint: " << Alert::Notation(constraint) << Alert::NewLine
                  << Alert::Raise;
    fObj << augmented << "\n";
    fObj.flush();

    if (topologicalStep)
    {
      Alert::Info() << "Topological optimization..." << Alert::Raise;
      Alert::Info() << "Inserting supporting region..." << Alert::Raise;
      TrialFunction s(dsfes);
      TestFunction  t(dsfes);
      Problem topo(s, t);
      topo = alpha * alpha * dt * dt * Integral(Grad(s), Grad(t))
           + Integral(s, t)
           + Integral((u.getSolution().T() * aniso * p.getSolution()).coeff(0, 0), t)
           + tgv * Integral(s, t).over(Support, Unsupported);
      Solver::CG(topo).solve();

      s.getSolution().save("Topo.gf");
      dsfes.getMesh().save("Topo.mesh");

      GridFunction raw(dsfes);
      raw = -(u.getSolution().T() * aniso * p.getSolution()).coeff(0, 0);
      raw.save("Raw.gf");
      dsfes.getMesh().save("Raw.mesh");

      Alert::Info() << "Computing nucleation locations..." << Alert::Raise;
      const Real tc = s.getSolution().max();
      std::vector<Point> cs;
      for (auto it = dOmega.getVertex(); !it.end(); ++it)
      {
        const Point p(*it, it->getTransformation(),
            Polytope::getVertex(0, Polytope::Type::Point), it->getCoordinates());
        const Real tp = s.getSolution()(p);
        if (tp > 0 && (tp / tc) > (1 - 1e-12))
        {
          cs.emplace_back(std::move(p));
          break;
        }
      }

      Alert::Info() << "Nucleating " << Alert::Notation(cs.size()) << " holes..."
                    << Alert::Raise;
      holes(radius, dist, cs);
    }
    else if (geometricStep)
    {
      Alert::Info() << "Geometrical optimization..." << Alert::Raise;

      Alert::Info() << "Computing conormal..." << Alert::Raise;
      GridFunction conormal(dvfes);
      conormal.projectOnFaces(Average(Grad(dist)));
      conormal.stableNormalize();

      conormal.getFiniteElementSpace().getMesh().save("Conormal.mesh");
      conormal.save("Conormal.gf");

      Alert::Info() << "Computing shape gradient..." << Alert::Raise;

      TrialFunction theta(dvfes);
      TestFunction  w(dvfes);
      Problem hilbert(theta, w);
      hilbert = alpha * alpha * dt * dt * Integral(Jacobian(theta), Jacobian(w))
              + Integral(theta, w)
              + 1.0 / epsilon * FaceIntegral(
                  Dot(u.getSolution(), p.getSolution()), Dot(conormal, w)).over(dSupport)
              + constraint * ellA * FaceIntegral(conormal, w).over(dSupport)
              + 0.5 * 2 * bA * constraint * FaceIntegral(conormal, w).over(dSupport)
              // + ellP * FaceIntegral(Div(conormal).traceOf(Support) * conormal, w).over(dSupport)
              // + tgv * Integral(theta, w).over(Unsupported)
              ;
      Solver::CG(hilbert).solve();

      GridFunction norm(dsfes);
      norm = Frobenius(theta.getSolution());
      theta.getSolution() /= norm.max();

      dOmega.save("Theta.mesh");
      theta.getSolution().save("Theta.gf");

      Alert::Info() << "Advecting support..." << Alert::Raise;
      MMG::Advect(dist, theta.getSolution()).step(dt);
    }

    dist.save("dist.gf");
    dist.getFiniteElementSpace().getMesh().save("dist.mesh");

    Alert::Info() << "Meshing the support..." << Alert::Raise;
    try
    {
      GridFunction workaround(sfes);
      workaround.projectOnBoundary(dist);
      mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                        .split(Support, { Support, Gamma })
                                        .split(Gamma, { Support, Gamma })
                                        .noSplit(Unsupported)
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
      hmax = 0.9 * hmax < 0.2 ? 0.5 : hmax * 0.9;
      continue;
    }

    dist.save("dist.gf");
    dist.getFiniteElementSpace().getMesh().save("dist.mesh");

    Alert::Info() << "Saving files..." << Alert::Raise;
    mesh.save("Omega.mesh", IO::FileFormat::MEDIT);
    dOmega.save("out/dOmega." + std::to_string(i) +  ".mesh", IO::FileFormat::MEDIT);
    dOmega.save("out/dOmega.mfem." + std::to_string(i) +  ".mesh", IO::FileFormat::MFEM);

    Alert::Success() << "Completed Iteration: " << i << '\n' << Alert::Raise;

    ellA += bA * constraint;
    if (i > 0 && !topologicalStep)
      bA *= 1 + alpha * alpha * abs(oldAugmented - augmented);
    bA = fmin(bTarget, bA);

    i++;
  }
  return 0;
}

