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

#include <Rodin/Models/Hilbert/H1a.h>

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// Parameters
static constexpr Geometry::Attribute Gamma = 1;
static constexpr Geometry::Attribute GammaT = 6;

static constexpr Geometry::Attribute Clamp = 11; // Clamp
static constexpr Geometry::Attribute Locator = 12; // Locator

static constexpr Geometry::Attribute dClamp = 111;
static constexpr Geometry::Attribute dLocator = 6;

static constexpr size_t maxIt = 1000;

static constexpr Scalar epsilon = 1e-6;
static constexpr Scalar ellClamp = 0.001;
static constexpr Scalar ellLocator = 0.001;
static constexpr Scalar radius = 0.2;
static constexpr Scalar tgv = std::numeric_limits<float>::max();
static const Scalar alpha = 4;

using ScalarFES = P1<Scalar, Context::Sequential>;
using VectorFES = P1<Math::Vector, Context::Sequential>;
using ScalarGridFunction = GridFunction<ScalarFES>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

size_t rmc(MeshBase& mesh);
void holes(ScalarGridFunction& dist, const std::vector<Point>& cs);

int main(int, char**)
{
  Eigen::initParallel();
  Eigen::setNbThreads(8);
  std::cout << Eigen::nbThreads() << std::endl;
  MMG::Mesh mesh;
  //mesh.load("Omega0.o.mesh", IO::FileFormat::MEDIT);
  mesh.load("Mechanical.mesh", IO::FileFormat::MEDIT);

  MMG::Mesh D1;
  D1.load("D1.medit.o.mesh", IO::FileFormat::MEDIT);

  P1 h1d1(D1, 3);
  GridFunction phi1(h1d1);
  GridFunction phi2(h1d1);
  GridFunction phi3(h1d1);
  phi1.load("Phi_1.gf");
  phi2.load("Phi_2.gf");
  phi3.load("Phi_3.gf");

  Math::Matrix aniso(3, 3);
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

  // auto philf = [](const Math::Vector& v) {  }

  D1.save("miaow.mesh", IO::FileFormat::MEDIT);
  phi1.save("miaow.sol", IO::FileFormat::MEDIT);

  Threads::getGlobalThreadPool().reset(6);

  Scalar hmax = 0.3;
  Scalar hmin = hmax / 10.0;
  Scalar hausd = 0.5 * hmin;
  Scalar hgrad = 1.2;

  // Alert::Info() << "Initializing clamp region..." << Alert::Raise;
  // {
  //   P1 vh(mesh);

  //   GridFunction dist(vh);
  //   dist = [&](const Point& p)
  //   {
  //     Math::SpatialVector c(3);
  //     c << -1.5, -2, -5;
  //     return (p - c).norm() - 1;
  //   };

  //   mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
  //                                     .split(Gamma, { Clamp, Gamma })
  //                                     .noSplit(Locator)
  //                                     .noSplit(GammaT)
  //                                     .setHMax(hmax)
  //                                     .setHMin(hmin)
  //                                     .setHausdorff(hausd)
  //                                     .surface()
  //                                     .discretize(dist);
  // }

  // Alert::Info() << "Initializing locator region..." << Alert::Raise;
  // {
  //   P1 vh(mesh);

  //   GridFunction dist(vh);
  //   dist = [&](const Point& p)
  //   {
  //     Math::SpatialVector c(3);
  //     c << -1.5, 2, -5;
  //     return (p - c).norm() - 1;
  //   };

  //   mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
  //                                     .split(Gamma, { Locator, Gamma })
  //                                     .noSplit(Clamp)
  //                                     .noSplit(GammaT)
  //                                     .setHMax(hmax)
  //                                     .setHMin(hmin)
  //                                     .setHausdorff(hausd)
  //                                     .surface()
  //                                     .discretize(dist);
  // }

  // Alert::Info() << "Initializing tool region..." << Alert::Raise;
  // {
  //   P1 vh(mesh);

  //   GridFunction dist(vh);
  //   dist = [&](const Point& p)
  //   {
  //     Math::SpatialVector c(3);
  //     c << -5.5, 0, -5;
  //     return (p - c).norm() - 0.5;
  //   };

  //   mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
  //                                     .split(Gamma, { GammaT, Gamma })
  //                                     .noSplit(Clamp)
  //                                     .noSplit(Locator)
  //                                     .setHMax(hmax)
  //                                     .setHMin(hmin)
  //                                     .setHausdorff(hausd)
  //                                     .surface()
  //                                     .discretize(dist);
  // }

  // mesh.save("Omega0.mesh", IO::FileFormat::MEDIT);

  std::ofstream fObj("obj.txt");
  size_t i = 0;
  size_t regionCount;
  Scalar objective = 0, oldObjective = 9999;
  while (i < maxIt)
  {
    bool topologicalStep = i > 2 && i < 100 && ((i % 10) == 0 || ((i - 1) % 10 == 0));
    bool geometricStep = !topologicalStep;
    bool clampStep = i % 2 && false;
    bool locatorStep = !clampStep;

    hmin = hmax / 10.0;
    hausd = 0.5 * hmin;
    hgrad = 1.2;

    const Scalar k = 0.5 * (hmax + hmin);
    const Scalar dt = 0.05 * alpha * k;
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
      hmax = 0.5 * hmax < 0.1 ? 0.1 : hmax * 0.9;
      continue;
    }

    Alert::Info() << "Computing required connectivity..." << Alert::Raise;
    mesh.getConnectivity().compute(2, 3); // Computes boundary
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(1, 2);

    // Alert::Info() << "RMC..." << Alert::Raise;
    // regionCount = rmc(mesh);

    // Alert::Info() << "Found " << Alert::Notation(regionCount) << " regions."
    //   << Alert::Raise;

    Alert::Info() << "Skinning mesh..." << Alert::Raise;
    auto dOmega = mesh.skin();
    dOmega.trace({
        {{ Clamp, Gamma }, dClamp },
        {{ Locator, Gamma }, dLocator }});

    dOmega.save("dOmega.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "Building finite element spaces..." << Alert::Raise;
    ScalarFES sfes(mesh);
    VectorFES vfes(mesh, mesh.getSpaceDimension());
    ScalarFES dsfes(dOmega);
    VectorFES dvfes(dOmega, dOmega.getSpaceDimension());

    Alert::Info() << "Distancing clamp and locator domains..." << Alert::Raise;
    auto distClamp =
      MMG::Distancer(dsfes).setInteriorDomain(Clamp)
                           .distance(dOmega);
    auto distLocator =
      MMG::Distancer(dsfes).setInteriorDomain(Locator)
                           .distance(dOmega);

    Solver::CG cg;

    // Parameters
    const Scalar lambda = 0.5769, mu = 0.3846;

    // VectorFunction f{0, -1, 0}; // Locator
    auto f = -BoundaryNormal(mesh);
    VectorFunction g{1, 0, 0}; // GammaT
    GridFunction nnn(dvfes);
    nnn.projectOnFaces(Average(f));
    nnn.save("nnn.gf");
    nnn.getFiniteElementSpace().getMesh().save("nnn.mesh");

    // Bump function
    auto h = [](Scalar r)
    {
      if (r < -1.0)
        return 1.0;
      else if (r > 1.0)
        return 0.0;
      else
        return 1.0 - 1.0 / (1.0 + std::exp(4 * r / (r * r - 1.0)));
    };

    ScalarFunction heClamp =
      [&](const Geometry::Point& p) { return h(distClamp(p) / epsilon) / epsilon; };

    Alert::Info() << "Solving state equation..." << Alert::Raise;
    TrialFunction u(vfes);
    TestFunction  v(vfes);
    Problem state(u, v);
    state = LinearElasticityIntegral(u, v)(lambda, mu)
          - BoundaryIntegral(g, v).over(GammaT)
          - BoundaryIntegral(f, v).over(Locator)
          + FaceIntegral(heClamp * u, v).over(Clamp)
          + BoundaryIntegral(tgv * u.z(), v.z()).over(7);
    state.assemble();
    Alert::Info() << "Solving state equation..." << Alert::Raise;
    state.solve(cg);

    u.getSolution().save("u.gf");
    mesh.save("u.mesh");

    std::exit(1);

    Alert::Info() << "Solving adjoint equation..." << Alert::Raise;
    TrialFunction p(vfes);
    TestFunction  q(vfes);
    Problem adjoint(p, q);
    adjoint = LinearElasticityIntegral(p, q)(lambda, mu)
            + Integral(u.getSolution() / mesh.getVolume(), q)
            + FaceIntegral(heClamp * p, q).over(Clamp);
    adjoint.solve(cg);


    p.getSolution().save("p.gf");
    mesh.save("p.mesh");

    Alert::Info() << "Computing objective..." << Alert::Raise;
    ScalarGridFunction j(sfes);
    j = Frobenius(u.getSolution()) / mesh.getVolume();
    j.setWeights();
    if (i > 0)
      oldObjective = objective;


    const Scalar J = Integral(j).compute();
    const Scalar pLocator = ellLocator * mesh.getPerimeter(Locator);
    const Scalar pClamp = ellClamp * mesh.getPerimeter(Clamp);
    objective = J + pClamp + pLocator;

    Alert::Info() << "Objective: " << Alert::Notation(objective) << Alert::NewLine
                  << "J: " << J << Alert::NewLine
                  << "pClamp: " << pClamp << Alert::NewLine
                  << "pLocator: " << pLocator
                  << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

    if (topologicalStep)
    {
      Alert::Info() << "Topological optimization..." << Alert::Raise;

      TrialFunction s(dsfes);
      TestFunction  t(dsfes);

      Problem topo(s, t);
      if (clampStep)
      {
        Alert::Info() << "Inserting clamp region..." << Alert::Raise;
        topo = Integral(alpha * Grad(s), Grad(t))
             + Integral(s, t)
             + Integral((u.getSolution().T() * aniso * p.getSolution()).coeff(0, 0), t)
             + tgv * Integral(s, t).over(Clamp, Locator, GammaT);
      }
      else if (locatorStep)
      {
        Alert::Info() << "Inserting locator region..." << Alert::Raise;
        topo = Integral(alpha * Grad(s), Grad(t))
             + Integral(s, t)
             - Integral(Dot(f, p.getSolution()), t)
             + tgv * Integral(s, t).over(Clamp, Locator, GammaT);
      }
      topo.solve(cg);

      s.getSolution().save("Topo.gf");
      dsfes.getMesh().save("Topo.mesh");

      Alert::Info() << "Computing nucleation locations..." << Alert::Raise;
      const Scalar tc = s.getSolution().max();
      std::vector<Point> cs;
      for (auto it = dOmega.getVertex(); !it.end(); ++it)
      {
        const Point p(*it, it->getTransformation(),
            Polytope::getVertex(0, Polytope::Type::Point), it->getCoordinates());
        const Scalar tp = s.getSolution()(p);
        if (tp > 0 && (tp / tc) > (1 - 1e-12))
        {
          cs.emplace_back(std::move(p));
          break;
        }
      }
      Alert::Info() << "Nucleating " << Alert::Notation(cs.size()) << " holes..."
                    << Alert::Raise;

      if (clampStep)
        holes(distClamp, cs);
      else if (locatorStep)
        holes(distLocator, cs);

      distClamp.save("distClamp.gf");
      distClamp.getFiniteElementSpace().getMesh().save("distClamp.mesh");
    }
    else if (geometricStep)
    {
      Alert::Info() << "Geometrical optimization..." << Alert::Raise;

      Alert::Info() << "Computing conormal..." << Alert::Raise;
      GridFunction conormalLocator(dvfes);
      conormalLocator = Grad(distLocator);
      conormalLocator.stableNormalize();

      GridFunction conormalClamp(dvfes);
      conormalClamp = Grad(distClamp);
      conormalClamp.stableNormalize();

      Alert::Info() << "Computing shape gradient..." << Alert::Raise;
      TrialFunction theta(dvfes);
      TestFunction  w(dvfes);
      Problem hilbert(theta, w);
      hilbert = Integral(alpha * Jacobian(theta), Jacobian(w))
              + Integral(theta, w)
              + (1.0 / epsilon) * FaceIntegral(
                Dot(u.getSolution(), p.getSolution()), Dot(conormalClamp, w)).over(dClamp)
              - FaceIntegral(
                  Dot(Average(f), p.getSolution()), Dot(conormalLocator, w)).over(dLocator)
              + ellClamp * FaceIntegral(conormalClamp, w).over(dClamp)
              + ellLocator * FaceIntegral(conormalLocator, w).over(dLocator)
              + tgv * Integral(theta, w).over(GammaT);
      hilbert.solve(cg);

      GridFunction norm(dsfes);
      norm = Frobenius(theta.getSolution());
      theta.getSolution() /= norm.max();

      dOmega.save("Theta.mesh");
      theta.getSolution().save("Theta.gf");

      if (clampStep)
      {
        Alert::Info() << "Advecting clamp..." << Alert::Raise;
        MMG::Advect(distClamp, theta.getSolution()).step(dt);
      }
      else if (locatorStep)
      {
        Alert::Info() << "Advecting locator..." << Alert::Raise;
        MMG::Advect(distLocator, theta.getSolution()).step(dt);
      }
    }


    if (clampStep)
    {
      Alert::Info() << "Meshing the clamp..." << Alert::Raise;
      try
      {
        GridFunction workaround(sfes);
        workaround.projectOnBoundary(distClamp);
        mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                          .split(Clamp, { Clamp, Gamma })
                                          .split(Gamma, { Clamp, Gamma })
                                          .noSplit(Locator)
                                          .noSplit(GammaT)
                                          .setHMax(hmax)
                                          .setHMin(hmin)
                                          .setGradation(hgrad)
                                          .setHausdorff(hausd)
                                          .surface()
                                          .discretize(workaround);
        hmax = 1.1 * hmax > 0.4 ? 0.4 : hmax * 1.1;
      }
      catch (Alert::Exception& e)
      {
        Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
        hmax = 0.9 * hmax < 0.1 ? 0.1 : hmax * 0.9;
        continue;
      }
    }
    else if (locatorStep)
    {
      Alert::Info() << "Meshing the locator..." << Alert::Raise;
      try
      {
        GridFunction workaround(sfes);
        workaround.projectOnBoundary(distLocator);
        mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                          .split(Locator, { Locator, Gamma })
                                          .split(Gamma, { Locator, Gamma })
                                          .noSplit(Clamp)
                                          .noSplit(GammaT)
                                          .setHMax(hmax)
                                          .setHMin(hmin)
                                          .setGradation(hgrad)
                                          .setHausdorff(hausd)
                                          .surface()
                                          .discretize(workaround);
        hmax = 1.1 * hmax > 0.4 ? 0.4 : hmax * 1.1;
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

size_t rmc(MeshBase& mesh)
{
  const size_t per = mesh.getPerimeter();
  const size_t D = mesh.getDimension();
  auto ccl = mesh.ccl(
      [](const Polytope& p1, const Polytope& p2)
      {
        return p1.getAttribute() == p2.getAttribute();
      }, D - 1, { Clamp, Locator } );
  size_t ccs = ccl.getCount();

  for (const auto& cc : ccl)
  {
    Scalar area = 0;
    for (const Index i : cc)
      area += mesh.getFace(i)->getMeasure();
    if ((area / per) < 1e-5)
    {
      for (const Index i : cc)
      {
        if (mesh.getFace(i)->getAttribute() == Clamp
            || mesh.getFace(i)->getAttribute() == Locator)
          mesh.setAttribute({ D - 1, i }, Gamma);
      }
      ccs--;
    }
  }
  return ccs;
}

void holes(ScalarGridFunction& dist, const std::vector<Point>& cs)
{
  if (cs.size())
  {
    auto insert =
      [&](const Point& v)
      {
        Scalar d = dist(v);
        for (const auto& c : cs)
        {
          const Scalar dd = (v - c).norm() - radius;
          d = std::min(d, dd);
        }
        return d;
      };
    dist = insert;
  }
}

