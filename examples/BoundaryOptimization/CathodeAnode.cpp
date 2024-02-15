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

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// Parameters
static constexpr Geometry::Attribute Anode = 2; // u = 1
static constexpr Geometry::Attribute Cathode = 3; // u = 0
static constexpr Geometry::Attribute Gamma = 4; // dn u = 0

static constexpr Geometry::Attribute dAnode = 3;
static constexpr Geometry::Attribute dCathode = 2;

static constexpr size_t maxIt = 10000;

static constexpr Scalar epsilon = 0.001;
static constexpr Scalar ellAnode = 1;
static constexpr Scalar ellCathode = 1;
static constexpr Scalar radius = 0.02;
static constexpr Scalar tgv = std::numeric_limits<float>::max();
static Scalar hmax = 1;
static const Scalar alpha = 4;

using ScalarFES = P1<Scalar, Context::Sequential>;
using VectorFES = P1<Math::Vector, Context::Sequential>;
using ScalarGridFunction = GridFunction<ScalarFES>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

size_t rmc(MeshBase& mesh);

int main(int, char**)
{
  MMG::Mesh mesh;
  mesh.load("Unnamed1-Fusion001.o.mesh", IO::FileFormat::MEDIT);


  Alert::Info() << "Initializing anode region..." << Alert::Raise;
  {
    P1 vh(mesh);

    GridFunction distAnode(vh);
    distAnode = [&](const Point& p)
    {
      Math::SpatialVector c(3);
      c << 2.5, 35, 7;
      return (p - c).norm() - 1;
    };

    mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                      .split(Anode, { Anode, Gamma })
                                      .split(Gamma, { Anode, Gamma })
                                      .noSplit(Cathode)
                                      .setHMax(hmax)
                                      .surface()
                                      .discretize(distAnode);
  }

  mesh.save("hole.mesh", IO::FileFormat::MEDIT);

  Alert::Info() << "Initializing cathode region..." << Alert::Raise;
  {
    P1 vh(mesh);

    GridFunction distCathode(vh);
    distCathode = [&](const Point& p)
    {
      Math::SpatialVector c(3);
      c << 2.5, 15, -2.5;
      return (p - c).norm() - 1;
    };

    mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                      .split(Cathode, { Cathode, Gamma })
                                      .split(Gamma, { Cathode, Gamma })
                                      .noSplit(Anode)
                                      .setHMax(hmax)
                                      .surface()
                                      .discretize(distCathode);

    mesh.save("hole.mesh", IO::FileFormat::MEDIT);
  }

  std::ofstream fObj("obj.txt");
  size_t i = 0;
  size_t regionCount, prevRegionCount;
  while (i < maxIt)
  {
    const Scalar hmin = hmax / 5.0;
    const Scalar hausd = hmax / 10.0;
    const Scalar hgrad = 1.2;
    const Scalar k = 0.5 * (hmax + hmin);
    const Scalar dt = k;

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
                      .setAngleDetection(false)
                      .optimize(mesh);
    }
    catch (Alert::Exception& e)
    {
      Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
      hmax = 0.9 * hmax < 0.02 ? 0.02 : hmax * 0.9;
      continue;
    }

    Alert::Info() << "Computing required connectivity..." << Alert::Raise;
    mesh.getConnectivity().compute(2, 3); // Computes boundary
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(1, 2);

    Alert::Info() << "RMC..." << Alert::Raise;
    prevRegionCount = regionCount;
    regionCount = rmc(mesh);

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
      MMG::Distancer(dsfes).setInteriorDomain(Anode)
                           .distance(dOmega);
    auto distCathode =
      MMG::Distancer(dsfes).setInteriorDomain(Cathode)
                           .distance(dOmega);

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
          - BoundaryIntegral(heAnode * uIn, v);
    Alert::Info() << "Solving state equation..." << Alert::Raise;
    state.solve(cg);
    u.getSolution().save("u.gf");

    Alert::Info() << "Solving adjoint equation..." << Alert::Raise;
    TrialFunction p(sfes);
    TestFunction  q(sfes);
    Problem adjoint(p, q);
    adjoint = Integral(gamma * Grad(p), Grad(q))
            + BoundaryIntegral((heAnode + heCathode) * p, q)
            - Integral(0.5 * Grad(u.getSolution()), Grad(q));
    adjoint.solve(cg);
    p.getSolution().save("p.gf");
    mesh.save("miaow.mesh");

    Alert::Info() << "Computing objective..." << Alert::Raise;
    GridFunction norm(sfes);
    norm = 0.5 * Frobenius(Grad(u.getSolution()));
    norm.setWeights();
    const Scalar objective =
      Integral(norm).compute() + ellAnode * mesh.getPerimeter(Anode) + ellCathode * mesh.getPerimeter(Cathode);
    Alert::Info() << "Objective: " << Alert::Notation(objective) << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

    Alert::Info() << "Computing conormal..." << Alert::Raise;
    GridFunction conormalCathode(dvfes);
    conormalCathode = Grad(distCathode);
    conormalCathode.stableNormalize();

    GridFunction conormalAnode(dvfes);
    conormalAnode = Grad(distAnode);
    conormalAnode.stableNormalize();

    TrialFunction theta(dvfes);
    TestFunction  w(dvfes);
    Problem hilbert(theta, w);
    dOmega.save("theta.mesh");
    if (i % 2 == 0)
    {
      Alert::Info() << "Computing shape gradient for the anode..." << Alert::Raise;
      hilbert = Integral(alpha * alpha * Jacobian(theta), Jacobian(w))
              + Integral(theta, w)
              - FaceIntegral(
                  1.0 / epsilon * u.getSolution() * p.getSolution() + ellAnode,
                  Dot(conormalAnode, w)).over(dAnode)
              - 1.0 / (alpha * alpha) * FaceIntegral(conormalCathode, w)
              ;
      hilbert.solve(cg);
      theta.getSolution().save("thetaAnode.gf");

      GridFunction norm(dsfes);
      norm = Frobenius(theta.getSolution());
      theta.getSolution() /= norm.max();

      MMG::Advect(distAnode, theta.getSolution()).step(dt);

      Alert::Info() << "Meshing the anode..." << Alert::Raise;
      try
      {
        GridFunction workaround(sfes);
        workaround.projectOnBoundary(distAnode);
        mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                          .split(Anode, { Anode, Gamma })
                                          .split(Gamma, { Anode, Gamma })
                                          .noSplit(Cathode)
                                          .setHMax(hmax)
                                          .surface()
                                          .discretize(distCathode);
        hmax = 1.1 * hmax > 0.05 ? 0.05 : hmax * 1.1;
      }
      catch (Alert::Exception& e)
      {
        Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
        hmax = 0.9 * hmax < 0.02 ? 0.02 : hmax * 0.9;
        continue;
      }
    }
    else
    {
      Alert::Info() << "Computing shape gradient cathode..." << Alert::Raise;
      hilbert = Integral(alpha * alpha * Jacobian(theta), Jacobian(w))
              + Integral(theta, w)
              + FaceIntegral(
                  1.0 / epsilon * (uIn - u.getSolution()) * p.getSolution() + ellCathode,
                  Dot(conormalCathode, w)).over(dCathode)
              - 1.0 / (alpha * alpha) * FaceIntegral(conormalAnode, w)
              ;
      hilbert.solve(cg);
      theta.getSolution().save("thetaCathode.gf");

      GridFunction norm(dsfes);
      norm = Frobenius(theta.getSolution());
      theta.getSolution() /= norm.max();

      MMG::Advect(distCathode, theta.getSolution()).step(dt);

      Alert::Info() << "Meshing the cathode..." << Alert::Raise;
      try
      {
        GridFunction workaround(sfes);
        workaround.projectOnBoundary(distCathode);
        mesh = MMG::ImplicitDomainMesher().setAngleDetection(false)
                                          .split(Cathode, { Cathode, Gamma })
                                          .split(Gamma, { Cathode, Gamma })
                                          .noSplit(Anode)
                                          .setHMax(hmax)
                                          .surface()
                                          .discretize(distCathode);
        hmax = 1.1 * hmax > 0.05 ? 0.05 : hmax * 1.1;
      }
      catch (Alert::Exception& e)
      {
        Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
        hmax = 0.9 * hmax < 0.02 ? 0.02 : hmax * 0.9;
        continue;
      }
    }
  }

  // Alert::Info() << "Saving files..." << Alert::Raise;
  // mesh.save("Omega.mesh", IO::FileFormat::MEDIT);
  // dOmega.save("out/dOmega." + std::to_string(i) +  ".mesh", IO::FileFormat::MEDIT);
  // dOmega.save("out/dOmega.mfem." + std::to_string(i) +  ".mesh", IO::FileFormat::MFEM);

  // Alert::Success() << "Completed Iteration: " << i << '\n' << Alert::Raise;
  i++;

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
      }, D - 1, { Anode, Cathode } );
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
        if (mesh.getFace(i)->getAttribute() == Cathode
            || mesh.getFace(i)->getAttribute() == Anode)
          mesh.setAttribute({ D - 1, i }, Gamma);
      }
      ccs--;
    }
  }
  return ccs;
}


