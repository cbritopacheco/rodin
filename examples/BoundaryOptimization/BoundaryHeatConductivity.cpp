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
static constexpr Geometry::Attribute Gamma = 6;
static constexpr Geometry::Attribute GammaD = 3;
static constexpr Geometry::Attribute GammaN = 2;

static constexpr Geometry::Attribute SigmaD = 3;
static constexpr Geometry::Attribute SigmaN = 2;

static constexpr size_t maxIt = 10000;

static constexpr Scalar epsilon = 0.001;
static constexpr Scalar ell = 1;
static constexpr Scalar radius = 0.02;
static constexpr Scalar tgv = std::numeric_limits<float>::max();

using ScalarFES = P1<Scalar, Context::Serial>;
using VectorFES = P1<Math::Vector, Context::Serial>;
using ScalarGridFunction = GridFunction<ScalarFES>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

inline
size_t rmc(MeshBase& mesh)
{
  const size_t per = mesh.getPerimeter();
  const size_t D = mesh.getDimension();
  auto ccl = mesh.ccl(
      [](const Polytope& p1, const Polytope& p2)
      {
        return p1.getAttribute() == p2.getAttribute();
      }, D - 1, GammaD);
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
        if (mesh.getFace(i)->getAttribute() == GammaD)
          mesh.setAttribute({ D - 1, i }, Gamma);
      }
      ccs--;
    }
  }
  return ccs;
}

int main(int, char**)
{
  const char* meshFile = "linkrods.mesh";
  //const char* meshFile = "Omega.mesh";

  // Load and build finite element spaces on the volumetric domain
  Scalar dc = 1;
  Scalar hmax = 0.05;
  size_t regionCount;
  MMG::Mesh Omega;
  Omega.load(meshFile, IO::FileFormat::MEDIT);

  {
    P1 vh(Omega);
    GridFunction dist(vh);
    dist = [&](const Point& p)
      {
        Math::Vector c1(3);
        c1(0) = 4.3 - p.x();
        c1(1) = 3.8 - p.y();
        c1(2) = 1.0 - p.z();
        double dd = c1.norm() - 0.1;
        return dd;
      };

    Omega = MMG::ImplicitDomainMesher().noSplit(GammaN)
                                       .setAngleDetection(false)
                                       .split(GammaD, {GammaD, Gamma})
                                       .split(Gamma, {GammaD, Gamma})
                                       .setHMax(hmax)
                                       .setHMin(hmax / 5.0)
                                       .setHausdorff(hmax / 10.0)
                                       .setGradation(1.2)
                                       .surface()
                                       .discretize(dist);
  }

  auto J = [&](const ScalarGridFunction& u)
  {
    return Integral(u).compute() + ell * Omega.getPerimeter(GammaD);
  };

  std::ofstream fObj("obj.txt");
  size_t i = 0;
  size_t prevRegionCount = 0;
  while (i < maxIt)
  {
    const Scalar hmin = hmax / 5.0;
    const Scalar hausd = hmax / 10.0;
    const Scalar hgrad = 1.2;
    const Scalar k = 0.5 * (hmax + hmin);
    const Scalar dt = dc * k;

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
                      .setGradation(hgrad)
                      .setHausdorff(hausd)
                      .setAngleDetection(false)
                      .optimize(Omega);
    }
    catch (Alert::Exception& e)
    {
      Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
      hmax = 0.9 * hmax < 0.02 ? 0.02 : hmax * 0.9;
      continue;
    }

    Alert::Info() << "Computing required connectivity..." << Alert::Raise;
    Omega.getConnectivity().compute(2, 3);
    Omega.getConnectivity().compute(2, 2);
    Omega.getConnectivity().compute(1, 2);

    Alert::Info() << "RMC..." << Alert::Raise;
    prevRegionCount = regionCount;
    regionCount = rmc(Omega);

    Alert::Info() << "Found " << Alert::Notation(regionCount) << " regions."
      << Alert::Raise;

    Alert::Info() << "Skinning mesh..." << Alert::Raise;
    auto dOmega = Omega.skin();
    dOmega.trace({{{GammaD, Gamma}, SigmaD}, {{GammaN, Gamma}, SigmaN}});

    Alert::Info() << "Building finite element spaces..." << Alert::Raise;
    ScalarFES sfes(Omega);
    VectorFES vfes(Omega, Omega.getSpaceDimension());
    ScalarFES dsfes(dOmega);
    VectorFES dvfes(dOmega, dOmega.getSpaceDimension());

    Alert::Info() << "Distancing domain..." << Alert::Raise;
    auto dist = MMG::Distancer(dsfes).setInteriorDomain(GammaD)
                                     .distance(dOmega);

    // Parameters
    Solver::CG cg;

    ScalarFunction f = 1;
    ScalarFunction g = -1.0;

    auto h = [](Scalar r)
    {
      if (r < -1.0)
        return 1.0;
      else if (r > 1.0)
        return 0.0;
      else
        return 1.0 - 1.0 / (1.0 + std::exp(4 * r / (r * r - 1.0)));
    };

    ScalarFunction he =
      [&](const Geometry::Point& p) { return h(dist(p) / epsilon) / epsilon; };

    Alert::Info() << "Solving state equation..." << Alert::Raise;
    TrialFunction u(sfes);
    TestFunction  v(sfes);
    Problem state(u, v);
    state = Integral(Grad(u), Grad(v))
          + FaceIntegral(he * u, v).over({Gamma, GammaD})
          - Integral(f, v);
    state.solve(cg);
    u.getSolution().save("u.gf");
    Omega.save("miaow.mesh");

    auto dj = -1.0 / Omega.getVolume();
    Alert::Info() << "Solving adjoint equation..." << Alert::Raise;
    TrialFunction p(sfes);
    TestFunction  q(sfes);
    Problem adjoint(p, q);
    adjoint = Integral(Grad(p), Grad(q))
            + BoundaryIntegral(he * p, q).over({Gamma, GammaD})
            - Integral(dj * q);
    adjoint.solve(cg);

    Alert::Info() << "Computing objective..." << Alert::Raise;
    const Scalar objective = J(u.getSolution());
    Alert::Info() << "Objective: " << Alert::Notation(objective) << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

    Alert::Info() << "Computing conormal to GammaD..." << Alert::Raise;
    GridFunction conormal(dvfes);
    conormal = Grad(dist);
    conormal.stableNormalize();

    Alert::Info() << "Computing shape gradient..." << Alert::Raise;
    auto hadamard = 1. / (epsilon * epsilon) * u.getSolution() * p.getSolution() + ell;
    TrialFunction theta(dsfes);
    TestFunction  w(dsfes);
    Problem hilbert(theta, w);
    hilbert = Integral(dc * dc * Grad(theta), Grad(w))
            + Integral(theta, w)
            + Integral(tgv * g, w).over(GammaN)
            - FaceIntegral(hadamard, w).over(SigmaD);
    hilbert.solve(cg);
    GridFunction grad(dvfes);
    grad = theta.getSolution() * conormal;
    grad *= -1.0;
    GridFunction norm(dsfes);
    norm = Frobenius(grad);
    grad /= norm.max();

    Alert::Info() << "Advecting the distance function." << Alert::Raise;
    MMG::Advect(dist, grad).step(dt);

    if (true)
    {
      Alert::Info() << "Computing topological sensitivity..." << Alert::Raise;
      GridFunction topo(dsfes);
      topo = u.getSolution() * p.getSolution();

      Alert::Info() << "Computing nucleation locations..." << Alert::Raise;
      const Scalar tc = topo.min();
      std::vector<Point> cs;
      for (auto it = dOmega.getVertex(); !it.end(); ++it)
      {
        const Point p(*it, it->getTransformation(),
            Polytope::getVertices(Polytope::Type::Point).col(0), it->getCoordinates());
        const Scalar tp = topo(p);
        if (Math::abs(1 - tc / tp) < 1e-5)
          cs.emplace_back(std::move(p));
      }

      Alert::Info() << "Nucleating " << Alert::Notation(cs.size()) << " holes..."
                    << Alert::Raise;
      if (cs.size())
      {
        auto holes =
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
        dist = holes;
      }
    }

    Alert::Info() << "Meshing the domain..." << Alert::Raise;
    GridFunction workaround(sfes);
    workaround.projectOnBoundary(dist);
    try
    {
      Omega = MMG::ImplicitDomainMesher().noSplit(GammaN)
                                         .split(GammaD, {GammaD, Gamma})
                                         .split(Gamma, {GammaD, Gamma})
                                         .setHMax(hmax)
                                         .setHMin(hmin)
                                         .setHausdorff(hausd)
                                         .setGradation(hgrad)
                                         .setAngleDetection(false)
                                         .surface()
                                         .discretize(workaround);
      hmax = 1.1 * hmax > 0.05 ? 0.05 : hmax * 1.1;
    }
    catch (Alert::Exception& e)
    {
      Alert::Warning() << "Meshing failed. Trying with new parameters." << Alert::Raise;
      hmax = 0.9 * hmax < 0.02 ? 0.02 : hmax * 0.9;
      continue;
    }

    Alert::Info() << "Saving files..." << Alert::Raise;
    Omega.save("Omega.mesh", IO::FileFormat::MEDIT);
    dOmega.save("out/dOmega." + std::to_string(i) +  ".mesh", IO::FileFormat::MEDIT);
    dOmega.save("out/dOmega.mfem." + std::to_string(i) +  ".mesh", IO::FileFormat::MFEM);

    Alert::Success() << "Completed Iteration: " << i << '\n' << Alert::Raise;
    i++;
  }

  return 0;
}

