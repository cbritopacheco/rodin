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
static constexpr Real epsilon = 0.001;
static constexpr Real ell = 1;
static constexpr Real hmax = 0.2 / 5.0;
static constexpr Real hmin = hmax / 5.0;
static constexpr Real hausd = hmax / 10.0;
static constexpr Real hgrad = 1.2;
static constexpr Real dx = 0.5 * (hmax + hmin);
static constexpr Real cfl = 0.1;
static constexpr Real dt = cfl * dx;
static constexpr Real tgv = std::numeric_limits<float>::max();

using RealFES = P1<Real>;
using VectorFES = P1<Math::Vector<Real>>;
using RealGridFunction = GridFunction<RealFES>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

Real J(const RealGridFunction& u, Real ell);

Real cutoff(Real r);

int main(int, char**)
{
  const char* meshFile = "linkrods.mesh";

  // Load and build finite element spaces on the volumetric domain
  MMG::Mesh mesh;
  mesh.load(meshFile, IO::FileFormat::MEDIT);

  Alert::Info() << "Loaded mesh." << Alert::Raise;

  {
    P1 vh(mesh);
    GridFunction dist(vh);
    dist = [&](const Point& p)
      {
        Math::Vector<Real> c1(3);
        c1(0) = 4.3 - p.x();
        c1(1) = 3.8 - p.y();
        c1(2) = 1.0 - p.z();
        double dd = c1.norm() - 0.1;
        return dd;
      };

    mesh = MMG::ImplicitDomainMesher().noSplit(GammaN)
                                      .setAngleDetection(false)
                                      .split(GammaD, {GammaD, Gamma})
                                      .split(Gamma, {GammaD, Gamma})
                                      .setHMax(hmax)
                                      .setHMin(hmax / 5.0)
                                      // .setHausdorff(hmax / 10.0)
                                      // .setGradation(1.2)
                                      .surface()
                                      .discretize(dist);
  }

  Alert::Info() << "Discretized mesh." << Alert::Raise;
  mesh.save("Omega0.mesh", IO::FileFormat::MEDIT);

  size_t i = 0;
  while (i < maxIt)
  {
    Alert::Info() << "Iteration: " << i                         << Alert::NewLine
                  << "HMax:      " << Alert::Notation(hmax)     << Alert::NewLine
                  << "HMin:      " << Alert::Notation(hmin)     << Alert::NewLine
                  << "Hausdorff: " << Alert::Notation(hausd)    << Alert::NewLine
                  << "HGrad:     " << Alert::Notation(hgrad)    << Alert::NewLine
                  << "dt:        " << Alert::Notation(dt)       << Alert::Raise;

    Alert::Info() << "Optimizing the domain..." << Alert::Raise;
    MMG::Optimizer().setHMax(hmax)
                    .setHMin(hmin)
                    .setGradation(hgrad)
                    .setHausdorff(hausd)
                    .setAngleDetection(false)
                    .optimize(mesh);

    mesh.save("Omega0.mesh", IO::FileFormat::MEDIT);
    std::exit(1);

    Alert::Info() << "Computing required connectivity..." << Alert::Raise;
    mesh.getConnectivity().compute(2, 3);
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(1, 2);

    Alert::Info() << "Skinning mesh..." << Alert::Raise;
    auto dOmega = mesh.skin();
    dOmega.trace({{{GammaD, Gamma}, SigmaD}, {{GammaN, Gamma}, SigmaN}});

    Alert::Info() << "Building finite element spaces..." << Alert::Raise;

    RealFES sfes(mesh);
    VectorFES vfes(mesh, mesh.getSpaceDimension());

    RealFES dsfes(dOmega);
    VectorFES dvfes(dOmega, dOmega.getSpaceDimension());

    Alert::Info() << "Distancing domain..." << Alert::Raise;
    auto dist = MMG::Distancer(dsfes).setInteriorDomain(GammaD)
                                     .distance(dOmega);

    RealFunction f = 1;
    RealFunction g = -1.0;

    RealFunction he =
      [&](const Geometry::Point& p) { return cutoff(dist(p) / epsilon) / epsilon; };

    Alert::Info() << "Solving state equation..." << Alert::Raise;

    TrialFunction u(sfes);
    TestFunction  v(sfes);
    Problem state(u, v);
    state = Integral(Grad(u), Grad(v))
          + FaceIntegral(he * u, v).over({Gamma, GammaD})
          - Integral(f, v);
    Solver::CG(state).solve();

    Alert::Info() << "Solving adjoint equation..." << Alert::Raise;

    auto dj = -1.0 / mesh.getVolume();
    TrialFunction p(sfes);
    TestFunction  q(sfes);
    Problem adjoint(p, q);
    adjoint = Integral(Grad(p), Grad(q))
            + BoundaryIntegral(he * p, q).over({Gamma, GammaD})
            - Integral(dj * q);
    Solver::CG(adjoint).solve();

    Alert::Info() << "Computing objective..." << Alert::Raise;
    const Real objective = J(u.getSolution(), ell);
    Alert::Info() << "Objective: " << Alert::Notation(objective) << Alert::Raise;
  }

}

Real J(const RealGridFunction& u, Real ell)
{
  const auto& fes = u.getFiniteElementSpace();
  const auto& mesh = fes.getMesh();
  return Integral(u).compute() + ell * mesh.getPerimeter(GammaD);
}

Real cutoff(Real r)
{
  if (r < -1.0)
    return 1.0;
  else if (r > 1.0)
    return 0.0;
  else
    return 1.0 - 1.0 / (1.0 + std::exp(4 * r / (r * r - 1.0)));
}
