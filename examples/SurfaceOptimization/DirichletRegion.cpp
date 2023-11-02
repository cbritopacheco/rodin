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

static constexpr size_t maxIt = 250;

static constexpr Scalar hmax = 0.04;
static constexpr Scalar hmin = 0.1 * hmax;
static constexpr Scalar hausd = 0.25 * hmin;
static constexpr Scalar alpha = 2.0;
static constexpr Scalar epsilon = hmin;
static constexpr Scalar ell = 0.01;
static constexpr Scalar radius = 1.1 * hmax;
static constexpr Scalar tgv = std::numeric_limits<float>::max();

using ScalarFES = P1<Scalar, Context::Serial>;
using VectorFES = P1<Math::Vector, Context::Serial>;
using ScalarGridFunction = GridFunction<ScalarFES>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

template <class Derived, class Solver>
ShapeGradient getShapeGradient(
  const VectorFES& vecFes, const ScalarGridFunction& dist,
  const FunctionBase<Derived>& expr, Solver& solver)
{
  TrialFunction d(vecFes);
  TestFunction  v(vecFes);

  auto gdist = Grad(dist);
  gdist.traceOf(GammaD);

  Problem conormal(d, v);
  conormal = Integral(alpha * Jacobian(d), Jacobian(v))
           + Integral(d, v)
           - FaceIntegral(gdist, v).over(SigmaD);
  conormal.solve(solver);

  const auto& cnd = d.getSolution();
  const auto cn = cnd / Frobenius(cnd);

  TrialFunction g(vecFes);
  Problem hilbert(g, v);
  hilbert = Integral(alpha * Jacobian(g), Jacobian(v))
          + Integral(g, v)
          + Integral(tgv * g, v).over(GammaN)
          - FaceIntegral(expr * cn, v).over(SigmaD);
  hilbert.solve(solver);

  return g.getSolution();
}

int main(int, char**)
{
  const char* meshFile = "cool.mesh";
  //const char* meshFile = "Omega.mesh";

  // Load and build finite element spaces on the volumetric domain
  MMG::Mesh Omega;
  Omega.load(meshFile, IO::FileFormat::MEDIT);
  // MMG::Optimizer().setHMin(0.1).optimize(Omega);

  auto J = [&](const ScalarGridFunction& u)
  {
    return Integral(u).compute() + ell * Omega.getPerimeter(GammaD);
  };

  std::ofstream fObj("obj.txt");
  for (size_t i = 0; i < maxIt; i++)
  {
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    // Skin the mesh, computing the borders of the new regions
    Alert::Info() << "   | Skinning mesh." << Alert::Raise;
    Omega.getConnectivity().compute(2, 3);
    Omega.getConnectivity().compute(1, 2);
    auto dOmega = Omega.skin();
    dOmega.trace({{{GammaD, Gamma}, SigmaD}, {{GammaN, Gamma}, SigmaN}});
    dOmega.save("dOmega.mesh", IO::FileFormat::MEDIT);

    // Build finite element spaces
    Alert::Info() << "   | Building finite element spaces." << Alert::Raise;
    ScalarFES sfes(Omega);
    VectorFES vfes(Omega, Omega.getSpaceDimension());

    ScalarFES dsfes(dOmega);
    VectorFES dvfes(dOmega, dOmega.getSpaceDimension());

    Alert::Info() << "   | Distancing domain." << Alert::Raise;
    auto dist = MMG::Distancer(dsfes).setInteriorDomain(GammaD)
                                     .distance(dOmega);

    Solver::CG cg;

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

    // State equation
    Alert::Info() << "   | Solving state equation." << Alert::Raise;
    ScalarFunction f = 1;
    ScalarFunction g = -1.0;

    TrialFunction u(sfes);
    TestFunction  v(sfes);
    Problem state(u, v);
    state = Integral(Grad(u), Grad(v))
          + FaceIntegral(he * u, v).over({Gamma, GammaD})
          - Integral(f, v)
          - FaceIntegral(g, v).over(GammaN);
    state.solve(cg);
    u.getSolution().save("u.gf");
    Omega.save("u.mesh");

    // Adjoint equation
    auto dj = -1.0 / Omega.getVolume();
    Alert::Info() << "   | Solving adjoint equation." << Alert::Raise;
    TrialFunction p(sfes);
    TestFunction  q(sfes);
    Problem adjoint(p, q);
    adjoint = Integral(Grad(p), Grad(q))
            + BoundaryIntegral(he * p, q).over({Gamma, GammaD})
            - Integral(dj * q);
    adjoint.solve(cg);
    p.getSolution().save("p.gf");

    const Scalar objective = J(u.getSolution());
    Alert::Info() << "   | Objective: " << objective
                  << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

    // Compute the shape gradient
    Alert::Info() << "   | Computing shape gradient." << Alert::Raise;
    auto hadamard = 1. / (epsilon * epsilon) * u.getSolution() * p.getSolution() + ell;
    auto grad = getShapeGradient(dvfes, dist, hadamard, cg);
    grad *= -1.0;
    grad.save("grad.gf");
    dOmega.save("grad.mesh");

    // Advect the distance function with the gradient
    Alert::Info() << "   | Advecting the distance function." << Alert::Raise;

    GridFunction norm(dsfes);
    norm = Frobenius(grad);

    const Scalar dt = 2 * hmax / norm.max();
    MMG::Advect(dist, grad).surface().step(dt);

    // Topological optimization
    Alert::Info() << "   | Computing topological sensitivity." << Alert::Raise;

    // Compute the topological sensitivity
    GridFunction topo(sfes);
    topo = M_PI * u.getSolution() * p.getSolution();

    Scalar tc = std::numeric_limits<Scalar>::min();
    std::optional<Point> c;
    for (auto it = Omega.getVertex(); !it.end(); ++it)
    {
      const Point p(*it, it->getTransformation(),
          Polytope::getVertices(Polytope::Type::Point).col(0), it->getCoordinates());
      const Scalar tp = topo(p);
      if (tp < tc)
      {
        tc = tp;
        c.emplace(std::move(p));
      }
    }

    auto holes =
      [&](const Point& v)
      {
        Scalar d = dist(v);
        const Scalar dd = (v - c.value()).norm() - radius;
        d = std::min(d, dd);
        return d;
      };
    dist = holes;

    // Mesh only the surface part
    GridFunction workaround(sfes);
    workaround.projectOnBoundary(dist);
    workaround.getFiniteElementSpace().getMesh().save("workaround.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "   | Meshing the domain." << Alert::Raise;
    Omega = MMG::ImplicitDomainMesher().noSplit(GammaN)
                                       .split(GammaD, {GammaD, Gamma})
                                       .split(Gamma, {GammaD, Gamma})
                                       .setHMax(hmax)
                                       .setHMin(hmin)
                                       .setAngleDetection(false)
                                       .surface()
                                       .discretize(workaround);

    Alert::Info() << "   | Optimizing the domain." << Alert::Raise;
    MMG::Optimizer().setHMax(hmax).setHMin(hmin)
      .setAngleDetection(false).optimize(Omega);

    dOmega.save("out/dOmega." + std::to_string(i) +  ".mesh", IO::FileFormat::MEDIT);
    Omega.save("Omega.mesh", IO::FileFormat::MEDIT);
  }

  return 0;
}

