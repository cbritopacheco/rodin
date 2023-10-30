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

static constexpr double hmax = 0.05;
static constexpr double alpha = 0.6;
static constexpr double epsilon = 0.01;
static constexpr double ell = 0.05;
static constexpr double tgv = std::numeric_limits<double>::max();

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
  const char* meshFile =
    "u.mesh";

  // Load and build finite element spaces on the volumetric domain
  MMG::Mesh Omega;
  Omega.load(meshFile);
  //MMG::Optimizer().setHMin(0.1).optimize(Omega);

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
    dOmega.getConnectivity().compute(1, 2);
    dOmega.save("miaow.mesh", IO::FileFormat::MEDIT);
    dOmega.trace({{{GammaD, Gamma}, SigmaD}, {{GammaN, Gamma}, SigmaN}});
    // dOmega.save("miaow.mesh", IO::FileFormat::MEDIT);

    // Build finite element spaces
    Alert::Info() << "   | Building finite element spaces." << Alert::Raise;
    ScalarFES sfes(Omega);
    VectorFES vfes(Omega, Omega.getSpaceDimension());

    ScalarFES dsfes(dOmega);
    VectorFES dvfes(dOmega, dOmega.getSpaceDimension());

    Alert::Info() << "   | Distancing domain." << Alert::Raise;
    auto dist = MMG::Distancer(dsfes).setInteriorDomain(GammaD)
                                     .distance(dOmega);

    Solver::CG solver;

    auto h = [](double r)
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
          + BoundaryIntegral(he * u, v).over({Gamma, GammaD})
          - Integral(f, v)
          - BoundaryIntegral(g, v).over(GammaN);
    state.solve(solver);

    // Adjoint equation
    auto dj = -u.getSolution() / Omega.getVolume();
    Alert::Info() << "   | Solving adjoint equation." << Alert::Raise;
    TrialFunction p(sfes);
    TestFunction  q(sfes);
    Problem adjoint(p, q);
    adjoint = Integral(Grad(p), Grad(q))
            + BoundaryIntegral(he * p, q).over({Gamma, GammaD})
            - Integral(dj, q);
    adjoint.solve(solver);

    double objective = J(u.getSolution());
    Alert::Info() << "   | Objective: " << objective
                  << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

    u.getSolution().save("u.gf");
    Omega.save("u.mesh");

    // Compute the shape gradient
    Alert::Info() << "   | Computing shape gradient." << Alert::Raise;
    auto hadamard = 1. / (epsilon * epsilon) * u.getSolution() * p.getSolution() + ell;
    auto grad = getShapeGradient(dvfes, dist, hadamard, solver);
    grad *= -1.0;

    // Advect the distance function with the gradient
    Alert::Info() << "   | Advecting the distance function." << Alert::Raise;
    GridFunction norm(sfes);
    norm = Frobenius(grad);
    double dt = 2 * hmax / norm.max();
    MMG::Advect(dist, grad).surface().step(dt);

    // // Topological optimization
    // if (i % topoPeriod == 0)
    // {
    //   Alert::Info() << "   | Computing topological sensitivity." << Alert::Raise;

    //   // Compute the topological sensitivity
    //   GridFunction topo(sfes);
    //   topo = M_PI * u.getSolution() * p.getSolution();

    //   // Geodesic distance
    //   auto gd = [&](const Point& x, const Point& c)
    //          {
    //           return std::acos((x(0) * c(0) + x(1) * c(1) + x(2) * c(2)));
    //          };

    //   // double tmin = topo.min();
    //   // double tmax = topo.max();
    //   // if (tmin < 0)
    //   // {
    //   //   // auto cs = topo.where(topo < tmin * (1 - 0.001) + (tmax - tmin) * 0.001);
    //   //   if (cs.size() > 0)
    //   //   {
    //   //    Alert::Info() << "   | " << cs.size() << " possible hole centers." << Alert::Raise;
    //   //    auto it = std::min_element(cs.begin(), cs.end(),
    //   //       [&](const Point& lhs, const Point& rhs) -> bool
    //   //       {
    //   //        return topo(lhs).scalar() < topo(rhs).scalar();
    //   //       });

    //   //    auto holes = ScalarFunction(
    //   //       [&](const Point& v) -> double
    //   //       {
    //   //        double d = dist(v);
    //   //        double dd = gd(v, *it) - radius;
    //   //        d = std::min(d, dd);
    //   //        return d;
    //   //       });
    //   //    dist = holes;
    //   //   }
    //   // }
    // }

   // Mesh only the surface part
   Alert::Info() << "   | Meshing the domain." << Alert::Raise;
   Omega = MMG::ImplicitDomainMesher().noSplit(GammaN)
                                      .split(GammaD, {GammaD, Gamma})
                                      .split(Gamma, {GammaD, Gamma})
                                      .setHMax(hmax)
                                      .surface()
                                      .discretize(dist);

   Alert::Info() << "   | Optimizing the domain." << Alert::Raise;
   MMG::Optimizer().setHMax(hmax).optimize(Omega);

   dOmega.save("out/dOmega." + std::to_string(i) +  ".mesh", IO::FileFormat::MEDIT);
   Omega.save("Omega.mesh", IO::FileFormat::MEDIT);
  }

  return 0;
}

