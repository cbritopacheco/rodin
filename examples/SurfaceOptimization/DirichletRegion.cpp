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
#include <Rodin/Models/Distance/SignedPoisson.h>

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

static constexpr size_t maxIt = 1000;

static constexpr Scalar epsilon = 0.001;
static constexpr Scalar ell = 1;
static constexpr Scalar radius = 0.02;
static constexpr Scalar tgv = std::numeric_limits<float>::max();

using ScalarFES = P1<Scalar, Context::Serial>;
using VectorFES = P1<Math::Vector, Context::Serial>;
using ScalarGridFunction = GridFunction<ScalarFES>;
using VectorGridFunction = GridFunction<VectorFES>;
using ShapeGradient = VectorGridFunction;

template <class Derived, class Solver>
ShapeGradient getShapeGradient(
  const VectorFES& vecFes, const ScalarGridFunction& dist,
  const FunctionBase<Derived>& expr, Solver& solver, Scalar alpha)
{
  TrialFunction d(vecFes);
  TestFunction  v(vecFes);

  auto gdist = Grad(dist);
  gdist.traceOf(GammaD);

  Problem conormal(d, v);
  conormal = Integral(alpha * alpha * Jacobian(d), Jacobian(v))
           + Integral(d, v)
           - FaceIntegral(gdist, v).over(SigmaD);
  conormal.solve(solver);

  const auto& cnd = d.getSolution();
  const auto cn = cnd / Frobenius(cnd);

  TrialFunction g(vecFes);
  Problem hilbert(g, v);
  hilbert = Integral(alpha * alpha * Jacobian(g), Jacobian(v))
          + Integral(g, v)
          + Integral(tgv * g, v).over(GammaN)
          - FaceIntegral(expr * cn, v).over(SigmaD);
  hilbert.solve(solver);
  return g.getSolution();
}

inline
void rmc(MeshBase& mesh)
{
  const size_t per = mesh.getPerimeter();
  const size_t D = mesh.getDimension();
  auto ccl = mesh.ccl(
      [](const Polytope& p1, const Polytope& p2)
      {
        return p1.getAttribute() == p2.getAttribute();
      }, D - 1, GammaD);

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
    }
  }
}

int main(int, char**)
{
  const char* meshFile = "linkrods.mesh";
  //const char* meshFile = "Omega.mesh";

  // Load and build finite element spaces on the volumetric domain
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
                                       .setHMax(0.05)
                                       .setHMin(0.005)
                                       .setHausdorff(0.005)
                                       .setGradation(1.2)
                                       .surface()
                                       .discretize(dist);
  }

  Omega.save("Omega0.mesh", IO::FileFormat::MEDIT);

  auto J = [&](const ScalarGridFunction& u)
  {
    return Integral(u).compute() + ell * Omega.getPerimeter(GammaD);
  };

  std::ofstream fObj("obj.txt");
  Scalar hmax = 0.05 ;
  Scalar dc = 3;
  size_t i = 0;
  while (i < maxIt)
  {
    Scalar hmin = 0.005;
    Scalar hausd = 0.005;
    Scalar hgrad = 1.2;

    Alert::Info() << "----- Iteration: " << i << Alert::NewLine
                  << "hmax: " << hmax
                  << Alert::Raise;

    Alert::Info() << "   | Optimizing the domain." << Alert::Raise;
    MMG::Optimizer().setHMax(hmax)
                    .setHMin(hmin)
                    .setGradation(hgrad)
                    .setHausdorff(hausd)
                    .setAngleDetection(false)
                    .optimize(Omega);

    Omega.save("Omega.mesh", IO::FileFormat::MEDIT);

    // Skin the mesh, computing the borders of the new regions
    Alert::Info() << "   | Skinning mesh." << Alert::Raise;
    Omega.getConnectivity().compute(2, 3);
    Omega.getConnectivity().compute(1, 2);

    Alert::Info() << "RMC." << Alert::Raise;
    Omega.getConnectivity().compute(2, 2);
    rmc(Omega);

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

    // Solver::CG cg;
    Eigen::ConjugateGradient<Math::SparseMatrix, Eigen::Upper | Eigen::Lower> ecg;
    Solver::EigenSolver cg(std::ref(ecg));

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
          - Integral(f, v);
          //- FaceIntegral(g, v).over(GammaN);
    state.solve(cg);

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

    const Scalar objective = J(u.getSolution());
    Alert::Info() << "   | Objective: " << objective
                  << Alert::Raise;
    fObj << objective << "\n";
    fObj.flush();

    // Compute the shape gradient
    Alert::Info() << "   | Computing shape gradient." << Alert::Raise;
    auto hadamard = 1. / (epsilon * epsilon) * u.getSolution() * p.getSolution() + ell;
    auto grad = getShapeGradient(dvfes, dist, hadamard, cg, dc);


    // Advect the distance function with the gradient
    Alert::Info() << "   | Advecting the distance function." << Alert::Raise;
    GridFunction norm(dsfes);
    norm = Frobenius(grad);
    grad *= -1.0;
    grad /= norm.max();

    Scalar k = 0.5 * (hmax + hmin);

    Alert::Info() << "k = " << k << Alert::NewLine << "dc = " << dc << Alert::Raise;
    const Scalar dt = dc * k;
    MMG::Advect(dist, grad).step(dt);

    // Topological optimization
    if (true)
    {
      Alert::Info() << "   | Computing topological sensitivity." << Alert::Raise;
      Scalar tc = 0;
      std::optional<Point> c;
      for (auto it = Omega.getVertex(); !it.end(); ++it)
      {
        const auto topo = u.getSolution() * p.getSolution();
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
    }

    // Mesh only the surface part
    GridFunction workaround(sfes);
    workaround.projectOnBoundary(dist);
    Omega.save("dist.mesh", IO::FileFormat::MEDIT);
    workaround.save("dist.sol", IO::FileFormat::MEDIT);

    Alert::Info() << "   | Meshing the domain." << Alert::Raise;
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

      // hmax = hmax > 0.1 ? 0.1 : hmax * 1.1;
    }
    catch (Alert::Exception& e)
    {
      hmax = hmax < 0.05 ? 0.05 : hmax * 0.8;
      continue;
    }

    dOmega.save("out/dOmega." + std::to_string(i) +  ".mesh", IO::FileFormat::MEDIT);
    dOmega.save("out/dOmega.mfem." + std::to_string(i) +  ".mesh", IO::FileFormat::MFEM);
    i++;
  }

  return 0;
}

