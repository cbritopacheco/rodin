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

static constexpr Geometry::Attribute SigmaD = 1;
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

  Problem conormal(d, v);
  conormal = Integral(alpha * Jacobian(d), Jacobian(v))
        + Integral(d, v)
        - FaceIntegral(Grad(dist).traceOf(GammaD), v).over(SigmaD);
  conormal.solve(solver);

  const auto& cnd = d.getSolution();
  const auto cn = cnd / Frobenius(cnd);

  TrialFunction g(vecFes);
  Problem hilbert(g, v);
  // hilbert = Integral(alpha * Jacobian(g), Jacobian(v))
  //      + Integral(g, v)
  //      + Integral(tgv * g, v).over(GammaN)
  //      - BoundaryIntegral(expr * cn, v).over(SigmaD);
  hilbert.solve(solver);

  return g.getSolution();
}


int main(int, char**)
{
  const char* meshFile = "../resources/mmg/dirichlet-region-example.mesh";

  // Load and build finite element spaces on the volumetric domain
  MMG::Mesh Omega;
  Omega.load(meshFile, IO::FileFormat::MEDIT);
  Omega.asSubMesh();

  auto J = [&](const ScalarGridFunction& u)
  {
    return Integral(u).compute() + ell * Omega.getPerimeter(GammaD);
  };
}
