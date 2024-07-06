/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Math.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Models/Distance/Poisson.h>
#include <Rodin/Models/Distance/SignedPoisson.h>
#include <Rodin/Models/Distance/SpaldingTucker.h>
#include <Rodin/Models/Distance/Rvachev.h>

#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static constexpr Attribute boundary = 2;

int main(int, char**)
{
  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 32, 32 });
  mesh.scale(1. / (31));
  mesh.getConnectivity().compute(2, 0);

  // Distance
  {
    P1 vh(mesh);
    GridFunction dist(vh);
    dist = [&](const Point& p)
    {
      Real d = (p - Math::SpatialVector<Real>{{0.75, 0.25}}).norm() - 0.05;
      d = std::min(d, (p - Math::SpatialVector<Real>{{0.25, 0.25}}).norm() - 0.25);
      d = std::min(d, (p - Math::SpatialVector<Real>{{0.75, 0.75}}).norm() - 0.1);
      return d;
    };

    mesh = MMG::ImplicitDomainMesher().setBoundaryReference(5)
      .setAngleDetection().setHMax(0.02).discretize(dist);
    mesh.save("implicit.mesh", IO::FileFormat::MEDIT);
  }

  std::cout << "implicit done\n";
  mesh.getConnectivity().compute(1, 2);

  P1 vh(mesh);
  auto ls = Models::Distance::SignedPoisson()(5, 3, vh);
  auto dist = Models::Distance::SpaldingTucker()(ls);

  dist.save("miaow.gf");
  mesh.save("miaow.mesh");

  return 0;
}
