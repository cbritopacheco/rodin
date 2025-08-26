/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Math.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Models/Eikonal/FMM.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  // Create a simple triangular mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 32, 32 });
  mesh.scale(1.0 / 31.0);
  mesh.getConnectivity().compute(2, 0);
  mesh.getConnectivity().compute(0, 0);

  // Create P1 finite element space
  P1 vh(mesh);
  GridFunction u(vh);

  std::cout << u.getSize() << std::endl;
  std::cout << u.getData().size() << std::endl;

  // Define speed function (constant speed = 1)
  auto f =
    [](const Math::SpatialVector<Real>& p) -> Real
    {
      return 1.0;
    };

  // Create FMM solver
  Models::Eikonal::FMM fmm(u, f);

  // Set interface: center point as source
  std::vector<Index> interface;
  for (auto it = mesh.getVertex(); !it.end(); ++it)
  {
    const auto& coord = mesh.getVertexCoordinates(it->getIndex());
    Real distance = (coord - Math::SpatialVector<Real>{{0.5, 0.5}}).norm();
    if (distance < 0.05)
      interface.push_back(it->getIndex());
  }

  fmm.setInterface(std::move(interface));

  // Solve the Eikonal equation
  fmm.solve();

  // Save results
  u.save("eikonal_distance.gf");
  mesh.save("eikonal_mesh.mesh");

  std::cout << "Fast Marching Method completed successfully!" << std::endl;
  std::cout << "Results saved to eikonal_distance.sol and eikonal_mesh.mesh" << std::endl;

  return 0;
}
