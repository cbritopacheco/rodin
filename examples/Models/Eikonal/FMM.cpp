/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Geometry/GeometryIndexed.h"
#include "Rodin/Geometry/Point.h"
#include <Rodin/Math.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Eikonal/FMM.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  std::cout << "=== Fast Marching Method - Enhanced Version ===" << std::endl;
  std::cout << "Supporting surface meshes and curved elements" << std::endl;

  // Example 1: 2D triangular mesh (surface in 2D space)
  std::cout << "\n--- Example 1: 2D Triangular Mesh ---" << std::endl;
  {
    Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 32, 32 });
    mesh.scale(1.0 / 31.0);
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(0, 0);

    P1 vh(mesh);
    GridFunction u(vh);

    std::cout << "Mesh dimension: " << mesh.getDimension() << std::endl;
    std::cout << "Space dimension: " << mesh.getSpaceDimension() << std::endl;
    std::cout << "Degrees of freedom: " << u.getSize() << std::endl;

    // Define speed function (constant speed = 1)
    auto f = [](const Geometry::Point& p) -> Real { return 1.0; };

    Eikonal::FMM fmm(u, f);

    // Set interface: center point as source
    std::vector<Index> interface;
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      const auto& coord = mesh.getVertexCoordinates(it->getIndex());
      Real distance = (coord - Math::SpatialVector<Real>{{0.5, 0.5}}).norm();
      if (distance < 0.05)
        interface.push_back(it->getIndex());
    }

    fmm.seed(interface).solve();

    u.save("eikonal_2d.gf");
    mesh.save("eikonal_2d.mesh");
    std::cout << "2D example completed successfully!" << std::endl;
  }

  // Example 2: Surface mesh with varying speed function
  std::cout << "\n--- Example 2: Surface Mesh with Varying Speed ---" << std::endl;
  {
    Mesh volume;
    // mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 24, 24 });
    // mesh.scale(2.0 / 23.0);
    // mesh.displace(VectorFunction{ -1.0, -1.0 });
    volume.load("../resources/mmg/MechanicalPiece_1.mesh", IO::FileFormat::MEDIT);

    volume.getConnectivity().compute(2, 3);
    volume.getConnectivity().compute(2, 0);
    volume.getConnectivity().compute(0, 0);
    volume.getConnectivity().compute(0, 2);
    volume.getConnectivity().compute(2, 1); // For edge-based computations
    volume.getConnectivity().compute(0, 1); // For edge-based computations

    Mesh mesh = volume.skin();
    mesh.save("skin.mesh", IO::FileFormat::MEDIT);

    P1 vh(mesh);
    GridFunction u(vh);

    std::cout << "Surface mesh with " << mesh.getVertexCount() << " vertices" << std::endl;

    // Define varying speed function (slower in center, faster at edges)
    auto varying_speed = [](const Geometry::Point& p) -> Real
    {
      Real r = p.norm();
      return 0.5 + 1.5 * std::exp(-2.0 * r * r); // Gaussian speed profile
    };

    Eikonal::FMM fmm(u, varying_speed);

    // Set interface: corner point as source
    std::set<Index> s;
    for (auto it = mesh.getFace(); !it.end(); ++it)
    {
      if (it->getAttribute() == 118)
      {
        for (const auto& v : it->getVertices())
          s.insert(v);
      }
    }

    std::vector<Index> interface;
    interface.reserve(s.size());
    for (const auto& v : s)
      interface.push_back(v);

    std::cout << "Interface size: " << interface.size() << std::endl;

    fmm.seed(interface).solve();

    u.save("eikonal_surface.gf");
    mesh.save("eikonal_surface.mesh");
    std::cout << "Surface example completed successfully!" << std::endl;
  }

  std::cout << "\n=== All examples completed successfully! ===" << std::endl;
  std::cout << "Results saved to:" << std::endl;
  std::cout << "  - eikonal_2d.gf, eikonal_2d.mesh" << std::endl;
  std::cout << "  - eikonal_surface.gf, eikonal_surface.mesh" << std::endl;

  return 0;
}
