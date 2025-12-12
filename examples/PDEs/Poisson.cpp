/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <type_traits>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

constexpr auto order = std::integral_constant<size_t, 3>{};

int main(int, char**)
{
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Wedge, { 2, 2, 2 });

  // mesh = mesh.Build().initialize(3).nodes(6).vertex({ 0.0, 0.0, 0.0 })
  //                             .vertex({ 1.0, 0.0, 0.0 })
  //                             .vertex({ 0.0, 1.0, 0.0 })
  //                             .vertex({ 0.0, 0.0, 1.0 })
  //                             .vertex({ 1.0, 0.0, 1.0 })
  //                             .vertex({ 0.0, 1.0, 1.0 })
  //                             .polytope(Polytope::Type::Wedge, { 0, 1, 2, 3, 4, 5 }).finalize();

  mesh.getConnectivity().compute(3, 2);
  mesh.getConnectivity().compute(2, 1);

  mesh.save("RODIN.mesh");

  H1 vh(order, mesh);


  for (size_t c = 0; c < mesh.getPolytopeCount(3); c++)
  {
    std::cout << "Wedge element DOFs:" << std::endl;
    for (const auto& dof : vh.getDOFs(3, 0))
    {
      std::cout << "DOF: " << dof << std::endl;
    }

    const auto& conn32 = mesh.getConnectivity().getIncidence({ 3, 2}, c);
    size_t localFaceIdx = 0;
    for (size_t i = 0; i < conn32.size(); i++)
    {
      std::cout << "Local face idx: " << localFaceIdx << std::endl;
      const auto& dofs = vh.getDOFs(2, conn32[i]);
      for (const auto& dof : dofs)
      {
        std::cout << "DOF: " << dof << std::endl;
      }
      localFaceIdx++;
    }
  }

  GridFunction u(vh);

  u = [](const Geometry::Point& p)
  {
    return p.x() * (1 - p.y()) + p.z();
  };

  u.save("RODIN.gf");

  return 0;
}
