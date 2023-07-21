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
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::External;

static constexpr Attribute interior = 1;
static constexpr Attribute exterior = 2;
static constexpr Attribute ball = 3;

static constexpr Scalar hmax = 0.01;
static constexpr Scalar hmin = 0.01 * hmax;

int main(int, char**)
{
  // Build a mesh
  MMG::Mesh mesh;
  mesh = mesh.build().initialize(2)
                     .nodes(5)
                     .vertex({ 0, 0 })
                     .vertex({ 1, 0 })
                     .vertex({ 1, 1 })
                     .vertex({ 0, 1 })
                     .vertex({ 0.5, 0.5 })
                     .polytope(Polytope::Geometry::Triangle, { 0, 1, 4 })
                     .polytope(Polytope::Geometry::Triangle, { 1, 2, 4 })
                     .polytope(Polytope::Geometry::Triangle, { 2, 3, 4 })
                     .polytope(Polytope::Geometry::Triangle, { 3, 0, 4 })
                     .polytope(Polytope::Geometry::Segment, { 0, 4 })
                     .polytope(Polytope::Geometry::Segment, { 1, 4 })
                     .polytope(Polytope::Geometry::Segment, { 0, 3 })
                     .polytope(Polytope::Geometry::Segment, { 1, 2 })
                     .finalize();
  mesh.getConnectivity().compute(1, 0);

  mesh.setCorner(0).setCorner(1).setCorner(2).setCorner(3).setCorner(4)
      .setRidge(2).setRidge(3);
  mesh.setAttribute({2, 0}, 1);
  mesh.setAttribute({2, 1}, 2);
  mesh.setAttribute({2, 2}, 2);
  mesh.setAttribute({2, 3}, 2);

  for (auto it = mesh.getFace(); !it.end(); ++it)
    mesh.setAttribute({ it->getDimension(), it->getIndex() }, 2);

  mesh.setAttribute({1, 0}, 3);
  mesh.setAttribute({1, 1}, 3);
  mesh.setAttribute({1, 2}, 4);
  mesh.setAttribute({1, 3}, 5);

  // MMG::Optimize().setAngleDetection(false).setHMin(hmin).setHMax(hmax).optimize(mesh);

  // // Refine around the interior corner
  // {
  //   P1 fes(mesh);
  //   MMG::ScalarGridFunction sizeMap(fes);
  //   sizeMap = [](const Geometry::Point& p)
  //   {
  //     const Scalar s = std::abs(std::sqrt((p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5)));
  //     if (s == 0)
  //       return hmin;
  //     else
  //       return s;
  //   };
  //   MMG::Adapt().setAngleDetection(false).setHMax(hmax).setHMin(hmin).adapt(mesh, sizeMap);
  //   MMG::Optimize().setAngleDetection(false).setHMin(hmin).setHMax(hmax).optimize(mesh);
  // }

  mesh.save("Q.medit.mesh", IO::FileFormat::MEDIT);
  mesh.save("Q.mfem.mesh", IO::FileFormat::MFEM);

  return 0;
}


