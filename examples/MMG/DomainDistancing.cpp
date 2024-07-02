/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static constexpr Attribute material = RODIN_DEFAULT_POLYTOPE_ATTRIBUTE;
static constexpr Attribute interior = 1;
static constexpr Attribute exterior = 2;
static constexpr Real hmax = 0.1;
static constexpr Real radius = 0.1;

int main(int, char**)
{
  // const size_t n = 16;
  // MMG::Mesh mesh;
  // mesh = mesh.UniformGrid(Polytope::Type::Triangle, n, n);
  // mesh.scale(1. / (n - 1));

  // P1 fes(mesh);
  // MMG::RealGridFunction gf(fes);
  // gf = [](const Geometry::Point& p)
  //      {
  //        return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - radius;
  //      };

  // gf.save("LevelSet.gf");
  // mesh.save("Domain.mesh");

  // MMG::Mesh implicit = MMG::ImplicitDomainMesher().split(material, { interior, exterior })
  //                                                 .setHMax(hmax)
  //                                                 .discretize(gf);
  // P1 miaow(implicit);

  // implicit.save("Discretized.mesh", IO::FileFormat::MEDIT);

  // auto dist = MMG::Distance(miaow).setInteriorDomain(interior).distance(implicit);

  // dist.save("Discretized.sol", IO::FileFormat::MEDIT);

  // MMG::Mesh implicit2 = MMG::ImplicitDomainMesher().split(material, { interior, exterior })
  //                                                  .split(exterior, { interior, exterior })
  //                                                  .setHMax(hmax)
  //                                                  .discretize(dist);

  // MMG::Optimize().setHMin(0.05).setHMax(hmax).optimize(implicit2);

  // implicit2.save("Discretized2.mesh", IO::FileFormat::MEDIT);

  MMG::Mesh mesh;
  mesh.load("degu.mesh", IO::FileFormat::MEDIT);
  mesh.save("implicit.mesh", IO::FileFormat::MFEM);
  P1 fes(mesh);
  // MMG::RealGridFunction gf(fes);

  auto dist = MMG::Distancer(fes).distance(mesh);
  dist.save("implicit.gf");


  return 0;
}
