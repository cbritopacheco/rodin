/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Geometry/Types.h"
#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

static std::set<Attribute> s_epicardium = { 1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 17, } ;
static std::set<Attribute> s_endocardium = { 3, 7 } ;

template <class MeshT>
bool isAnnulus(const MeshT& mesh, Index edgeIdx)
{
  const auto& inc = mesh.getConnectivity().getIncidence({ 1, 2 }, edgeIdx);

  bool boundary = false;
  bool endocardium = false;
  bool different = false;

  for (const auto& j : inc)
  {
    if (!mesh.isBoundary(j))
      break;                 // keep your "must be boundary" convention
    boundary = true;

    const auto attr = mesh.getAttribute(2, j); // Optional<Attribute>
    if (!attr)
      continue; // no attribute => cannot be in endocardium set, also cannot mark different

    if (s_endocardium.contains(*attr))
      endocardium = true;
    else
      different = true;
  }

  return boundary && endocardium && different;
}

int main(int, char**)
{
  Mesh mesh;
  mesh.load("BiEllipsoid.mesh", IO::FileFormat::MEDIT);
  mesh.scale(10);
  mesh.save("BiEllipsoid_scaled.mesh", IO::FileFormat::MEDIT);
  std::exit(1);

  mesh.getConnectivity().compute(2, 3);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(1, 1);

  auto ccl = mesh.ccl(
      1,
      [&](const Polytope& a, const Polytope& b)
      {
        return isAnnulus(mesh, a.getIndex()) && isAnnulus(mesh, b.getIndex());
      },
      [&](const Polytope& p)
      {
      return isAnnulus(mesh, p.getIndex());
      });

  std::cout << "Found " << ccl.getCount() << " annuli." << std::endl;

  Attribute label = 50;
  for (const auto& cc : ccl)
  {
    for (const auto& i : cc)
    {
      if (!s_epicardium.contains(i))
        mesh.setAttribute({ 1, i }, label);
    }
    label++;
  }

  mesh.save("Annulus.mesh", IO::FileFormat::MEDIT);

  return 0;
}

