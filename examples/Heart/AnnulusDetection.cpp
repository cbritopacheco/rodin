/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

static constexpr Attribute Endocardium = 1;
static constexpr Attribute Annulus = 50;

template <class Mesh>
bool isAnnulus(const Mesh& mesh, Index edgeIdx)
{
  const auto& inc = mesh.getConnectivity().getIncidence({ 1, 2 }, edgeIdx);
  bool boundary = false;
  bool endocardium = false;
  bool different = false;
  for (const auto& j : inc)
  {
    if (mesh.isBoundary(j))
      boundary = true;
    else
      break;
    if (mesh.getAttribute(2, j) == Endocardium)
      endocardium = true;
    if (mesh.getAttribute(2, j) != Endocardium)
      different = true;
  }
  return boundary && endocardium && different;
}

int main(int, char**)
{
  Mesh mesh;
  mesh.load("LeftAtriumTetra.mesh", IO::FileFormat::MEDIT);
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

  size_t label = Annulus;
  for (const auto& cc : ccl)
  {
    for (const auto& i : cc)
    {
      if (label != 1)
        mesh.setAttribute({ 1, i }, label);
    }
    label++;
  }

  mesh.save("Annulus.mesh", IO::FileFormat::MEDIT);

  return 0;
}

