/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

/**
 * Edge is a boundary-interface edge iff:
 *  - all incident faces are boundary faces
 *  - at least two incident faces exist
 *  - at least two of them have different attributes
 */
template <class Mesh>
bool isBoundaryInterfaceEdge(const Mesh& mesh, Index edgeIdx)
{
  const auto& inc = mesh.getConnectivity().getIncidence({ 1, 2 }, edgeIdx);

  if (inc.size() < 2)
    return false;

  Optional<Attribute> refAttr;

  for (Index f : inc)
  {
    if (!mesh.isBoundary(f))
      return false;

    const Attribute a = mesh.getAttribute(2, f);
    if (!refAttr)
    {
      refAttr = a;
    }
    else if (*refAttr != a)
    {
      return true; // attribute jump on boundary
    }
  }

  return false;
}

int main(int, char**)
{
  Mesh mesh;
  mesh.load("CoronaryArtery_Fluid.mesh", IO::FileFormat::MEDIT);

  // Required incidences
  mesh.getConnectivity().compute(2, 3); // face -> cell (for isBoundary)
  mesh.getConnectivity().compute(1, 2); // edge -> face
  mesh.getConnectivity().compute(1, 1); // edge adjacency (for CCL)

  auto ccl = mesh.ccl(
    1,
    [&](const Polytope& a, const Polytope& b)
    {
      return isBoundaryInterfaceEdge(mesh, a.getIndex()) &&
             isBoundaryInterfaceEdge(mesh, b.getIndex());
    },
    [&](const Polytope& p)
    {
      return isBoundaryInterfaceEdge(mesh, p.getIndex());
    });

  std::cout << "Found " << ccl.getCount()
            << " boundary interface edge components." << std::endl;

  Attribute label = 200;
  for (const auto& cc : ccl)
  {
    for (Index e : cc)
      mesh.setAttribute({ 1, e }, label);
    ++label;
  }

  mesh.save("BoundaryInterfaceEdges.mesh", IO::FileFormat::MEDIT);
  return 0;
}
