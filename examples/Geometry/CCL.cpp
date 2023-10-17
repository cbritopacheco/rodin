/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>
#include <Rodin/Alert.h>

using namespace Rodin;
using namespace Geometry;

int main(int, char**)
{
  const char* filename = "../resources/examples/Geometry/CCL.mfem.mesh";
  Mesh mesh;
  mesh.load(filename);
  mesh.getConnectivity().compute(2, 2);

  Alert::Info() << "Performing CCL on mesh attributes..." << Alert::Raise;

  auto ccl = mesh.ccl(
      [](const Polytope& el1, const Polytope& el2)
      {
        return el1.getAttribute() == el2.getAttribute();
      });

  Alert::Info() << ccl.getComponents().size() << " components found." << Alert::Raise;

  size_t cci = 1;
  for (const auto& cc : ccl.getComponents())
  {
    for (const Index i : cc)
      mesh.setAttribute({ mesh.getDimension(), i }, cci);
    cci++;
  }

  Alert::Info() << "Saved mesh to CCL.mesh" << Alert::Raise;

  mesh.save("CCL.mesh");

  return 0;
}
