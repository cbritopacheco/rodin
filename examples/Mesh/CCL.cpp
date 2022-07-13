/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Mesh.h>
#include <Rodin/Alert.h>

using namespace Rodin;

int main(int, char**)
{
  Mesh mesh;
  mesh.load("Omega.mesh");

  Alert::Info() << "Performing CCL on mesh attributes..." << Alert::Raise;

  auto ccs = mesh.ccl(
      [](const Element& el1, const Element& el2)
      {
        return el1.getAttribute() == el2.getAttribute();
      });

  Alert::Info() << ccs.size() << " components found." << Alert::Raise;

  for (size_t i = 0; i < ccs.size(); i++)
  {
    const auto& cc = ccs[i];
    mesh.edit(cc,
        [&](ElementView element)
        {
          element.setAttribute(i + 1);
        });
  }

  Alert::Info() << "Saved mesh to CCL.mesh" << Alert::Raise;

  mesh.save("CCL.mesh");

  return 0;
}
