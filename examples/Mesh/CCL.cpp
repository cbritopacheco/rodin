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
  const char* filename = "../resources/mfem/ccl-2d-example.mesh";
  Mesh mesh;
  mesh.load(filename);

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
    mesh.edit(
        [&](ElementView element)
        {
          element.setAttribute(i + 1);
        }, cc);
  }

  Alert::Info() << "Saved mesh to ccl.mesh" << Alert::Raise;

  mesh.save("ccl.mesh");

  return 0;
}
