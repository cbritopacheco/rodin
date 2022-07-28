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

  Alert::Info() << "Trimming mesh..." << Alert::Raise;

  auto trimmed = mesh.trim(2);

  Alert::Info() << "Saved mesh to trimmed.mesh" << Alert::Raise;

  trimmed.save("trimmed.mesh");
}
