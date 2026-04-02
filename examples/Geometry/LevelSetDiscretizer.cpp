/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Types.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Geometry/LevelSetDiscretizerTetrahedra.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  using MeshType = Mesh<Context::Local>;

  // --------------------------------------------------------------------------
  // 1. Background mesh: uniform tetrahedral grid on the unit cube
  // --------------------------------------------------------------------------
  MeshType mesh;
  mesh = mesh.UniformGrid(
    Polytope::Type::Tetrahedron,
    { 24, 24, 24 } // refine as needed
  );
  mesh.scale(1.0 / (24 - 1)); // unit cube

  mesh.getConnectivity().compute(3, 2);
  mesh.getConnectivity().compute(2, 3);
  mesh.getConnectivity().compute(3, 1);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);

  // Give all cells a base attribute so split() has something to remap.
  static constexpr Attribute bulk = 1;
  for (auto c = mesh.getCell(); !c.end(); ++c)
    mesh.setAttribute({3, c->getIndex()}, bulk);

  // --------------------------------------------------------------------------
  // 2. P1 level set for a sphere
  //    Use phi = |x-c|^2 - R^2 to avoid a sqrt.
  // --------------------------------------------------------------------------
  P1 phiSpace(mesh);
  GridFunction phi(phiSpace);

  const Real cx = 0.5;
  const Real cy = 0.5;
  const Real cz = 0.5;
  const Real R  = 0.28;

  for (Index i = 0; i < mesh.getVertexCount(); ++i)
  {
    const auto& x = mesh.getVertexCoordinates(i);
    const Real dx = x(0) - cx;
    const Real dy = x(1) - cy;
    const Real dz = x(2) - cz;
    phi[i] = dx * dx + dy * dy + dz * dz - R * R;
  }

  // Optional: save the background mesh + level set
  {
    IO::XDMF xdmf("SphereLevelSetBackground");
    xdmf.grid().setMesh(mesh).add("phi", phi, IO::XDMF::Center::Node);
    xdmf.write();
    xdmf.close();
  }

  // --------------------------------------------------------------------------
  // 3. Conforming split
  // --------------------------------------------------------------------------
  static constexpr Attribute inside  = 10;
  static constexpr Attribute outside = 20;
  static constexpr Attribute gamma   = 30;

  LevelSetDiscretizerTetrahedra mt(phi);
  mt
    .setSignTolerance(1e-14)
    .setSnapTolerance(1e-14)
    .setInterface(2, gamma)   // mark interface faces
    .setInterface(1, gamma)   // mark interface edges
    .split(3, bulk, {inside, outside})
    ;

  MeshType out = mt.discretize();
  out.save("split.mesh", IO::FileFormat::MEDIT);

  // Build useful connectivities on the output
  out.getConnectivity().compute(2, 3);
  out.getConnectivity().compute(1, 3);
  out.getConnectivity().compute(1, 2);
  out.getConnectivity().compute(0, 0);

  // --------------------------------------------------------------------------
  // 4. Save result
  // --------------------------------------------------------------------------
  {
    IO::XDMF xdmf("SphereLevelSetSplit");
    xdmf.grid().setMesh(out);
    xdmf.write();
    xdmf.close();
  }

  return 0;
}
