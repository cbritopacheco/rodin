/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include <Rodin/Types.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MMG.h>
#include <Rodin/Geometry/LevelSetDiscretizerTetrahedra.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  using MeshType = MMG::Mesh;

  Alert::Info() << "Initializing evolving sphere example." << Alert::Raise;

  // --------------------------------------------------------------------------
  // Parameters
  // --------------------------------------------------------------------------
  static constexpr Attribute bulk    = 1;
  static constexpr Attribute inside  = 10;
  static constexpr Attribute outside = 20;
  static constexpr Attribute gamma   = 30;

  const Real cx = 0.5;
  const Real cy = 0.5;
  const Real cz = 0.5;

  const Real R0    = 0.18;
  const Real speed = 0.015;  // radius change per unit time
  const Real dt    = 0.05;
  const size_t nt  = 20;

  const Real hmax  = 0.06;
  const Real hmin  = 0.015;
  const Real hausd = 0.005;

  // Background mesh resolution
  const size_t nx = 24;
  const size_t ny = 24;
  const size_t nz = 24;

  Alert::Info() << "Parameters set." << Alert::Raise;

  // --------------------------------------------------------------------------
  // Initial background mesh
  // --------------------------------------------------------------------------
  Alert::Info() << "Building initial background mesh." << Alert::Raise;
  MeshType th;
  th = th.UniformGrid(Polytope::Type::Tetrahedron, { nx, ny, nz });
  th.scale(1.0 / (nx - 1)); // unit cube

  Alert::Info() << "Assigning bulk attribute to initial cells." << Alert::Raise;
  for (auto c = th.getCell(); !c.end(); ++c)
    th.setAttribute({3, c->getIndex()}, bulk);

  Alert::Info() << "Saving initial mesh." << Alert::Raise;
  th.save("out/Initial.mesh", IO::FileFormat::MEDIT);

  // --------------------------------------------------------------------------
  // XDMF output
  // --------------------------------------------------------------------------
  Alert::Info() << "Initializing XDMF output." << Alert::Raise;
  IO::XDMF xdmf("out/EvolvingSphere");

  auto domain = xdmf.grid("domain");
  domain.setMesh(th, IO::XDMF::MeshPolicy::Transient);

  // --------------------------------------------------------------------------
  // Time loop
  // --------------------------------------------------------------------------
  for (size_t it = 0; it < nt; ++it)
  {
    const Real t = it * dt;
    const Real R = R0 + speed * t;

    Alert::Info() << "----- Iteration: " << it
                  << " | time = " << t
                  << " | radius = " << R
                  << Alert::Raise;

    // ------------------------------------------------------------------------
    // Build/connectivity on current mesh
    // ------------------------------------------------------------------------
    Alert::Info() << "   | Computing connectivity on current mesh." << Alert::Raise;
    auto& conn = th.getConnectivity();
    conn.compute(3, 2);
    conn.compute(2, 3);
    conn.compute(3, 1);
    conn.compute(2, 1);
    conn.compute(1, 0);

    // ------------------------------------------------------------------------
    // Level set on current mesh:
    // phi(x) = |x-c|^2 - R^2
    // negative inside sphere, positive outside
    // ------------------------------------------------------------------------
    Alert::Info() << "   | Building P1 level-set field on current mesh." << Alert::Raise;
    P1 phiSpace(th);
    GridFunction phi(phiSpace);

    Alert::Info() << "   | Evaluating level-set values." << Alert::Raise;
    for (Index i = 0; i < th.getVertexCount(); ++i)
    {
      const auto& x = th.getVertexCoordinates(i);
      const Real dx = x(0) - cx;
      const Real dy = x(1) - cy;
      const Real dz = x(2) - cz;
      phi[i] = dx * dx + dy * dy + dz * dz - R * R;
    }

    Alert::Info() << "   | Saving current level-set field." << Alert::Raise;
    {
      std::ostringstream oss;
      oss << "SpherePhi_" << std::setw(4) << std::setfill('0') << it;
      IO::XDMF phiOut("out/" + oss.str());
      phiOut.grid().setMesh(th).add("phi", phi, IO::XDMF::Center::Node);
      phiOut.write();
      phiOut.close();
    }

    // ------------------------------------------------------------------------
    // Conforming level-set discretization
    // ------------------------------------------------------------------------
    Alert::Info() << "   | Configuring level-set discretizer." << Alert::Raise;
    LevelSetDiscretizerTetrahedra lsd(phi);
    lsd
      .setSignTolerance(1e-14)
      .setSnapTolerance(1e-14)
      .setInterface(2, gamma)
      .split(3, bulk,    {inside, outside})
      .split(3, inside,  {inside, outside})
      .split(3, outside, {inside, outside});

    Alert::Info() << "   | Discretizing level-set conformingly." << Alert::Raise;
    MeshType split = lsd.discretize();

    Alert::Info() << "   | Saving split mesh." << Alert::Raise;
    split.save(
      "out/Split_" + std::to_string(it) + ".mesh",
      IO::FileFormat::MEDIT
    );

    // ------------------------------------------------------------------------
    // Recompute connectivity on the split mesh
    // ------------------------------------------------------------------------
    Alert::Info() << "   | Computing connectivity on split mesh." << Alert::Raise;
    auto& sconn = split.getConnectivity();
    sconn.compute(3, 2);
    sconn.compute(2, 3);
    sconn.compute(3, 1);
    sconn.compute(2, 1);
    sconn.compute(1, 0);

    // ------------------------------------------------------------------------
    // MMG optimization / parametrization of the fitted mesh
    // ------------------------------------------------------------------------
    Alert::Info() << "   | Optimizing split mesh with MMG." << Alert::Raise;
    try
    {
      MMG::Optimizer()
        .setHMax(hmax)
        .setHMin(hmin)
        .setHausdorff(hausd)
        .setAngleDetection(false)
        .optimize(split);
    }
    catch (const Alert::Exception& e)
    {
      Alert::Warning()
        << "MMG optimization failed at iteration " << it
        << ". Keeping the split mesh without optimization."
        << Alert::Raise;
    }

    Alert::Info() << "   | Saving optimized mesh." << Alert::Raise;
    split.save(
      "out/Optimized." + std::to_string(it) + ".mesh",
      IO::FileFormat::MEDIT
    );

    // ------------------------------------------------------------------------
    // Transient XDMF output
    // ------------------------------------------------------------------------
    Alert::Info() << "   | Recomputing level-set field on optimized mesh." << Alert::Raise;
    P1 phiSplitSpace(split);
    GridFunction phiSplit(phiSplitSpace);

    for (Index i = 0; i < split.getVertexCount(); ++i)
    {
      const auto& x = split.getVertexCoordinates(i);
      const Real dx = x(0) - cx;
      const Real dy = x(1) - cy;
      const Real dz = x(2) - cz;
      phiSplit[i] = dx * dx + dy * dy + dz * dz - R * R;
    }

    Alert::Info() << "   | Writing transient XDMF snapshot." << Alert::Raise;
    domain.clear();
    domain.setMesh(split, IO::XDMF::MeshPolicy::Transient);
    domain.add("phi", phiSplit, IO::XDMF::Center::Node);
    xdmf.write(it).flush();

    Alert::Info() << "   | Updating mesh for next iteration." << Alert::Raise;
    th = std::move(split);
  }

  Alert::Info() << "Closing XDMF output." << Alert::Raise;
  xdmf.close();

  Alert::Success() << "Finished evolving sphere example." << Alert::Raise;
  return 0;
}
