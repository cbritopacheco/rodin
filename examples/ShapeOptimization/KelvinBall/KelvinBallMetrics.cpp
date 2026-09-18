/* Validate Kelvin-ball resistance metrics on a chamber and its reconstruction. */
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include <Rodin/Alert/Raise.h>
#include <Rodin/Alert/Success.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>

#include "Metrics.h"
#include "SewedOutput.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  void print(const char* name, const KelvinBall::Values& values)
  {
    std::cout << name << ": k=" << values.k << ", c=" << values.c << ", q=" << values.q
              << ", rho=" << values.rho;
    if (values.nitscheJump > 0)
      std::cout << ", nitsche_jump=" << values.nitscheJump;
    std::cout << '\n';
  }

  Real relativeDifference(Real chamber, Real full)
  {
    return std::abs(chamber - full) /
      std::max({std::abs(chamber), std::abs(full), Real(1e-30)});
  }

  void checkReconstruction(KelvinBall::Mesh& mesh)
  {
    mesh.getConnectivity().compute(2, 3);
    size_t openCuts = 0;
    for (auto face = mesh.getBoundary(); face; ++face)
    {
      const auto attribute = face->getAttribute();
      if (!attribute ||
        (*attribute != KelvinBall::Gamma && *attribute != KelvinBall::Outer))
        ++openCuts;
    }
    if (openCuts != 0)
      throw std::runtime_error(
        "The reconstructed mesh has nonconforming internal chamber cuts.");
  }
}

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::cerr << "usage: KelvinBallMetrics chamber.mesh h [penalty] [stabilization]\n";
    return 1;
  }

  KelvinBall::Parameters parameters;
  parameters.h = std::stod(argv[2]);
  if (argc > 3)
    parameters.nitschePenalty = std::stod(argv[3]);
  if (argc > 4)
    parameters.stabilizationFactor = std::stod(argv[4]);
  if (!(parameters.h > 0 && parameters.nitschePenalty > 0 &&
        parameters.stabilizationFactor >= 0))
    throw std::runtime_error("Invalid Kelvin-ball metric parameters.");

  KelvinBall::Mesh chamber;
  chamber.load(argv[1], IO::FileFormat::MEDIT);
  KelvinBall::prepare(chamber);
  SubMesh chamberFluid = chamber.trim(KelvinBall::Obstacle);
  const KelvinBall::Values chamberValues =
    KelvinBall::evaluateChamber(chamberFluid, parameters);

  auto reconstructed = KelvinBall::sew(
    chamberFluid, FlatSet<Attribute>{KelvinBall::Gamma, KelvinBall::Outer});
  checkReconstruction(reconstructed.mesh);
  const KelvinBall::Values fullValues =
    KelvinBall::evaluateFull(reconstructed.mesh, parameters);

  print("chamber", chamberValues);
  print("full", fullValues);
  std::cout << "relative differences: k="
            << relativeDifference(chamberValues.k, fullValues.k)
            << ", c=" << relativeDifference(chamberValues.c, fullValues.c)
            << ", q=" << relativeDifference(chamberValues.q, fullValues.q)
            << ", rho=" << relativeDifference(chamberValues.rho, fullValues.rho) << '\n';

  IO::XDMF output("KelvinBallMetrics");
  output.setMesh(reconstructed.mesh);
  output.write(0).flush();
  output.close();

  Alert::Success() << "Wrote KelvinBallMetrics.xdmf" << Alert::Raise;
  return 0;
}
