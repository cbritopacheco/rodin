/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/* Validate the spherical Kelvin ball on its chamber and sewn domains. */
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include <Rodin/Alert/Raise.h>
#include <Rodin/Alert/Success.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>

#include "Metrics.h"
#include "SewedOutput.h"
#include "Sphere.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace KelvinBall
{
  class SphereValidation
  {
    public:
      SphereValidation(int argc, char** argv)
        : m_argc(argc),
          m_argv(argv)
      {}

      int run();

    private:
      void print(const char* name, const Values& values) const
      {
        std::cout << name << ": k=" << values.k << ", c=" << values.c
                  << ", q=" << values.q << ", rho=" << values.rho;
        if (values.nitscheJump > 0)
          std::cout << ", nitsche_jump=" << values.nitscheJump;
        std::cout << '\n';
      }

      Real relativeDifference(Real chamber, Real full) const
      {
        return std::abs(chamber - full) /
          std::max({std::abs(chamber), std::abs(full), Real(1e-30)});
      }

      void checkReconstruction(Mesh& mesh) const
      {
        mesh.getConnectivity().compute(2, 3);
        size_t openCuts = 0;
        for (auto face = mesh.getBoundary(); face; ++face)
        {
          const auto attribute = face->getAttribute();
          if (!attribute || (*attribute != Gamma && *attribute != Outer))
            ++openCuts;
        }
        if (openCuts != 0)
          throw std::runtime_error(
            "The reconstructed mesh has nonconforming internal chamber cuts.");
      }

      int m_argc;
      char** m_argv;
  };
}

int KelvinBall::SphereValidation::run()
{
  const int argc = m_argc;
  char** argv = m_argv;
  Configuration configuration;
  bool saveMesh = false;
  for (int argument = 1; argument < argc; ++argument)
  {
    const std::string_view option(argv[argument]);
    if (configuration.parse(option))
      continue;
    if (option == "--save-mesh")
      saveMesh = true;
    else if (option == "--help")
    {
      std::cout
        << "Usage: " << argv[0] << " [options]\n"
        << "  --n=<points>              Background points per edge (default: 13).\n"
        << "  --h=<size>                Requested mesh size; alternative to --n.\n"
        << "  --outer-radius=<value>    Chamber outer radius (default: 2).\n"
        << "  --penalty=<value>         Rotational Nitsche penalty (default: 320).\n"
        << "  --stabilization=<value>   P1--P1 stabilization factor (default: 0.05).\n"
        << "  --save-mesh               Write KelvinBallSphere.mesh.\n";
      return 0;
    }
    else
      throw std::runtime_error("Unknown KelvinBallSphere option: " + std::string(option));
  }
  configuration.finalize();

  const Parameters parameters{configuration.getH(), configuration.nitschePenalty,
    configuration.stabilizationFactor};

  auto sphere = Sphere(configuration).discretize(true);
  auto& chamber = sphere.mesh;
  if (saveMesh)
    chamber.save("KelvinBallSphere.mesh", IO::FileFormat::MEDIT);
  prepare(chamber);
  SubMesh chamberFluid = chamber.trim(Obstacle);
  Metrics metrics(parameters);
  const Values chamberValues = metrics.evaluateChamber(chamberFluid);

  SewedOutput reconstructed(chamberFluid, FlatSet<Attribute>{Gamma, Outer});
  checkReconstruction(reconstructed.getMesh());
  const Values fullValues = metrics.evaluateSewed(reconstructed.getMesh());

  print("chamber", chamberValues);
  print("full", fullValues);
  std::cout << "relative differences: k="
            << relativeDifference(chamberValues.k, fullValues.k)
            << ", c=" << relativeDifference(chamberValues.c, fullValues.c)
            << ", q=" << relativeDifference(chamberValues.q, fullValues.q)
            << ", rho=" << relativeDifference(chamberValues.rho, fullValues.rho) << '\n';

  IO::XDMF output("KelvinBallSphere");
  output.setMesh(reconstructed.getMesh());
  output.write(0).flush();
  output.close();

  Alert::Success() << "Wrote KelvinBallSphere.xdmf" << Alert::Raise;
  return 0;
}

int main(int argc, char** argv)
{
  return KelvinBall::SphereValidation(argc, argv).run();
}
