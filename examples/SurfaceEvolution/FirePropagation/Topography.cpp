#include <Rodin/Math.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

#include "FastNoiseLite.h"

using namespace Rodin;
using namespace Geometry;
using namespace Variational;
using namespace External;

const char* meshfile = "../resources/examples/SurfaceEvolution/FirePropagation/GroundPlane.mfem.mesh";

struct Octave
{
  Real elevation;
  Real period;
};

int main(int, char**)
{
  // constexpr int width = 5e4, height = 5e4; // meters
  constexpr Real flatness = 4.0;

  std::vector<Octave> octaves = {
   {500.0, 500},
   {1000.0, 1000},
   {2000.0, 1000},
   {3000.0, 2000},
   {4000.0, 3000},
  };

  Mesh topography;
  topography.load(meshfile);

  srand(std::time(0));
  FastNoiseLite gen(rand());
  gen.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

  auto noise = [&](Real x, Real y) { return gen.GetNoise(x, y) / 2.0 + 0.5; };

  Real avgElementVolume = 0;
  for (size_t i = 0; i < topography.getCellCount(); i++)
    avgElementVolume += topography.getCell(i)->getMeasure();
  avgElementVolume /= topography.getCellCount();

  Real maxElevation =
    std::max_element(octaves.begin(), octaves.end(),
        [](auto&& lhs, auto& rhs) { return lhs.elevation < rhs.elevation; })->elevation;

  Real totalElevation = 0.0;
  for (const auto& octave : octaves)
   totalElevation += octave.elevation;

  auto elevation =
    [&](const Point& p)
    {
      Real nx = p.x(), ny = p.y();
      Real e = 0.0;
      for (const auto& octave : octaves)
      {
       Real f = avgElementVolume / (octave.period * octave.period);
       e += octave.elevation * noise(f * nx, f * ny);
      }
      e /= totalElevation;
      assert(e < 1.0 + std::numeric_limits<Real>::epsilon());
      return maxElevation * std::pow(e, flatness);
    };

  P1 vh(topography, topography.getSpaceDimension());
  GridFunction displacement(vh);
  displacement = VectorFunction{0, 0, elevation};
  topography.displace(displacement);

  P1 sh(topography);
  GridFunction h(sh);
  h = [](const Point& p) { return p.z(); };

  topography.save("Topography.mesh");
  h.save("Topography.gf");

  return 0;
}
