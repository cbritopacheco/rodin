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
  Scalar elevation;
  Scalar period;
};

int main(int, char**)
{
  // constexpr int width = 5e4, height = 5e4; // meters
  constexpr Scalar flatness = 4.0;

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

  auto noise = [&](Scalar x, Scalar y) { return gen.GetNoise(x, y) / 2.0 + 0.5; };

  Scalar avgElementVolume = 0;
  for (size_t i = 0; i < topography.getElementCount(); i++)
    avgElementVolume += topography.getElement(i)->getMeasure();
  avgElementVolume /= topography.getElementCount();

  Scalar maxElevation =
    std::max_element(octaves.begin(), octaves.end(),
        [](auto&& lhs, auto& rhs) { return lhs.elevation < rhs.elevation; })->elevation;

  Scalar totalElevation = 0.0;
  for (const auto& octave : octaves)
   totalElevation += octave.elevation;

  auto elevation =
    [&](const Point& p)
    {
      Scalar nx = p.x(), ny = p.y();
      Scalar e = 0.0;
      for (const auto& octave : octaves)
      {
       Scalar f = avgElementVolume / (octave.period * octave.period);
       e += octave.elevation * noise(f * nx, f * ny);
      }
      e /= totalElevation;
      assert(e < 1.0 + std::numeric_limits<Scalar>::epsilon());
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
