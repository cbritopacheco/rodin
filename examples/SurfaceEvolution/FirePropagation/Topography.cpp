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

struct Octave
{
  double elevation;
  double period;
};

int main(int, char**)
{
  constexpr int width = 5e4, height = 5e4; // meters
  constexpr double flatness = 4.0;

  std::vector<Octave> octaves = {
   {500.0, 500},
   {1000.0, 1000},
   {2000.0, 1000},
   {3000.0, 2000},
   {4000.0, 3000},
  };

  MMG::Mesh topography;
  topography.initialize(2, 3);
  topography.vertex({0, 0, 0});
  topography.vertex({width, 0, 0});
  topography.vertex({0, height, 0});
  topography.vertex({width, height, 0});
  topography.element(
    Geometry::Type::Triangle, {0, 1, 2});
  topography.element(
    Geometry::Type::Triangle, {1, 2, 3});
  topography.finalize();

  // for (int i = 0; i < 8; i++)
  //  topography.refine();

  srand(std::time(0));
  FastNoiseLite gen(rand());
  gen.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

  auto noise = [&](double x, double y) { return gen.GetNoise(x, y) / 2.0 + 0.5; };

  double elVol = 0;
  // for (int i = 0; i < topography.count<Element>(); i++)
  //  elVol += topography.get<Element>(i).getVolume();
  // elVol /= topography.count<Element>();

  double maxElevation = std::max_element(
    octaves.begin(), octaves.end(),
    [](auto&& lhs, auto& rhs)
    { return lhs.elevation < rhs.elevation; })->elevation;

  double totalElevation = 0.0;
  for (const auto& octave : octaves)
   totalElevation += octave.elevation;

  H1 vh(topography, topography.getSpaceDimension());
  GridFunction displacement(vh);
  displacement =
   VectorFunction{0, 0,
    [&](const Point& p)
    {
      double nx = p.x(), ny = p.y();
      double e = 0.0;
      for (const auto& octave : octaves)
      {
       double f = elVol / (octave.period * octave.period);
       e += octave.elevation * noise(f * nx, f * ny);
      }
      e /= totalElevation;
      assert(e < 1.0 + std::numeric_limits<double>::epsilon());
      return maxElevation * std::pow(e, flatness);
    }};
  topography.displace(displacement);

  H1 sh(topography);
  GridFunction h(sh);
  h = [](const Point& p) { return p.z(); };

  topography.save("topography.mesh");
  h.save("topography.gf");

  return 0;
}
