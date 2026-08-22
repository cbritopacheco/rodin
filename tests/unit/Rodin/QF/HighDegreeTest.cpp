#include <gtest/gtest.h>
#include <chrono>
#include "Rodin/QF/NodeElimination.h"
#include "QuadratureInvariants.h"
using namespace Rodin; using namespace Rodin::QF;
using namespace Rodin::Geometry; using namespace Rodin::Tests::QF;
TEST(Sweep, TriangleToTwelve)
{
  std::printf("%4s %5s %6s %6s %10s %4s %7s\n","deg","seed","final","bound","err","ok","t");
  for (size_t p = 1; p <= 12; ++p)
  {
    const auto seed = NodeElimination::productSeed(Polytope::Type::Triangle, p);
    const auto t0 = std::chrono::steady_clock::now();
    const auto out = NodeElimination::reduce(Polytope::Type::Triangle, p, seed);
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    const auto rep = exactnessSweep(out, Polytope::Type::Triangle, p);
    const size_t bound = ((p+1)*(p+2)/2 + 2) / 3;
    std::printf("%4zu %5zu %6zu %6zu %10.1e %4d %6.1fs\n", p, seed.getSize(),
      out.getSize(), bound, rep.worstRelativeError,
      (int)(allWeightsPositive(out) && allPointsInside(out, Polytope::Type::Triangle)), s);
    std::fflush(stdout);
  }
}
