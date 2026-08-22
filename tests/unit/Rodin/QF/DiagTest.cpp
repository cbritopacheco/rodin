#include <gtest/gtest.h>
#include <chrono>
#include "Rodin/QF/NodeElimination.h"
#include "QuadratureInvariants.h"
using namespace Rodin; using namespace Rodin::QF;
using namespace Rodin::Geometry; using namespace Rodin::Tests::QF;
TEST(Diag, ReduceStillWorks)
{
  for (size_t p = 1; p <= 10; ++p)
  {
    const auto seed = NodeElimination::productSeed(Polytope::Type::Triangle, p);
    const auto t0 = std::chrono::steady_clock::now();
    const auto out = NodeElimination::reduce(Polytope::Type::Triangle, p, seed);
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    const auto rep = exactnessSweep(out, Polytope::Type::Triangle, p);
    std::printf("d%-2zu %3zu -> %3zu  err %.1e  ok %d  %.1fs\n", p, seed.getSize(),
      out.getSize(), rep.worstRelativeError,
      (int)(allWeightsPositive(out) && allPointsInside(out, Polytope::Type::Triangle)), s);
    std::fflush(stdout);
  }
}
