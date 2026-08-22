#include <gtest/gtest.h>
#include <chrono>
#include "Rodin/QF/NodeElimination.h"
#include "QuadratureInvariants.h"
using namespace Rodin; using namespace Rodin::QF;
using namespace Rodin::Geometry; using namespace Rodin::Tests::QF;
TEST(Diag, TetLow)
{
  for (size_t p = 3; p <= 6; ++p)
  {
    const auto seed = NodeElimination::productSeed(Polytope::Type::Tetrahedron, p);
    const auto t0 = std::chrono::steady_clock::now();
    const auto out = NodeElimination::reduce(Polytope::Type::Tetrahedron, p, seed);
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    const auto rep = exactnessSweep(out, Polytope::Type::Tetrahedron, p);
    std::printf("tet d%zu %3zu -> %3zu  err %.1e ok %d  %.0fs\n", p, seed.getSize(),
      out.getSize(), rep.worstRelativeError,
      (int)(allWeightsPositive(out) && allPointsInside(out, Polytope::Type::Tetrahedron)), s);
    std::fflush(stdout);
  }
}
