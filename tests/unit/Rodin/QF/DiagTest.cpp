#include <gtest/gtest.h>
#include <chrono>
#include "Rodin/QF/NodeElimination.h"
#include "QuadratureInvariants.h"
using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;
TEST(Diag, TetD5Hard)
{
  // At the counting bound the system is square, so the basin of the isolated
  // solution is small. Many manifold moves give many starting points.
  for (size_t div : {size_t(50), size_t(400)})
  {
    const auto seed = NodeElimination::productSeed(Polytope::Type::Tetrahedron, 5);
    const auto t0 = std::chrono::steady_clock::now();
    const auto out = NodeElimination::reduce(
      Polytope::Type::Tetrahedron, 5, seed, 0, nullptr, 200, 8000, true, div);
    const double s =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const auto rep = exactnessSweep(out, Polytope::Type::Tetrahedron, 5);
    std::printf("diversifications %4zu -> %2zu points (bound 14) err %.1e ok %d %.0fs\n",
      div, out.getSize(), rep.worstRelativeError,
      (int)(allWeightsPositive(out) && allPointsInside(out, Polytope::Type::Tetrahedron)),
      s);
    std::fflush(stdout);
  }
}
