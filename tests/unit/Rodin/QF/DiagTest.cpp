#include <gtest/gtest.h>
#include "Rodin/QF/NodeElimination.h"
using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
TEST(Diag, WhereDoesDegree15Fail)
{
  for (size_t p : {14u, 15u, 16u})
  {
    auto seed = NodeElimination::productSeed(Polytope::Type::Triangle, p);
    NodeElimination::Diagnostics d;
    const auto out = NodeElimination::reduce(Polytope::Type::Triangle, p, seed, 0, &d);
    std::printf(
      "d%zu: %zu -> %zu | tried %zu, notConv %zu, notAdm %zu (negw %zu, outside %zu)\n",
      p, seed.getSize(), out.getSize(), d.candidatesTried, d.notConverged,
      d.notAdmissible, d.worstNegativeWeights, d.worstExteriorPoints);
    std::fflush(stdout);
  }
}
