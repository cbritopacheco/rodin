#include <gtest/gtest.h>
#include <chrono>
#include "Rodin/QF/NodeElimination.h"
using namespace Rodin; using namespace Rodin::QF; using namespace Rodin::Geometry;
TEST(Diag, DoesRefineConvergeAtHighDegree)
{
  for (size_t p : {12u, 15u, 18u, 20u})
  {
    auto seed = NodeElimination::productSeed(Polytope::Type::Triangle, p);
    auto trial = seed;
    trial.points.erase(trial.points.begin());
    trial.weights.erase(trial.weights.begin());
    const auto t0 = std::chrono::steady_clock::now();
    const bool conv = NodeElimination::refine(Polytope::Type::Triangle, p, trial);
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    const bool adm = NodeElimination::isAdmissible(Polytope::Type::Triangle, trial);
    std::printf("deg %2zu: n=%3zu conv=%d adm=%d (%.1fs)\n",
      p, trial.getSize(), (int)conv, (int)adm, s);
    std::fflush(stdout);
  }
}
