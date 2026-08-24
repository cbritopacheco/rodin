#include <gtest/gtest.h>
#include "Rodin/QF/SymmetricRuleSolver.h"
using namespace Rodin; using namespace Rodin::QF; using namespace Rodin::Geometry;
TEST(Diag, RestartSensitivity)
{
  for (size_t restarts : {size_t(256), size_t(2000)})
  {
    const auto r = SymmetricRuleSolver::solve(Polytope::Type::Triangle, 3,
      {{2, 1}, {2, 1}}, restarts);
    std::printf("restarts %4zu: conv=%d adm=%d res=%.3e used=%zu\n",
      restarts, (int)r.converged, (int)r.admissible, r.residual, r.restarts);
    std::fflush(stdout);
  }
  // and via the search, which passed
  SymmetricRuleSolver::Configuration cfg;
  const auto s = SymmetricRuleSolver::search(Polytope::Type::Triangle, 3, 24, 96, &cfg);
  std::printf("search: conv=%d adm=%d n=%zu\n", (int)s.converged, (int)s.admissible,
    SymmetricRuleSolver::configurationSize(cfg));
}
