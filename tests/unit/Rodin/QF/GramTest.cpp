#include <gtest/gtest.h>
#include "Rodin/QF/SymmetricRuleSolver.h"
using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
TEST(G, GramCondition)
{
  for (size_t p : {6u, 10u, 12u, 14u, 16u, 20u})
  {
    SymmetricRuleSolver::MomentData md(Polytope::Type::Triangle, 2, p);
    std::printf("deg %2zu: monomials %3zu  cholesky %s\n", p, md.alphas.size(),
      md.chol.size() > 0 ? "OK" : "FAILED");
  }
}
