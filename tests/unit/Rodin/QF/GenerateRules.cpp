// Development driver: prints the rules the solver finds, for inlining.
#include <cstdio>
#include "Rodin/QF/SymmetricRuleSolver.h"
using namespace Rodin;
using namespace Rodin::QF;
int main(int argc, char** argv)
{
  const auto g = (argc > 1 && std::string(argv[1]) == "tet")
    ? Geometry::Polytope::Type::Tetrahedron : Geometry::Polytope::Type::Triangle;
  const size_t maxDeg = (argc > 2) ? std::stoul(argv[2]) : 8;
  for (size_t p = 1; p <= maxDeg; ++p)
  {
    SymmetricRuleSolver::Configuration cfg;
    const auto r = SymmetricRuleSolver::search(g, p, 40, 96, &cfg);
    if (!(r.converged && r.admissible)) { std::printf("d%zu: none (res %.2e)\n", p, r.residual); continue; }
    std::printf("d%zu n=%zu :", p, SymmetricRuleSolver::configurationSize(cfg));
    for (const auto& o : r.orbits)
    {
      auto b = o.getBarycentric();
      std::sort(b.begin(), b.end());
      std::printf(" {");
      for (size_t i = 0; i + 1 < b.size(); ++i) std::printf("%.17g,", b[i]);
      std::printf("%.17g|w=%.17g}", b.back(), o.getWeight());
    }
    std::printf("\n");
  }
  return 0;
}
