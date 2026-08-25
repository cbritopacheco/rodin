/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/**
 * @file
 * @brief Development driver writing the WitherdenVincent coefficient tables.
 *
 * Not a test. Prints one array of WitherdenVincent.cpp, as pasteable C++, for
 * the named element.
 *
 * @par Usage
 * @code
 *   RodinWitherdenVincent <element> [maxDegree] [seconds] [from] [jobs]
 * @endcode
 * where @c element is one of tri, quad, tet, wedge, pyr, hex. Only the element
 * is required, and it names the array written:
 * @code
 *   RodinWitherdenVincent wedge > wedge.txt      # s_wedge
 *   RodinWitherdenVincent tet 8 > tetrahedron.txt
 * @endcode
 * These are the fully symmetric rules, built from orbits of the element's
 * symmetry group, and they cover every element. The asymmetric rules of
 * XiaoGimbutas.cpp are a different construction reaching different counts, and
 * are written by RodinXiaoGimbutas.
 *
 * @par Parameters that matter
 * - @c seconds is a per-strength deadline and the only argument that changes
 *   what comes out: it truncates the search, so a strength that would have
 *   been solved comes back larger or not at all. Zero, the default, means no
 *   deadline. Set it only for an unattended sweep meant to cover ground rather
 *   than finish a particular strength.
 * - @c jobs defaults to the machine. It changes how long a search takes and
 *   not what it returns, since the ordering of decompositions and the seeding
 *   are fixed --- but it matters through the deadline. The hexahedron at
 *   strengths six to nine was once recorded as unreachable, having been given
 *   two threads and a fifteen minute deadline; it is found in about two
 *   minutes with the machine to itself.
 *
 * @par What to expect
 * Two cost regimes, and confusing them misdirects the effort. The
 * two-dimensional elements are slow per decomposition and few to enumerate:
 * the triangle at strength twenty spends about ten seconds on each of 756
 * candidates. The three-dimensional elements are the reverse --- the prism at
 * strength six enumerates 30 806 candidates at about a tenth of a second each,
 * and its 28-point rule lies well down that list.
 */
#include <string>

#include "PublishedCounts.h"
#include "TableEmitter.h"
#include "Rodin/QF/SymmetricRuleGenerator.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Tests::QF;

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::fprintf(stderr,
      "usage: RodinWitherdenVincent <tri|quad|tet|wedge|pyr|hex>"
      " [maxDegree] [seconds] [from] [jobs]\n");
    return 2;
  }
  const auto g = elementOf(argv[1]);
  const auto& published = witherdenVincentCounts(g);
  const size_t maxDegree = (argc > 2)
    ? std::stoul(argv[2])
    : (published.empty() ? 10 : published.rbegin()->first);
  const double seconds = (argc > 3) ? std::stod(argv[3]) : 0.0;
  const size_t from = (argc > 4) ? std::stoul(argv[4]) : 1;
  const size_t jobs = (argc > 5) ? std::stoul(argv[5]) : 0;

  for (size_t degree = from; degree <= maxDegree; ++degree)
  {
    const auto rule =
      SymmetricRuleGenerator::search(g, degree, 0, 64, 1e-12, jobs, seconds);
    if (!rule.converged || !rule.admissible)
    {
      std::printf("    // degree %zu: NOT FOUND (residual %.2e)%s\n", degree,
        rule.residual, (seconds > 0) ? " -- deadline reached; retry without one" : "");
      std::fflush(stdout);
      continue;
    }
    emitEntry(rule, g, degree);
  }
  return 0;
}
