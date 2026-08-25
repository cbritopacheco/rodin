/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/**
 * @file
 * @brief Development driver emitting inlinable quadrature coefficients.
 *
 * Not a test. Prints one coefficient table, as pasteable C++, for a named
 * family and element. Which table is produced is exactly the pair of
 * arguments: the family selects the file, the element selects the array
 * inside it.
 *
 * @par Usage
 * @code
 *   RodinGenerateTables <family> <element> [maxDegree] [seconds] [from] [jobs]
 * @endcode
 *
 * | family | file it feeds        | elements                            |
 * |--------|----------------------|-------------------------------------|
 * | @c wv  | WitherdenVincent.cpp | tri quad tet wedge pyr hex          |
 * | @c xg  | XiaoGimbutas.cpp     | tri tet                             |
 *
 * The families are different constructions and are not interchangeable. @c wv
 * builds fully symmetric rules from orbits of the element's symmetry group,
 * and covers every element. @c xg removes nodes from a larger rule, which
 * destroys symmetry but reaches lower counts, and applies only to simplices
 * --- which is also where Xiao and Gimbutas publish. Their published counts
 * differ for the same element and strength, so a table is only ever compared
 * against its own family.
 *
 * @par Examples
 * @code
 *   RodinGenerateTables wv wedge          # s_wedge of WitherdenVincent.cpp
 *   RodinGenerateTables xg tri            # s_triangle of XiaoGimbutas.cpp
 *   RodinGenerateTables wv tet 8          # strengths one to eight only
 * @endcode
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
 * Every rule is checked against the moment oracle before being printed, so a
 * rule that is not exact, positive and interior never reaches a table.
 */
#include <cstdio>
#include <string>
#include <vector>

#include "PublishedCounts.h"
#include "QuadratureInvariants.h"
#include "Rodin/QF/NodeElimination.h"
#include "Rodin/QF/SymmetricRuleGenerator.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Tests::QF;

namespace
{
  Geometry::Polytope::Type elementOf(const std::string& name)
  {
    if (name == "tet")
      return Geometry::Polytope::Type::Tetrahedron;
    if (name == "wedge")
      return Geometry::Polytope::Type::Wedge;
    if (name == "pyr")
      return Geometry::Polytope::Type::Pyramid;
    if (name == "quad")
      return Geometry::Polytope::Type::Quadrilateral;
    if (name == "hex")
      return Geometry::Polytope::Type::Hexahedron;
    return Geometry::Polytope::Type::Triangle;
  }

  /// @brief How far each family's table is worth taking by default: as far as
  /// that family publishes.
  size_t publishedThrough(const std::string& family, Geometry::Polytope::Type g)
  {
    const auto& table =
      (family == "xg") ? xiaoGimbutasCounts(g) : witherdenVincentCounts(g);
    return table.empty() ? 10 : table.rbegin()->first;
  }

  /// @brief Prints one degree as a table entry.
  template <class Rule>
  void emit(const Rule& rule, size_t degree, size_t d, Real error)
  {
    std::printf(
      "    // degree %zu: %zu points, oracle %.2e\n", degree, rule.getSize(), error);
    std::printf("      {\n");
    for (size_t i = 0; i < rule.getSize(); ++i)
    {
      std::printf("        ");
      for (size_t k = 0; k < d; ++k)
        std::printf("%.17g, ", rule.getPoint(i)[static_cast<Eigen::Index>(k)]);
      std::printf("%.17g,\n", rule.getWeight(i));
    }
    std::printf("      },\n");
    std::fflush(stdout);
  }
}

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    std::fprintf(stderr,
      "usage: RodinGenerateTables <wv|xg> <tri|quad|tet|wedge|pyr|hex>"
      " [maxDegree] [seconds] [from] [jobs]\n"
      "  wv  symmetric rules, every element      -> WitherdenVincent.cpp\n"
      "  xg  node elimination, simplices only    -> XiaoGimbutas.cpp\n");
    return 2;
  }

  const std::string family = argv[1];
  const std::string which = argv[2];
  if (family != "wv" && family != "xg")
  {
    std::fprintf(stderr, "unknown family '%s': expected wv or xg\n", family.c_str());
    return 2;
  }
  const auto g = elementOf(which);
  if (family == "xg" && g != Geometry::Polytope::Type::Triangle &&
    g != Geometry::Polytope::Type::Tetrahedron)
  {
    std::fprintf(stderr,
      "the xg family covers only the simplices, which is where it is"
      " published; use wv for '%s'\n",
      which.c_str());
    return 2;
  }

  const size_t maxDegree = (argc > 3) ? std::stoul(argv[3]) : publishedThrough(family, g);
  const double seconds = (argc > 4) ? std::stod(argv[4]) : 0.0;
  const size_t from = (argc > 5) ? std::stoul(argv[5]) : 1;
  const size_t jobs = (argc > 6) ? std::stoul(argv[6]) : 0;
  const size_t d = Geometry::Polytope::Traits(g).getDimension();

  for (size_t degree = from; degree <= maxDegree; ++degree)
  {
    if (family == "xg")
    {
      const auto seed = NodeElimination::productSeed(g, degree);
      const auto rule = NodeElimination::reduce(g, degree, seed);
      const auto report = exactnessSweep(rule, g, degree);
      if (report.worstRelativeError >= 1e-12 || !allWeightsPositive(rule) ||
        !allPointsInside(rule, g, 1e-12))
      {
        std::printf("    // degree %zu: REJECTED (oracle %.2e)\n", degree,
          report.worstRelativeError);
        std::fflush(stdout);
        continue;
      }
      emit(rule, degree, d, report.worstRelativeError);
      continue;
    }

    const auto rule =
      SymmetricRuleGenerator::search(g, degree, 0, 64, 1e-12, jobs, seconds);
    if (!rule.converged || !rule.admissible)
    {
      std::printf("    // degree %zu: NOT FOUND (residual %.2e)%s\n", degree,
        rule.residual, (seconds > 0) ? " -- deadline reached; retry without one" : "");
      std::fflush(stdout);
      continue;
    }
    const auto report = exactnessSweep(rule, g, degree);
    if (report.worstRelativeError >= 1e-12 || !allWeightsPositive(rule) ||
      !allPointsInside(rule, g, 1e-12))
    {
      std::printf(
        "    // degree %zu: REJECTED (oracle %.2e)\n", degree, report.worstRelativeError);
      std::fflush(stdout);
      continue;
    }
    emit(rule, degree, d, report.worstRelativeError);
  }
  return 0;
}
