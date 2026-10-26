/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/// @file @brief Published samples remain readable during static destruction.
#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "Rodin/QF/Centroid.h"
#include "Rodin/QF/GrundmannMoller.h"
#include "Rodin/QF/WitherdenVincent.h"
#include "Rodin/QF/XiaoGimbutas.h"

using namespace Rodin;

namespace
{
  const QF::QuadratureFormula* formula = nullptr;
  Real expected = 0;

  /// @brief Reads published samples after ordinary static coefficient tables are destroyed.
  struct FinalReader
  {
    ~FinalReader()
    {
      if (formula)
      {
        const Real actual = formula->getWeight(0);
        if (!std::isfinite(actual) || actual != expected)
        {
          std::fprintf(
            stderr, "Published quadrature coefficients expired at shutdown.\n");
          std::abort();
        }
      }
    }
  };

#if defined(__GNUC__) || defined(__clang__)
  // Register before the ordinary-priority coefficient tables. Destruction is
  // reversed, so ASan observes any sample storage released before this reader.
  FinalReader reader __attribute__((init_priority(101)));
#else
  FinalReader reader;
#endif

  template <class Family, class... Arguments>
  void retain(Arguments... arguments)
  {
    formula = &QF::QuadratureFormula::get<Family>(arguments...);
    expected = formula->getWeight(0);
    EXPECT_GT(expected, 0);
  }
}

TEST(QuadratureLifetime, XiaoGimbutas)
{
  retain<QF::XiaoGimbutas>(size_t(1), Geometry::Polytope::Type::Triangle);
}

TEST(QuadratureLifetime, WitherdenVincent)
{
  retain<QF::WitherdenVincent>(size_t(1), Geometry::Polytope::Type::Triangle);
}

TEST(QuadratureLifetime, GrundmannMoller)
{
  retain<QF::GrundmannMoller>(size_t(0), Geometry::Polytope::Type::Triangle);
}

TEST(QuadratureLifetime, Centroid)
{
  retain<QF::Centroid>(Geometry::Polytope::Type::Triangle);
}

/// @brief Borrowed tensor-product weights remain readable during shutdown.
TEST(QuadratureLifetime, WitherdenVincentQuadrilateral)
{
  retain<QF::WitherdenVincent>(size_t(1), Geometry::Polytope::Type::Quadrilateral);
}

/// @brief Borrowed wedge weights remain readable during shutdown.
TEST(QuadratureLifetime, WitherdenVincentWedge)
{
  retain<QF::WitherdenVincent>(size_t(1), Geometry::Polytope::Type::Wedge);
}

/// @brief Borrowed hexahedron weights remain readable during shutdown.
TEST(QuadratureLifetime, WitherdenVincentHexahedron)
{
  retain<QF::WitherdenVincent>(size_t(1), Geometry::Polytope::Type::Hexahedron);
}
