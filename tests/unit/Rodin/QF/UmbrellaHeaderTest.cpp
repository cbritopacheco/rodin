/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */

/**
 * @file
 * @brief Checks that the QF umbrella header exposes what it advertises.
 *
 * Nothing else in the tree includes @c Rodin/QF.h, so without this test the
 * umbrella is never compiled and a formula can be added to the module without
 * ever becoming reachable through the documented entry point -- which is how
 * the three vendored families came to be missing from it. This file therefore
 * includes the umbrella and nothing else on purpose.
 */
#include <Rodin/QF.h>

#include <gtest/gtest.h>

using namespace Rodin;
using namespace Rodin::Geometry;

/// @brief Every formula the umbrella documents can be named through it.
TEST(QFUmbrellaHeaderTest, AdvertisedFormulaeAreReachable)
{
  const QF::Centroid centroid(Polytope::Type::Triangle);
  EXPECT_EQ(centroid.getSize(), 1u);

  const QF::GaussLegendre line(Polytope::Type::Segment, 2);
  EXPECT_EQ(line.getSize(), 2u);

  const QF::GrundmannMoller gm(1, Polytope::Type::Triangle);
  EXPECT_GT(gm.getSize(), 0u);

  const QF::XiaoGimbutas xg(3, Polytope::Type::Triangle);
  EXPECT_GT(xg.getSize(), 0u);

  const QF::WitherdenVincent wv(3, Polytope::Type::Hexahedron);
  EXPECT_GT(wv.getSize(), 0u);

  const QF::TensorProduct tensor(Polytope::Type::Quadrilateral, line, line);
  EXPECT_EQ(tensor.getSize(), line.getSize() * line.getSize());

  const auto& dispatched =
    QF::PolytopeQuadratureFormula::get(3, Polytope::Type::Tetrahedron);
  EXPECT_GT(dispatched.getSize(), 0u);
}
