/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <utility>
#include "Distancer2D.h"
#include "ScalarSolution2D.h"

namespace Rodin::External::MMG
{
  ScalarSolution2D<> Distancer2D::distance(Mesh2D& box, Mesh2D& contour) const
  {
    auto boxp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    box.save(boxp);

    auto contourp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    contour.save(contourp);

    m_mshdist.run({ boxp.string(), contourp.string() });

    auto res = ScalarSolution2D<>::load(boxp.replace_extension(".sol")).setMesh(box);
    return res;
  }
}
