/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Redistancer2D.h"
#include "ScalarSolution2D.h"

namespace Rodin::External::MMG
{
  void Redistancer2D::redistance(ScalarSolution2D<>& ls) const
  {
    auto meshp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    ls.getMesh().save(meshp);

    std::filesystem::path solp(meshp);
    solp.replace_extension(".sol");
    ls.save(solp);

    m_mshdist.run({ solp });

    auto res = ScalarSolution2D<>::load(solp).setMesh(ls.getMesh());
    ls = std::move(res);
  }
}
