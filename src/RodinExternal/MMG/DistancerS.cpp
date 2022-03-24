/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <thread>
#include <fstream>
#include <utility>

#include "DistancerS.h"
#include "ScalarSolutionS.h"

namespace Rodin::External::MMG
{
  DistancerS::DistancerS()
    : m_scale(true),
      m_activeBorder(false),
      m_ncpu(std::thread::hardware_concurrency()),
      m_mshdist(MSHDIST_EXECUTABLE)
  {}

  DistancerS& DistancerS::setInteriorDomain(
      const std::set<MaterialReference>& refs)
  {
    m_interiorDomains = refs;
    return *this;
  }

  DistancerS& DistancerS::setActiveBorders()
  {
    m_activeBorder = true;
    return *this;
  }

  DistancerS& DistancerS::setActiveBorder(const std::set<MaterialReference>& refs)
  {
    assert(refs.size() > 0);
    m_activeBorders = refs;
    m_activeBorder = true;
    return *this;
  }

  ScalarSolutionS DistancerS::distance(MeshS& box)
  {
    auto boxp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    box.save(boxp);

    if (m_interiorDomains)
    {
      auto paramp(boxp);
      paramp.replace_extension(".mshdist");
      std::ofstream paramf(paramp, std::ios::trunc);
      if (m_interiorDomains->size() > 0)
      {
        paramf << "InteriorDomains\n"
               << m_interiorDomains->size() << "\n\n";
        for (const auto& ref : *m_interiorDomains)
          paramf << ref << "\n";
      }
      if (m_activeBorder && m_activeBorders->size() > 0)
      {
        paramf << "ActiveBorders\n"
               << m_activeBorders->size() << "\n\n";
        for (const auto& ref : *m_activeBorders)
          paramf << ref << "\n";
      }
    }

    m_mshdist.run(
        boxp,
        "-dom",
        m_activeBorder ? "-activeBorder" : "",
        "-ncpu", m_ncpu,
        "-v 0");

    auto res = ScalarSolutionS::load(boxp.replace_extension(".sol")).setMesh(box);
    return res;
  }

  ScalarSolutionS DistancerS::distance(MeshS& box, MeshS& contour)
  {
    auto boxp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    box.save(boxp);

    auto contourp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    contour.save(contourp);

    if (m_scale)
      m_mshdist.run(boxp, contourp, "-ncpu", m_ncpu, "-v 0");
    else
      m_mshdist.run(boxp, contourp, "-noscale", "-ncpu", m_ncpu, "-v 0");

    auto res = ScalarSolutionS::load(boxp.replace_extension(".sol")).setMesh(box);
    return res;
  }

  void DistancerS::redistance(ScalarSolutionS& ls)
  {
    auto meshp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    ls.getMesh().save(meshp);

    std::filesystem::path solp(meshp);
    solp.replace_extension(".sol");
    ls.save(solp);

    m_mshdist.run(solp , "-ncpu", m_ncpu, "-v 0");

    auto res = ScalarSolutionS::load(solp).setMesh(ls.getMesh());
    ls = std::move(res);
  }

  DistancerS& DistancerS::setCPUs(unsigned int ncpu)
  {
    m_ncpu = ncpu;
    return *this;
  }

  DistancerS& DistancerS::enableScaling(bool b)
  {
    m_scale = b;
    return *this;
  }

  unsigned int DistancerS::getCPUs() const
  {
    return m_ncpu;
  }
}

