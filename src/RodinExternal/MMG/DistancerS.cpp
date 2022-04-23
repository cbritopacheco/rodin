/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <thread>
#include <fstream>
#include <utility>

#include "Rodin/Alert.h"

#include "DistancerS.h"
#include "ScalarSolutionS.h"

namespace Rodin::External::MMG
{
  DistancerS::DistancerS()
    : m_scale(true),
      m_ncpu(std::thread::hardware_concurrency()),
      m_mshdist(MSHDIST_EXECUTABLE)
  {}

  DistancerS& DistancerS::setInteriorDomain(
      const std::set<MaterialReference>& refs)
  {
    m_interiorDomains = refs;
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
      std::ofstream paramf(paramp.string(), std::ios::trunc);
      if (m_interiorDomains->size() > 0)
      {
        paramf << "InteriorDomains\n"
               << m_interiorDomains->size() << "\n\n";
        for (const auto& ref : *m_interiorDomains)
          paramf << ref << "\n";
      }
    }

    int retcode = m_mshdist.run(
        boxp.string(),
        "-fmm",
        "-surf",
        "-dom",
        "-ncpu", m_ncpu,
        "-v 0");

    if (retcode != 0)
      Alert::Exception("ISCD::Mshdist invocation failed.").raise();

    ScalarSolutionS res(box);
    res.load(boxp.replace_extension(".sol"));
    return res;
  }

  ScalarSolutionS DistancerS::distance(MeshS& box, MeshS& contour)
  {
    auto boxp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    box.save(boxp);

    auto contourp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    contour.save(contourp);

    int retcode;
    if (m_scale)
      retcode = m_mshdist.run(
          boxp.string(), contourp.string(), "-surf", "-ncpu", m_ncpu, "-v 0");
    else
      retcode = m_mshdist.run(
          boxp.string(), contourp.string(), "-surf", "-noscale", "-ncpu", m_ncpu, "-v 0");

    if (retcode != 0)
      Alert::Exception("ISCD::Mshdist invocation failed.").raise();

    ScalarSolutionS res(box);
    res.load(boxp.replace_extension(".sol"));
    return res;
  }

  void DistancerS::redistance(ScalarSolutionS& ls)
  {
    auto meshp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    ls.getMesh().save(meshp);

    boost::filesystem::path solp(meshp);
    solp.replace_extension(".sol");
    ls.save(solp);

    auto name = solp;
    name.replace_extension();

    int retcode = m_mshdist.run(name.string(), "-fmm", "-surf", "-ncpu", m_ncpu, "-v 0");
    if (retcode != 0)
      Alert::Exception("ISCD::Mshdist invocation failed.").raise();

    ScalarSolutionS res(ls.getMesh());
    res.load(solp);
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

