/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <thread>
#include <fstream>
#include <utility>

#include "Distancer2D.h"
#include "ScalarSolution2D.h"

namespace Rodin::External::MMG
{
  Distancer2D::Distancer2D()
    : m_scale(true),
      m_ncpu(std::thread::hardware_concurrency()),
      m_mshdist(MSHDIST_EXECUTABLE)
  {}

  Distancer2D& Distancer2D::setInteriorDomain(
      const std::set<MaterialReference>& refs)
  {
    m_interiorDomains = refs;
    return *this;
  }

  ScalarSolution2D Distancer2D::distance(Mesh2D& box)
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
        "-dom",
        "-ncpu", m_ncpu,
        "-v 0");

    if (retcode != 0)
      Alert::Exception("ISCD::Mshdist invocation failed.").raise();

    auto res = ScalarSolution2D::load(boxp.replace_extension(".sol")).setMesh(box);
    return res;
  }

  ScalarSolution2D Distancer2D::distance(Mesh2D& box, Mesh2D& contour)
  {
    auto boxp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    box.save(boxp);

    auto contourp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    contour.save(contourp);

    int retcode;
    if (m_scale)
      retcode = m_mshdist.run(boxp.string(), contourp.string(), "-ncpu", m_ncpu, "-v 0");
    else
      retcode = m_mshdist.run(boxp.string(), contourp.string(), "-noscale", "-ncpu", m_ncpu, "-v 0");

    if (retcode != 0)
      Alert::Exception("ISCD::Mshdist invocation failed.").raise();

    auto res = ScalarSolution2D::load(boxp.replace_extension(".sol")).setMesh(box);
    return res;
  }

  void Distancer2D::redistance(ScalarSolution2D& ls)
  {
    auto meshp = m_mshdist.tmpnam(".mesh", "RodinMMG");
    ls.getMesh().save(meshp);

    boost::filesystem::path solp(meshp);
    solp.replace_extension(".sol");
    ls.save(solp);

    auto name = solp;
    name.replace_extension();

    int retcode = m_mshdist.run(name.string(), "-ncpu", m_ncpu, "-v 0");
    if (retcode != 0)
      Alert::Exception("ISCD::Mshdist invocation failed.").raise();

    auto res = ScalarSolution2D::load(solp).setMesh(ls.getMesh());
    ls = std::move(res);
  }

  Distancer2D& Distancer2D::setCPUs(unsigned int ncpu)
  {
    m_ncpu = ncpu;
    return *this;
  }

  Distancer2D& Distancer2D::enableScaling(bool b)
  {
    m_scale = b;
    return *this;
  }

  unsigned int Distancer2D::getCPUs() const
  {
    return m_ncpu;
  }
}
