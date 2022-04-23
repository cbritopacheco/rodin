#include "Rodin/Alert.h"

#include "ScalarSolutionS.h"
#include "VectorSolutionS.h"

#include "AdvectS.h"

namespace Rodin::External::MMG
{
  AdvectS::AdvectS(ScalarSolutionS& ls, VectorSolutionS& disp)
    : m_t(0),
      m_ex(true),
      m_ls(ls),
      m_disp(disp),
      m_advect(ADVECTION_EXECUTABLE)
  {}

  AdvectS& AdvectS::enableExtrapolation(bool ex)
  {
    m_ex = ex;
    return *this;
  }

  void AdvectS::step(double dt)
  {
    assert(!std::isnan(dt) && !std::isinf(dt));

    auto& mesh = m_ls.getMesh();

    auto meshp = m_advect.tmpnam(".mesh", "RodinMMG");
    mesh.save(meshp);

    auto solp = m_advect.tmpnam(".sol", "RodinMMG");
    m_ls.save(solp);

    auto dispp = m_advect.tmpnam(".sol", "RodinMMG");
    m_disp.save(dispp);

    auto outp = m_advect.tmpnam(".sol", "RodinMMG");

    int retcode = m_advect.run(
      meshp.string(),
      "-dt", std::to_string(dt),
      m_ex ? "" : "-noex",
      "-c", solp.string(),
      "-s", dispp.string(),
      "-o", outp.string(),
      "-nocfl");

    if (retcode != 0)
      Alert::Exception("ISCD::Avection invocation failed.").raise();

    m_ls = ScalarSolutionS(mesh).load(outp);
    m_t += dt;
  }
}
