#include "ScalarSolutionS.h"
#include "VectorSolutionS.h"

#include "AdvectS.h"

namespace Rodin::External::MMG
{
  AdvectS::AdvectS(ScalarSolutionS& ls, VectorSolutionS& disp)
    : m_t(0),
      m_avoidTrunc(false),
      m_ex(true),
      m_ls(ls),
      m_disp(disp),
      m_advect(ADVECTION_EXECUTABLE)
  {}

  AdvectS& AdvectS::avoidTimeTruncation(bool avoidTrunc)
  {
    m_avoidTrunc = avoidTrunc;
    return *this;
  }

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

    std::vector<std::string> args{
      meshp,
      "-dt", std::to_string(dt),
      m_avoidTrunc ? "-nocfl" : "",
      m_ex ? "" : "-noex",
      "-c", solp,
      "-s", dispp,
      "-o", outp,
      "-v"
    };

    m_advect.run(args);

    m_ls = ScalarSolutionS::load(outp).setMesh(mesh);

    m_t += dt;
  }
}
