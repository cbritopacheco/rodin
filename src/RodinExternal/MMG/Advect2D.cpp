#include <cmath>

#include "ScalarSolution2D.h"
#include "VectorSolution2D.h"

#include "Advect2D.h"

namespace Rodin::External::MMG
{
  Advect2D::Advect2D(ScalarSolution2D& ls, VectorSolution2D& disp)
    : m_t(0),
      m_avoidTrunc(false),
      m_ex(true),
      m_ls(ls),
      m_disp(disp),
      m_advect(ADVECTION_EXECUTABLE)
  {}

  Advect2D& Advect2D::avoidTimeTruncation(bool avoidTrunc)
  {
    m_avoidTrunc = avoidTrunc;
    return *this;
  }

  Advect2D& Advect2D::enableExtrapolation(bool ex)
  {
    m_ex = ex;
    return *this;
  }

  void Advect2D::step(double dt)
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
      meshp.string(),
      "-dt", std::to_string(dt),
      m_avoidTrunc ? "-nocfl" : "",
      m_ex ? "" : "-noex",
      "-c", solp.string(),
      "-s", dispp.string(),
      "-o", outp.string(),
      "-v"
    };

    int retcode = m_advect.run(args);
    if (retcode != 0)
      Alert::Exception("ISCD::Advection invocation failed.").raise();

    m_ls = ScalarSolution2D::load(outp).setMesh(mesh);

    m_t += dt;
  }
}
