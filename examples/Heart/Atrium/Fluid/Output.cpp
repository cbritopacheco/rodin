/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Output.h"

#include <limits>
#include <filesystem>
#include <Rodin/Alert.h>

namespace Rodin::Examples::Heart::Atrium
{
  Fluid::Output::Output(Fluid& atrium, const Context::MPI& context)
    : m_atrium(atrium),
      m_xdmf(context.getCommunicator(), (atrium.m_run.outputDirectory + "/Fluid"))
  {
    std::filesystem::create_directories(atrium.m_run.outputDirectory);
  }

  void Fluid::Output::bind()
  {
    m_atrium.m_u.setName("Velocity");
    m_atrium.m_p.setName("Pressure");
    m_atrium.m_tawss.setName("TimeAveragedWallShearStress");
    m_atrium.m_osi.setName("OscillatoryShearIndex");
    m_atrium.m_wss.setName("WallShearStress");

    m_xdmf.setMesh(m_atrium.m_mesh,
      m_atrium.m_coupled ? IO::XDMF::MeshPolicy::Transient
                         : IO::XDMF::MeshPolicy::Static);
    m_xdmf.add("Velocity", m_atrium.m_u.getSolution());
    m_xdmf.add("Pressure", m_atrium.m_p.getSolution());
    m_xdmf.add("TimeAveragedWallShearStress", m_atrium.m_tawss);
    m_xdmf.add("OscillatoryShearIndex", m_atrium.m_osi);
    m_xdmf.add("WallShearStress", m_atrium.m_wss);
    if (m_atrium.m_coupled)
    {
      m_xdmf.add("MeshDisplacement", m_atrium.m_meshDisplacement);
      m_xdmf.add("MeshVelocity", m_atrium.m_meshVelocity);
    }
    if (m_atrium.m_cfg.solveKinetics)
    {
      m_atrium.m_thCur.setName("Thrombin");
      m_atrium.m_fgCur.setName("Fibrinogen");
      m_atrium.m_fnCur.setName("Fibrin");
      m_atrium.m_activation.setName("ThrombosisActivation");
      m_xdmf.add("Thrombin", m_atrium.m_thCur);
      m_xdmf.add("Fibrinogen", m_atrium.m_fgCur);
      m_xdmf.add("Fibrin", m_atrium.m_fnCur);
      m_xdmf.add("ThrombosisActivation", m_atrium.m_activation);
    }
  }

  void Fluid::Output::collectDiagnostics(int cycle)
  {
    const auto& s = m_atrium.m_model.getState();
    const Real tawss = m_atrium.m_tawss.max();
    const Real osi = m_atrium.m_osi.max();
    const Real activation = m_atrium.m_activation.max();
    const Real th = m_atrium.m_thCur.max();
    const Real fg = m_atrium.m_fgCur.min();
    const Real fn = m_atrium.m_fnCur.max();

    m_atrium.m_diagnostics = Diagnostics{cycle, m_atrium.m_pIn,
      m_atrium.m_cfg.useVentricularModel ? s.pv : m_atrium.m_pOut,
      m_atrium.m_cfg.useVentricularModel ? Optional<Real>(s.par) : Optional<Real>{},
      m_atrium.m_cfg.useVentricularModel ? Optional<Real>(s.pd) : Optional<Real>{},
      m_atrium.m_valveOpen, m_atrium.m_qIn, m_atrium.m_qOut, m_atrium.m_speed, tawss, osi,
      activation, th, fg, fn, m_atrium.m_linearIterations,
      m_atrium.m_lastPreconditionerRefreshStep
        ? m_atrium.m_step - *m_atrium.m_lastPreconditionerRefreshStep
        : 0,
      m_atrium.m_preconditionerRefreshes, m_atrium.m_preconditionerRefreshed,
      m_atrium.m_preconditionerRetried};
  }

  void Fluid::Output::write(Real time)
  {
    m_xdmf.write(time).flush();
  }

  void Fluid::Output::close()
  {
    m_xdmf.close();
  }
}
