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
      m_atrium.m_pCur.setName("Prothrombin");
      m_atrium.m_tCur.setName("Thrombin");
      m_atrium.m_aCur.setName("Antithrombin");
      m_atrium.m_gCur.setName("Fibrinogen");
      m_atrium.m_fibrin.setName("Fibrin");
      m_atrium.m_tat.setName("ThrombinAntithrombin");
      m_atrium.m_drag.setName("BrinkmanDrag");
      m_atrium.m_drift.setName("InvariantDrift");
      m_atrium.m_activation.setName("ThrombosisActivation");
      m_atrium.m_layer.setName("ActivationLayer");
      m_xdmf.add("Prothrombin", m_atrium.m_pCur);
      m_xdmf.add("Thrombin", m_atrium.m_tCur);
      m_xdmf.add("Antithrombin", m_atrium.m_aCur);
      m_xdmf.add("Fibrinogen", m_atrium.m_gCur);
      m_xdmf.add("Fibrin", m_atrium.m_fibrin);
      m_xdmf.add("ThrombinAntithrombin", m_atrium.m_tat);
      m_xdmf.add("BrinkmanDrag", m_atrium.m_drag);
      m_xdmf.add("InvariantDrift", m_atrium.m_drift);
      m_xdmf.add("ThrombosisActivation", m_atrium.m_activation);
      m_xdmf.add("ActivationLayer", m_atrium.m_layer);
    }
  }

  void Fluid::Output::collectDiagnostics(int cycle)
  {
    const auto& s = m_atrium.m_model.getState();
    const Real tawss = m_atrium.m_tawss.max();
    const Real osi = m_atrium.m_osi.max();
    const Real activation = m_atrium.m_activation.max();
    const Real th = m_atrium.m_tCur.max();
    const Real fg = m_atrium.m_gCur.min();
    const Real fn = m_atrium.m_fibrin.max();

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
