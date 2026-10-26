/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Output.h"

#include <limits>
#include <Rodin/Alert.h>

namespace Rodin::Examples::Heart
{
  AtriumFluid::Output::Output(AtriumFluid& atrium, const Context::MPI& context)
    : m_atrium(atrium),
      m_xdmf(context.getCommunicator(), atrium.m_cfg.xdmfBasename)
  {}

  void AtriumFluid::Output::bind()
  {
    m_atrium.m_u.setName("velocity");
    m_atrium.m_p.setName("pressure");
    m_atrium.m_thCur.setName("thrombin");
    m_atrium.m_fgCur.setName("fibrinogen");
    m_atrium.m_fnCur.setName("fibrin");
    m_atrium.m_tawss.setName("TAWSS");
    m_atrium.m_osi.setName("OSI");
    m_atrium.m_activation.setName("activation");
    m_atrium.m_wss.setName("shearStress");

    m_xdmf.setMesh(m_atrium.m_mesh);
    m_xdmf.add("velocity", m_atrium.m_u.getSolution());
    m_xdmf.add("pressure", m_atrium.m_p.getSolution());
    m_xdmf.add("thrombin", m_atrium.m_thCur);
    m_xdmf.add("fibrinogen", m_atrium.m_fgCur);
    m_xdmf.add("fibrin", m_atrium.m_fnCur);
    m_xdmf.add("TAWSS", m_atrium.m_tawss);
    m_xdmf.add("OSI", m_atrium.m_osi);
    m_xdmf.add("activation", m_atrium.m_activation);
    m_xdmf.add("shearStress", m_atrium.m_wss);
  }

  void AtriumFluid::Output::openCSV()
  {
    if (m_atrium.isRoot())
    {
      m_csv.open(m_atrium.m_cfg.csvPath);
      if (!m_csv)
      {
        Alert::Exception() << "Failed to open " << m_atrium.m_cfg.csvPath << Alert::Raise;
      }
      writeCSVHeader();
    }
  }

  void AtriumFluid::Output::writeCSVHeader()
  {
    m_csv << "t,cycle,pIn,pv,par,pd,valveOpen,qIn,qOut,maxU,"
          << "maxTAWSS,maxOSI,maxActivation,maxThrombin,minFibrinogen,maxFibrin\n";
  }

  void AtriumFluid::Output::writeCSVRow(int cycle)
  {
    const auto& s = m_atrium.m_model.getState();
    const Real tawss = m_atrium.m_tawss.max();
    const Real osi = m_atrium.m_osi.max();
    const Real activation = m_atrium.m_activation.max();
    const Real th = m_atrium.m_thCur.max();
    const Real fg = m_atrium.m_fgCur.min();
    const Real fn = m_atrium.m_fnCur.max();

    if (!m_atrium.isRoot())
      return;

    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const Real outletPressure =
      m_atrium.m_cfg.useVentricularModel ? s.pv : m_atrium.m_pOut;
    const Real arterialPressure = m_atrium.m_cfg.useVentricularModel ? s.par : nan;
    const Real distalPressure = m_atrium.m_cfg.useVentricularModel ? s.pd : nan;
    m_csv << m_atrium.m_t << ',' << cycle << ',' << m_atrium.m_pIn << ','
          << outletPressure << ',' << arterialPressure << ',' << distalPressure << ','
          << (m_atrium.m_valveOpen ? 1 : 0) << ',' << m_atrium.m_qIn << ','
          << m_atrium.m_qOut << ',' << m_atrium.m_speed << ',' << tawss << ',' << osi
          << ',' << activation << ',' << th << ',' << fg << ',' << fn << '\n';
    m_csv.flush();
  }

  void AtriumFluid::Output::write(Real time)
  {
    m_xdmf.write(time).flush();
  }

  void AtriumFluid::Output::close()
  {
    m_xdmf.close();
    if (m_atrium.isRoot())
      m_csv.close();
  }
}
