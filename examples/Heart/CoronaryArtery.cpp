/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file CoupledLV0D_Coronary3D.cpp
 * @brief One-way coupled LV-0D to 3D coronary Navier–Stokes with RCR outlets.
 *
 * Coupling structure
 * ------------------
 *   prescribed atrial pressure p_at(t)
 *          -> LV 0D model
 *          -> arterial pressure p_ar(t)
 *          -> 3D coronary inlet pressure traction
 *          -> one RCR Windkessel per 3D outlet
 *
 * Notes
 * -----
 * 1) This is a partitioned one-way coupling:
 *      LV-0D -> 3D coronary.
 *    The 3D coronary flow does NOT yet feed back into the 0D arterial balance.
 *
 * 2) The 3D model is Newtonian and uses the same Oseen / Picard-style
 *    semi-implicit linearization pattern as your current working example.
 *
 * 3) Outlet RCR parameters below are placeholders and should be calibrated.
 *
 * Suggested PETSc options
 * -----------------------
 * Debug / direct:
 *   -ksp_type preonly
 *   -pc_type lu
 *   -pc_factor_mat_solver_type mumps
 *   -pc_factor_shift_type nonzero
 *   -pc_factor_shift_amount 1e-10
 *   -ksp_converged_reason
 *
 * Iterative:
 *   -ksp_type fgmres
 *   -pc_type fieldsplit
 *   -pc_fieldsplit_type schur
 *   -pc_fieldsplit_schur_fact_type lower
 *   -pc_fieldsplit_schur_precondition selfp
 *   -fieldsplit_u_ksp_type preonly
 *   -fieldsplit_u_pc_type gamg
 *   -fieldsplit_p_ksp_type preonly
 *   -fieldsplit_p_pc_type jacobi
 *   -ksp_rtol 1e-6
 *   -ksp_monitor
 *   -ksp_converged_reason
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <numbers>
#include <type_traits>
#include <algorithm>
#include <array>
#include <map>

#include "Rodin/Heart/CCMLC2014.h"

#include "Rodin/Geometry/Types.h"
#include "Rodin/IO/ForwardDecls.h"
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>

using Real = Rodin::Real;
using Model = Rodin::Heart::CCMLC2014T<>;

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  struct RCR
  {
    Real Rp;
    Real C;
    Real Rd;
    Real pd;
    Real pc;
    Real pout;
  };

  static void updateRCR(RCR& bc, Real Q, Real dt)
  {
    // Backward Euler update of:
    //   C (pc^{n+1} - pc^n)/dt + (pc^{n+1} - pd)/Rd = Q^{n+1}
    //
    // and
    //   pout^{n+1} = pc^{n+1} + Rp * Q^{n+1}
    //
    // Sign convention:
    //   Q = ∫_{Γ_out} u·n
    // with n outward, so Q > 0 means outflow.

    const Real a = bc.C / dt;
    bc.pc = (a * bc.pc + Q + bc.pd / bc.Rd) / (a + 1.0 / bc.Rd);
    bc.pout = bc.pc + bc.Rp * Q;
  }

  static Real periodic_activation(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);

    if (tau < 0.13)  return 0.0;
    if (tau < 0.141) return 35.0 * ((tau - 0.13) / 0.011);
    if (tau < 0.281) return 35.0;
    if (tau < 0.361) return 35.0 - 47.0 * ((tau - 0.281) / 0.08);
    if (tau < 0.45)  return -12.0;
    return 0.0;
  }

  static Real atrial_pressure(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);

    const Real min_value = 500.0;
    const Real max_value = 1000.0;
    const Real second_threshold = 1250.0;

    const Real t1 = 0.02;
    const Real t2 = 0.15;
    const Real t3 = 0.17;
    const Real t4 = 0.56;
    const Real t5 = 0.62;
    const Real t6 = 0.85;

    Real alpha = 0.0;
    Real value = min_value;

    if (tau < t1)
    {
      alpha = -(tau - t1) / t1;
      value = alpha * min_value + (1.0 - alpha) * max_value;
    }
    else if (tau < t2)
    {
      value = max_value;
    }
    else if (tau < t3)
    {
      alpha = -(tau - t3) / (t3 - t2);
      value = alpha * max_value + (1.0 - alpha) * min_value;
    }
    else if (tau < t4)
    {
      alpha = -(tau - t4) / (t4 - t3);
      value = alpha * min_value + (1.0 - alpha) * second_threshold;
    }
    else if (tau < t5)
    {
      value = second_threshold;
    }
    else if (tau < t6)
    {
      alpha = -(tau - t6) / (t6 - t5);
      value = alpha * second_threshold + (1.0 - alpha) * min_value;
    }
    else
    {
      value = min_value;
    }

    return value;
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  try
  {
    //--------------------------------------------------------------------------
    // 0D LV model setup
    //--------------------------------------------------------------------------

    Model::Input in;

    // Geometry / inertia
    in.rho = 1.0e3;
    in.R0 = 2.36e-2;
    in.d0 = 1.42e-2;

    // Active law parameters
    in.Es = 3.0e5;
    in.mu = 70.0;
    in.eta = 70.0;
    in.alpha = 3.0;
    in.k0 = 1.0e5;
    in.sigma0 = 5.0e5;

    // 0D Windkessel inside the LV model
    in.Rp = 8.0e6;
    in.Cp = 5.0e-9;
    in.Rd = 1.0e8;
    in.Cd = 1.0e-8;

    // Valve parameters
    in.Kat = 8.0e-7;
    in.Kp  = 5.0e-10;
    in.Kar = 1.3e-5;

    in.cavityCapacity = 5.0e-12;

    // Local active solve
    in.localTolerance = 1e-12;
    in.localMaxIterations = 50;
    in.localDamping = 1.0;
    in.absRegularization = 1e-14;

    // Initial local active state
    in.initFibDef = 0.0;
    in.initActiveStiffness = 0.0;
    in.initActiveStress = 0.0;

    in.pSv = [](Real) { return 1.0e3; };
    in.pAt = atrial_pressure;
    in.u = periodic_activation;

    {
      using PassiveEnergy = std::decay_t<decltype(in.passiveEnergy)>;
      typename PassiveEnergy::Parameters hp;
      hp.mu1 = 0.0;
      hp.mu2 = 0.0;
      hp.C0 = 1.9e3;
      hp.C1 = 1.1e-1;
      hp.C2 = 1.9e3;
      hp.C3 = 1.1e-1;
      in.passiveEnergy = PassiveEnergy(hp);
    }

    Model model(in);
    model.setMaxIterations(200)
         .setAbsoluteTolerance(1e-8)
         .setRelativeTolerance(1e-8)
         .setStepTolerance(1e-10)
         .setDampingFactor(1.0);

    Model::State s0;
    s0.t = 0.0;

    s0.y   = 0.0;
    s0.v   = 0.0;
    s0.pv  = in.pAt(0.0) - 100.0;
    s0.par = 11000.0;
    s0.pd  = 10000.0;

    s0.ec = in.initFibDef;
    s0.gamma = std::sqrt(std::max<Real>(in.initActiveStiffness, 0.0));
    s0.beta = (s0.gamma > 0.0) ? (in.initActiveStress / s0.gamma) : 0.0;
    s0.kc = s0.gamma * s0.gamma;
    s0.tauc = s0.gamma * s0.beta;

    model.initialize(s0);

    {
      const auto& s = model.getState();
      std::cout << "Initial 0D state:\n"
                << "  y     = " << s.y << '\n'
                << "  v     = " << s.v << '\n'
                << "  pv    = " << s.pv << '\n'
                << "  par   = " << s.par << '\n'
                << "  pd    = " << s.pd << '\n'
                << "  ec    = " << s.ec << '\n'
                << "  gamma = " << s.gamma << '\n'
                << "  beta  = " << s.beta << '\n'
                << "  kc    = " << s.kc << '\n'
                << "  tauc  = " << s.tauc << '\n';
    }

    //--------------------------------------------------------------------------
    // 3D coronary model setup
    //--------------------------------------------------------------------------

    const Attribute wall = 2;
    const Attribute inlet = 3;
    const std::array<Attribute, 6> outlets{{4, 5, 6, 7, 8, 9}};

    const Real eps = 1e-12;

    // Newtonian fluid parameters for the 3D coronary solve
    const Real rho = 1.0;
    const Real mu  = 0.01;

    // Shared time step with the 0D model
    const Real dt = 1e-3;
    const int nsteps = 3 * static_cast<int>(0.85 / dt);

    Mesh mesh;
    mesh.load("CoronaryArtery.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "Computing connectivity for CoronaryArtery.mesh ..."
                  << Alert::Raise;

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(2, 3);

    Alert::Info() << "Setting up CoronaryArtery.xdmf ..."
                  << Alert::Raise;

    IO::XDMF xdmf("CoronaryArtery");
    xdmf.setMesh(mesh);

    const size_t dim = mesh.getSpaceDimension();
    if (dim != 3)
      throw std::runtime_error("Expected a 3D coronary mesh.");

    H1 uh(std::integral_constant<size_t, 2>{}, mesh, dim);
    H1 ph(std::integral_constant<size_t, 1>{}, mesh);

    PETSc::Variational::TrialFunction u(uh); u.setName("u");
    PETSc::Variational::TrialFunction p(ph); p.setName("p");
    PETSc::Variational::TestFunction  v(uh);
    PETSc::Variational::TestFunction  q(ph);

    PETSc::Variational::GridFunction u_old(uh);
    PETSc::Variational::GridFunction p_old(ph);

    xdmf.add("velocity", u.getSolution());
    xdmf.add("pressure", p.getSolution());

    u_old = Math::Vector<Real>{{0.0, 0.0, 0.0}};
    p_old = 0.0;

    const auto n = BoundaryNormal(mesh);

    // Helper objects for flux computation
    PETSc::Variational::GridFunction one(ph);
    one = 1.0;
    PETSc::Variational::TestFunction qFlux(ph);
    LinearForm flux(qFlux);

    // One RCR Windkessel per outlet
    //
    // Placeholder parameters. Calibrate later.
    std::map<Attribute, RCR> wk;
    for (const Attribute tag : outlets)
    {
      RCR bc;
      bc.Rp   = 100.0;
      bc.C    = 1.0;
      bc.Rd   = 300.0;
      bc.pd   = 0.0;
      bc.pc   = bc.pd;
      bc.pout = bc.pd;
      wk.emplace(tag, bc);
    }

    //--------------------------------------------------------------------------
    // Diagnostics
    //--------------------------------------------------------------------------

    std::ofstream out0d("ccmlc2014_0d_cycle.csv");
    out0d << "t,y,v,pv,par,pd,ec,gamma,beta,kc,tauc,V,Q,pat\n";

    std::ofstream inletPressureFile("coronary_inlet_pressure.txt");
    std::ofstream inletFluxFile("coronary_inlet_flux.txt");

    std::map<Attribute, std::ofstream> outletFluxFiles;
    std::map<Attribute, std::ofstream> outletPressureFiles;
    for (const Attribute tag : outlets)
    {
      outletFluxFiles.emplace(tag, std::ofstream("coronary_flux_outlet_" + std::to_string(tag) + ".txt"));
      outletPressureFiles.emplace(tag, std::ofstream("coronary_pressure_outlet_" + std::to_string(tag) + ".txt"));
    }

    //--------------------------------------------------------------------------
    // Shared 0D / 3D time loop
    //--------------------------------------------------------------------------

    for (int i = 0; i < nsteps; ++i)
    {
      std::cout << "Step " << i << ": t = " << model.getState().t << "\n";

      //----------------------------------------------------------------------
      // 1) Advance the 0D LV model
      //----------------------------------------------------------------------

      const auto rep = model.step(dt);
      std::cout << "  0D Newton step: "
                << (rep.converged ? "converged" : "not converged")
                << ", iterations = " << rep.iterations
                << ", final residual = " << rep.finalResidual
                << ", final step norm = " << rep.finalStepNorm
                << '\n';

      if (!rep.converged)
      {
        std::cerr << "0D solver failed to converge at step "
                  << i << ", t = " << model.getState().t << "\n";
        break;
      }

      const auto& s = model.getState();

      const Real R = in.R0 + s.y;
      const Real V = (4.0 / 3.0) * std::numbers::pi_v<Real> * R * R * R;
      const Real Qlv = 4.0 * std::numbers::pi_v<Real> * R * R * s.v;
      const Real pat = in.pAt(s.t);

      out0d << s.t << ","
            << s.y << ","
            << s.v << ","
            << s.pv << ","
            << s.par << ","
            << s.pd << ","
            << s.ec << ","
            << s.gamma << ","
            << s.beta << ","
            << s.kc << ","
            << s.tauc << ","
            << V << ","
            << Qlv << ","
            << pat << "\n";

      //----------------------------------------------------------------------
      // 2) Use 0D arterial pressure as the 3D coronary inlet pressure
      //----------------------------------------------------------------------

      const Real pin = s.par;

      const auto conv_u = Mult(Jacobian(u), u_old);
      const auto div_u_old = Div(u_old);
      const auto beta = Max(-Dot(u_old, n), 0.0);

      Problem flow(u, p, v, q);
      flow =
          (rho / dt) * Integral(u, v)
        - (rho / dt) * Integral(u_old, v)
        + rho * Integral(Dot(conv_u, v))
        + 0.5 * rho * Integral(div_u_old * Dot(u, v))
        + mu * Integral(Jacobian(u), Jacobian(v))
        - Integral(p, Div(v))
        + Integral(Div(u), q)
        + eps * Integral(p, q)

        // Inlet pressure traction from the 0D arterial pressure
        + BoundaryIntegral(pin * Dot(v, n)).over(inlet)

        // One RCR pressure traction per outlet
        - BoundaryIntegral(wk[4].pout * Dot(v, n)).over(4)
        - BoundaryIntegral(wk[5].pout * Dot(v, n)).over(5)
        - BoundaryIntegral(wk[6].pout * Dot(v, n)).over(6)
        - BoundaryIntegral(wk[7].pout * Dot(v, n)).over(7)
        - BoundaryIntegral(wk[8].pout * Dot(v, n)).over(8)
        - BoundaryIntegral(wk[9].pout * Dot(v, n)).over(9)

        // Outlet backflow stabilization
        + BoundaryIntegral(0.5 * rho * beta * Dot(u, v)).over(4)
        + BoundaryIntegral(0.5 * rho * beta * Dot(u, v)).over(5)
        + BoundaryIntegral(0.5 * rho * beta * Dot(u, v)).over(6)
        + BoundaryIntegral(0.5 * rho * beta * Dot(u, v)).over(7)
        + BoundaryIntegral(0.5 * rho * beta * Dot(u, v)).over(8)
        + BoundaryIntegral(0.5 * rho * beta * Dot(u, v)).over(9)

        + DirichletBC(u, Zero(dim)).on(wall);

      Alert::Info() << "Assembling 3D time step " << i + 1 << " / " << nsteps << "..."
                    << Alert::Raise;

      flow.assemble().setFieldSplits();

      Alert::Info() << "Solving 3D time step " << i + 1 << " / " << nsteps << "..."
                    << Alert::Raise;

      Solver::KSP(flow).solve();

      //----------------------------------------------------------------------
      // 3) Diagnostics: inlet flux
      //----------------------------------------------------------------------

      flux = BoundaryIntegral(Dot(u.getSolution(), n), qFlux).over(inlet);
      flux.assemble();
      const Real Qin = flux(one);

      inletPressureFile << s.t << " " << pin << "\n";
      inletFluxFile << s.t << " " << Qin << "\n";

      //----------------------------------------------------------------------
      // 4) Update outlet Windkessels from outlet fluxes
      //----------------------------------------------------------------------

      for (const Attribute tag : outlets)
      {
        flux = BoundaryIntegral(Dot(u.getSolution(), n), qFlux).over(tag);
        flux.assemble();
        const Real Qout = flux(one);

        updateRCR(wk[tag], Qout, dt);

        outletFluxFiles[tag] << s.t << " " << Qout << "\n";
        outletPressureFiles[tag] << s.t << " " << wk[tag].pout << "\n";
      }

      //----------------------------------------------------------------------
      // 5) Advance 3D history and write output
      //----------------------------------------------------------------------

      u_old.setData(u.getSolution().getData());
      p_old.setData(p.getSolution().getData());

      xdmf.write(s.t).flush();
    }

    xdmf.close();
  }
  catch (const std::exception& e)
  {
    std::cerr << "Fatal error: " << e.what() << "\n";
    PetscFinalize();
    return 1;
  }

  PetscFinalize();
  return 0;
}
