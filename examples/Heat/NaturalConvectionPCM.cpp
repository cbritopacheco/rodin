// NaturalConvectionPCM.cpp
//
// Melting/solidification of a phase change material (PCM) with natural
// convection in the liquid, solved monolithically in (u, p, T) with Newton
// (PETSc SNES) at every backward-Euler step. Built on ForcedConvectionNewton:
// same Newton/OSGS machinery, plus the enthalpy-porosity model of
// Voller & Prakash (1987) as written in Elsevier_Felipe.pdf, Sec. 1:
//
//   div u = 0,
//   rho_l (du/dt + (grad u) u) - mu lap u + grad p
//       + A(Y) u = -rho_l beta (T - T_ref) g,                    (Eq. 2)
//   dH(T)/dt + H'(T) u.grad T - div(kappa(T) grad T) = 0,        (Eq. 3)
//
// with
//   Y(T)   = clamp((T - T_s)/(T_l - T_s), 0, 1),                 (Eq. 4)
//   H(T)   = int_{T_s}^T rho c_p dT' + rho(T) Y(T) L,            (Eq. 5)
//   Psi(T) = Psi_s + (Psi_l - Psi_s) Y(T),  Psi in {rho, kappa, c_p},  (Eq. 6)
//   A(Y)   = C (1 - Y)^2 / (Y^3 + eps)    (Carman-Kozeny).
//
// Eq. (3) is the paper's rho c_p (dT/dt + u.grad T) + rho L (dY/dt + u.grad Y)
// written in conservative enthalpy form: the time derivative is
// (H(T^{n+1}) - H(T^n))/dt, so latent heat is never skipped when a node jumps
// across the mushy interval within one step. Momentum is Boussinesq: rho_l in
// inertia and buoyancy (u ~ 0 wherever rho != rho_l, by the Darcy term);
// mu (grad u + grad u^T) reduces to mu lap u because div u = 0 and every wall
// is no-slip.
//
// Enthalpy rate with lumped mass: sum_i m_i (H(T_i) - H(T_i^n))/dt, m_i =
// (1, phi_i), added nodally to the T rows of the residual and m_i H'(T_i)/dt to
// the diagonal of the Jacobian (wrapping Rodin's SNES callbacks). Lumping makes
// that Jacobian exact and the capacity matrix diagonal; with a consistent mass
// H' jumps by ~80x between neighbouring nodes at the front and Newton cycles.
//
// Other nonlinear coefficients. H', kappa, kappa', A, A' are P1 interpolants
// of the current Newton iterate (group finite elements), refreshed by the SNES
// state update before each assembly.
//
// Phase-change truncation. H is convex at T_s and concave at T_l, so a plain
// Newton step taken from the solid (slope rho c_p, no latent heat) overshoots
// to the liquid, and from there back to the solid: the iteration cycles and no
// monotone line search can fix it. A line-search post-check moves every node
// whose update crosses the whole interval [T_s, T_l] to T_m, where the slope
// carries the latent heat; from there the next step is exact while the node
// stays mushy. Full (basic) steps are the default so the truncated iterate is
// the one evaluated.
//
// Stabilisation: split OSGS exactly as ForcedConvectionNewton, with the Darcy
// coefficient in tau_1 (Brinkman-type) and the lagged sensible capacity in the
// energy term:
//   tau_1 = [4 nu/h^2 + 2|u^n|/h + A(T^n)/rho_l]^{-1},  tau_2 = rho_l h^2/(4 tau_1),
//   tau_T = [4 alpha(T^n)/h^2 + 2|u^n|/h]^{-1},          weight rho c_p(T^n) tau_T.
//
// Boundary conditions (Eqs. 9-12 of the paper), outward flux q_n = -kappa dT/dn:
//   top    : q_n = h_a(t) (T - T_a) + sigma eps (T^4 - T_sky^4) - alpha_s I(t),
//            T_sky = 0.0552 T_a^1.5 (Swinbank), h_a = 8.91 + 2 w (Loveday-Taki);
//   bottom : q_n = h_in (T - T_in)   (convection only, indoor side);
//   sides  : adiabatic;  u = 0 on every wall.
// Solar irradiance I enters with a + sign as a gain (Eq. 9 in the PDF writes
// it on the loss side; the sign there should be checked).
//
// Inclination. The mesh stays a horizontal rectangle; the roof slope theta
// enters only through gravity in the local frame,
//   g = -|g| (sin theta, cos theta),
// so theta = 0 is a flat roof heated from above.
//
// Roof test (default configuration): 2 m x 0.10 m layer of RT42 paraffin
// (Rubitherm; melting range 311-315 K centred on 40 C), exposed on top to a
// typical January day in Calama (El Loa, -22.50, -68.89, 2326 m) and on the
// bottom to a room at 24 C. Weather is periodic in 24 h of solar time; replace
// it by measured hourly data with -pcm_weather file.csv (columns:
// hour, I [W/m^2], T_a [C], w [m/s]; hour in [0, 24)).
//
// Run from the build directory:
//   mpirun -n 4 ./examples/Heat/NaturalConvectionPCM -pcm_dt 30 -pcm_tend 86400 \
//     -pcm_angle 15 -snes_converged_reason
//
// Options: -pcm_nx, -pcm_ny, -pcm_dt, -pcm_tend, -pcm_angle (deg), -pcm_mushy (T_l - T_s, K;
// PETSc option names are case-insensitive, so not -pcm_dT), -pcm_round (K, 0: Eq. 4),
// -pcm_absorptivity, -pcm_hin, -pcm_tin (C), -pcm_weather, -pcm_output_every,
// -pcm_print_every; Newton via -snes_*, its linear solver via -ksp_*/-pc_*
// (MUMPS), projections under -pcm_proj_ (CG).
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include "OrthogonalProjection.h"

namespace Rodin::Examples::Heat
{
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  /// Kelvin offset.
  constexpr Real T0 = 273.15;

  /**
   * Phase change material: properties of each phase, liquid fraction (Eq. 4),
   * volumetric enthalpy (Eqs. 5-6) and Carman-Kozeny resistance, with the
   * derivatives the Newton Jacobian needs. Temperatures in K.
   */
  struct PCM
  {
    // RT42 (Rubitherm): Mohammed et al., Appl. Sci. 14 (2024) 3212, Table 1;
    // solid density from the Rubitherm data sheet (0.88 kg/l).
    Real rhoS = 880.0, rhoL = 760.0;    ///< kg/m^3
    Real cpS = 2000.0, cpL = 2000.0;    ///< J/(kg K)
    Real kS = 0.2, kL = 0.2;            ///< W/(m K)
    Real L = 165.0e3;                   ///< J/kg
    Real mu = 0.02351;                  ///< Pa s (liquid)
    Real beta = 5.0e-4;                 ///< 1/K
    Real Tm = T0 + 40.0;                ///< K, centre of the mushy interval
    Real dT = 1.0;                      ///< K, T_l - T_s (Biwole et al. 2013)
    Real C = 1.0e5;                     ///< kg/(m^3 s), mushy-zone constant
    Real eps = 1.0e-3;                  ///< Carman-Kozeny regularisation
    Real round = 0.2;                   ///< K, C^1 rounding half-width at T_s, T_l (0: Eq. 4)

    Real Ts() const { return Tm - 0.5 * dT; }
    Real Tl() const { return Tm + 0.5 * dT; }

    /// Eq. (4) with its two kinks rounded: on [T_s - d, T_s + d] the ramp is
    /// replaced by the parabola (T - T_s + d)^2 / (4 d dT), tangent to 0 and to
    /// the ramp, and symmetrically at T_l. Y is C^1, so Newton stays locally
    /// quadratic; a plain kink makes it cycle between the two slopes.
    Real Y(Real T) const
    {
      const Real d = std::min(round, 0.25 * dT);
      if (d <= 0.0)
        return std::clamp((T - Ts()) / dT, Real(0), Real(1));
      if (T <= Ts() - d)
        return 0.0;
      if (T < Ts() + d)
        return (T - Ts() + d) * (T - Ts() + d) / (4.0 * d * dT);
      if (T <= Tl() - d)
        return (T - Ts()) / dT;
      if (T < Tl() + d)
        return 1.0 - (Tl() + d - T) * (Tl() + d - T) / (4.0 * d * dT);
      return 1.0;
    }

    Real dY(Real T) const
    {
      const Real d = std::min(round, 0.25 * dT);
      if (d <= 0.0)
        return (T > Ts() && T < Tl()) ? 1.0 / dT : 0.0;
      if (T <= Ts() - d || T >= Tl() + d)
        return 0.0;
      if (T < Ts() + d)
        return (T - Ts() + d) / (2.0 * d * dT);
      if (T > Tl() - d)
        return (Tl() + d - T) / (2.0 * d * dT);
      return 1.0 / dT;
    }

    /// Lower and upper ends of the phase-change interval (Y in (0, 1)).
    Real Tlow() const { return Ts() - std::min(round, 0.25 * dT); }
    Real Thigh() const { return Tl() + std::min(round, 0.25 * dT); }

    Real rho(Real T) const { return rhoS + (rhoL - rhoS) * Y(T); }
    Real cp(Real T) const { return cpS + (cpL - cpS) * Y(T); }
    Real kappa(Real T) const { return kS + (kL - kS) * Y(T); }
    Real dKappa(Real T) const { return (kL - kS) * dY(T); }

    /// H(T) = int_{T_low}^T rho c_p dT' + rho(T) Y(T) L (Eqs. 5-6). Between
    /// the breakpoints rho c_p is a polynomial of degree <= 4 in T, so 3-point
    /// Gauss-Legendre on each piece is exact.
    Real H(Real T) const
    {
      const Real d = std::min(round, 0.25 * dT);
      const Real b[4] = { Ts() - d, Ts() + d, Tl() - d, Tl() + d };
      const auto a = [this](Real t) { return rho(t) * cp(t); };
      Real h = (T < b[0]) ? rhoS * cpS * (T - b[0]) : 0.0;
      static const Real x[3] = { -0.7745966692414834, 0.0, 0.7745966692414834 };
      static const Real w[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
      for (int k = 0; k < 3; ++k)
      {
        const Real lo = b[k], hi = std::clamp(T, b[k], b[k + 1]);
        const Real c = 0.5 * (lo + hi), r = 0.5 * (hi - lo);
        for (int q = 0; q < 3; ++q)
          h += w[q] * r * a(c + r * x[q]);
      }
      if (T > b[3])
        h += rhoL * cpL * (T - b[3]);
      return h + rho(T) * Y(T) * L;
    }

    /// H'(T) = rho c_p + L (rho_s + 2 (rho_l - rho_s) Y) Y'.
    Real dH(Real T) const
    {
      const Real y = Y(T);
      return rho(T) * cp(T) + L * (rhoS + 2.0 * (rhoL - rhoS) * y) * dY(T);
    }

    /// A(Y) = C (1 - Y)^2 / (Y^3 + eps).
    Real A(Real T) const
    {
      const Real y = Y(T);
      return C * (1.0 - y) * (1.0 - y) / (y * y * y + eps);
    }

    Real dA(Real T) const
    {
      const Real y = Y(T), d = y * y * y + eps;
      return C * dY(T) * (-2.0 * (1.0 - y) * d - 3.0 * y * y * (1.0 - y) * (1.0 - y)) / (d * d);
    }
  };

  /**
   * Hourly weather, periodic in 24 h of solar time, linearly interpolated.
   * Default: typical January day in Calama (El Loa Ad., DMC 220002).
   *  - I: clear-sky global horizontal irradiance for -22.50 deg, 2326 m,
   *    15 January (Laue air-mass model), scaled to 8.3 kWh/m^2/day; peak
   *    1051 W/m^2 at solar noon.
   *  - T_a: asymmetric cosine between T_min = 6 C (06 h) and T_max = 26 C
   *    (15 h), DMC January means for the station.
   *  - w: afternoon wind regime of the Loa basin, 2 m/s at night to 9 m/s
   *    at 15-16 h (approximate).
   */
  class Weather
  {
    public:
      struct Sample
      {
        Real hour, I, Ta, w;  ///< h, W/m^2, C, m/s
      };

      Weather()
      {
        const Real I[24] = { 0, 0, 0, 0, 0, 0, 91, 313, 544, 751, 913, 1016,
                             1051, 1016, 913, 751, 544, 313, 91, 0, 0, 0, 0, 0 };
        const Real Ta[24] = { 12.9, 11.0, 9.3, 7.9, 6.9, 6.2, 6.0, 6.6, 8.3, 11.0, 14.3, 17.7,
                              21.0, 23.7, 25.4, 26.0, 25.8, 25.1, 24.1, 22.7, 21.0, 19.1, 17.0, 15.0 };
        const Real w[24] = { 2.5, 2.3, 2.2, 2.0, 2.0, 2.0, 2.0, 2.2, 2.8, 3.5, 4.5, 5.8,
                             7.0, 8.0, 8.8, 9.2, 9.0, 8.2, 6.8, 5.2, 4.0, 3.3, 2.9, 2.7 };
        for (int h = 0; h < 24; ++h)
          m_samples.push_back({ Real(h), I[h], Ta[h], w[h] });
      }

      explicit Weather(const std::string& csv)
      {
        std::ifstream in(csv);
        if (!in)
          throw std::runtime_error("Cannot open weather file " + csv);
        std::string line;
        while (std::getline(in, line))
        {
          std::replace(line.begin(), line.end(), ',', ' ');
          std::istringstream ss(line);
          Sample s;
          if (ss >> s.hour >> s.I >> s.Ta >> s.w)
            m_samples.push_back(s);
        }
        if (m_samples.size() < 2)
          throw std::runtime_error("Weather file needs at least two rows: " + csv);
        std::sort(m_samples.begin(), m_samples.end(),
                  [](const Sample& a, const Sample& b) { return a.hour < b.hour; });
      }

      /// Interpolated sample at time t (s), t = 0 at solar midnight.
      Sample operator()(Real t) const
      {
        const Real hour = std::fmod(t / 3600.0, 24.0);
        const size_t n = m_samples.size();
        size_t i = n - 1;
        for (size_t k = 0; k < n; ++k)
          if (m_samples[k].hour <= hour)
            i = k;
        const Sample& a = m_samples[i];
        const Sample& b = m_samples[(i + 1) % n];
        Real span = b.hour - a.hour, x = hour - a.hour;
        if (span <= 0.0)
          span += 24.0;
        if (x < 0.0)
          x += 24.0;
        const Real s = x / span;
        return { hour, a.I + s * (b.I - a.I), a.Ta + s * (b.Ta - a.Ta), a.w + s * (b.w - a.w) };
      }

    private:
      std::vector<Sample> m_samples;
  };

  class NaturalConvectionPCM
  {
    public:
      using MeshType = Mesh<Context::MPI>;
      using VectorFES = H1<1, Math::SpatialVector<Real>, MeshType>;
      using ScalarFES = H1<1, Real, MeshType>;
      using VectorGF = PETSc::Variational::GridFunction<VectorFES>;
      using ScalarGF = PETSc::Variational::GridFunction<ScalarFES>;
      using VectorTrial = PETSc::Variational::TrialFunction<VectorGF, VectorFES>;
      using ScalarTrial = PETSc::Variational::TrialFunction<ScalarGF, ScalarFES>;
      using VectorTest = PETSc::Variational::TestFunction<VectorFES>;
      using ScalarTest = PETSc::Variational::TestFunction<ScalarFES>;
      using System = PETSc::Math::LinearSystem;
      using CoupledProblem = Problem<System, VectorTrial, ScalarTrial, ScalarTrial,
                                     VectorTest, ScalarTest, ScalarTest>;
      using CellFES = P0<Real, MeshType>;
      using CellGF = PETSc::Variational::GridFunction<CellFES>;
      using Clock = std::chrono::steady_clock;

      struct Labels
      {
        Attribute bottom = 1, top = 2, sides = 3;
      };

      struct Config
      {
        std::string output = "NaturalConvectionPCM";
        Labels labels;

        Real length = 2.0;     ///< m, roof span (x)
        Real thickness = 0.005; ///< m, PCM layer (y)
        size_t nx = 400;       ///< h = 5 mm
        size_t ny = 20;

        PCM pcm;
        Real gravity = 9.81;   ///< m/s^2
        Real angle = 0.0;      ///< deg, roof slope

        Real absorptivity = 1.0;       ///< solar absorptivity of the top face (Eq. 9: 1)
        Real emissivity = 0.9;         ///< long-wave emissivity of the top face
        Real hIn = 6.0;                ///< W/(m^2 K), 1/R_si for downward flux (ISO 6946)
        Real tIn = T0 + 24.0;          ///< K, room air
        std::string weatherFile;       ///< empty: Calama, January

        Real pressurePenalty = 1.0e-7; ///< 1/(Pa s): eps (p, q) fixes the mean pressure

        // Time step. With h = 5 mm and RT42: cell diffusion time
        // h^2/alpha = 190 s (Fourier alpha dt/h^2 = 0.16 at 30 s); a cell melts
        // in rho L h / q'' ~ 800 s at q'' ~ 800 W/m^2 (~25 steps); the weather
        // varies on the hour; liquid velocities stay below ~1 mm/s
        // (CFL <= 6, implicit + OSGS). dt = 30 s resolves all of these.
        Real dt = 30.0;        ///< s
        Real tEnd = 86400.0;   ///< s, one day
        int outputEvery = 120; ///< one hour at dt = 30 s
        int printEvery = 30;
      };

      static constexpr Real Sigma = 5.67e-8;  ///< W/(m^2 K^4)

      NaturalConvectionPCM(const Context::MPI& context, const Config& cfg)
        : m_cfg(cfg),
          m_weather(cfg.weatherFile.empty() ? Weather() : Weather(cfg.weatherFile)),
          m_mesh(makeMesh(context, m_cfg)),
          m_xdmf(context.getCommunicator(), m_cfg.output),
          m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
          m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
          m_ch(m_mesh),
          m_u(m_vh), m_p(m_sh), m_T(m_sh), m_v(m_vh), m_q(m_sh), m_w(m_sh),
          m_uOld(m_vh), m_TOld(m_sh),
          m_C(m_sh), m_kappa(m_sh), m_dKappa(m_sh),
          m_A(m_sh), m_dA(m_sh), m_Y(m_sh),
          m_ha(m_sh), m_Ta(m_sh), m_Tsky4(m_sh), m_I(m_sh),
          m_piConv(m_vh, "pcm_proj_"), m_piGradP(m_vh, "pcm_proj_"),
          m_piDiv(m_sh, "pcm_proj_"), m_piAdv(m_sh, "pcm_proj_"),
          m_tau1(m_ch), m_tau2(m_ch), m_rhoCpTauT(m_ch),
          m_problem(m_u, m_p, m_T, m_v, m_q, m_w), m_ksp(m_problem), m_snes(m_ksp),
          m_z(m_sh), m_one(m_sh), m_flux(m_z)
      {
        const size_t dim = m_mesh.getSpaceDimension();
        Math::SpatialVector<Real> zero(dim);
        for (size_t i = 0; i < dim; ++i)
          zero(i) = 0.0;
        for (auto* g : { &m_uOld, &m_u.getSolution(), &m_piConv.get(), &m_piGradP.get() })
          *g = zero;
        for (auto* g : { &m_p.getSolution(), &m_piDiv.get(), &m_piAdv.get() })
          *g = Real(0);
        m_one = Real(1);

        // Eq. (8): T = T_a(0), u = 0, fully solid.
        const Real Tinit = T0 + m_weather(0.0).Ta;
        m_TOld = Tinit;
        m_T.getSolution() = Tinit;
        updateCoefficients();
        updateWeather(0.0);

        const auto& L = m_cfg.labels;
        m_topArea = integral(m_one, L.top);
        m_bottomArea = integral(m_one, L.bottom);
        m_volume = integral(m_one);

        setupProblem();

        const size_t nu = m_vh.getSize(), np = m_sh.getSize();
        m_snes.setTolerances(1e-10, 1e-8, 1e-12, 30, 10000)
          .setStateUpdate([this, nu, np](const PETSc::Math::Vector& x) {
            m_u.getSolution().setData(x, 0);
            m_p.getSolution().setData(x, nu);
            m_T.getSolution().setData(x, nu + np);
            updateCoefficients();
          });
        setPhaseChangeTruncation(nu + np, np);
        setLumpedEnthalpy(nu + np);

        m_xdmf.setMesh(m_mesh);
        m_xdmf.add("velocity", m_u.getSolution());
        m_xdmf.add("pressure", m_p.getSolution());
        m_xdmf.add("temperature", m_T.getSolution());
        m_xdmf.add("liquid_fraction", m_Y);

        const auto cells = m_mesh.getCellCount(), vertices = m_mesh.getVertexCount();  // collective
        if (isRoot())
        {
          m_csv.open(m_cfg.output + ".csv");
          m_csv << "t,hour,newtonIterations,maxU,minT_C,maxT_C,topT_C,bottomT_C,"
                   "liquidFraction,qRoom,I,Ta_C,ha\n";
          const auto& m = m_cfg.pcm;
          Alert::Info() << "[mesh] cells=" << cells << " vertices=" << vertices
                        << "  " << m_cfg.length << " m x " << m_cfg.thickness << " m"
                        << "  angle=" << m_cfg.angle << " deg" << Alert::Raise;
          Alert::Info() << "[pcm] Ts=" << m.Ts() - T0 << " C  Tl=" << m.Tl() - T0
                        << " C  L=" << m.L << " J/kg  T(0)=" << Tinit - T0 << " C"
                        << Alert::Raise;
        }
      }

      int run()
      {
        const int steps = static_cast<int>(m_cfg.tEnd / m_cfg.dt + 0.5);
        m_problem.assemble();  // seeds the SNES iterate from the trial solutions
        m_xdmf.write(0.0).flush();
        for (int n = 1; n <= steps; ++n)
        {
          m_t = n * m_cfg.dt;

          const auto t0 = Clock::now();
          updateWeather(m_t);     // boundary data at t^{n+1}
          updateStabilization();  // tau and Pi at (u^n, p^n, T^n)
          m_jacobianState = -1;   // Rodin reassembles at the start of a solve
          m_snes.solve();
          const Real elapsed = std::chrono::duration<Real>(Clock::now() - t0).count();
          if (!m_snes.converged())
            throw std::runtime_error("Newton did not converge at t = " + std::to_string(m_t));

          updateCoefficients();  // at the converged T^{n+1}
          m_uOld.setData(m_u.getSolution().getData());
          m_TOld.setData(m_T.getSolution().getData());

          report(n, steps, elapsed);
          if (n % m_cfg.outputEvery == 0 || n == steps)
            m_xdmf.write(m_t).flush();
        }
        m_xdmf.close();
        return 0;
      }

    private:
      /// Structured triangulation of [0, length] x [0, thickness] built on
      /// rank 0, boundary facets tagged by position, then partitioned.
      static MeshType makeMesh(const Context::MPI& context, const Config& cfg)
      {
        const auto& comm = context.getCommunicator();
        Rodin::MPI::Sharder sharder(context);
        if (comm.rank() == 0)
        {
          auto mesh = Mesh<Context::Local>::UniformGrid(
              Polytope::Type::Triangle, { cfg.nx + 1, cfg.ny + 1 });
          const Real hx = cfg.length / static_cast<Real>(cfg.nx);
          const Real hy = cfg.thickness / static_cast<Real>(cfg.ny);
          for (auto it = mesh.getVertex(); it; ++it)
          {
            Math::SpatialPoint x = mesh.getVertexCoordinates(it->getIndex());
            x(0) *= hx;
            x(1) *= hy;
            mesh.setVertexCoordinates(it->getIndex(), x);
          }
          mesh.flush();
          connect(mesh);

          const Real tol = 1e-3 * std::min(hx, hy);
          const auto& L = cfg.labels;
          std::vector<std::pair<Index, Attribute>> tags;
          for (auto it = mesh.getBoundary(); it; ++it)
          {
            Math::SpatialPoint c(2);
            c.setZero();
            for (const auto& v : it->getVertices())
              c += mesh.getVertexCoordinates(v);
            c /= static_cast<Real>(it->getVertices().size());
            Attribute a = L.sides;
            if (c(1) < tol)
              a = L.bottom;
            else if (c(1) > cfg.thickness - tol)
              a = L.top;
            tags.emplace_back(it->getIndex(), a);
          }
          for (const auto& [f, a] : tags)
            mesh.setAttribute({ 1, f }, a);

#ifdef RODIN_USE_SCOTCH
          Scotch::Partitioner partitioner(mesh);
#else
          BalancedCompactPartitioner partitioner(mesh);
#endif
          partitioner.partition(static_cast<size_t>(comm.size()));
          sharder.shard(partitioner);
          sharder.scatter(0);
        }
        MeshType mesh = sharder.gather(0);
        connect(mesh);
        mesh.reconcile(1);
        return mesh;
      }

      template <class M>
      static void connect(M& mesh)
      {
        const size_t D = mesh.getDimension();
        mesh.getConnectivity().compute(D, D);
        mesh.getConnectivity().compute(D, 0);
        mesh.getConnectivity().compute(D, D - 1);
        mesh.getConnectivity().compute(D - 1, D);
        mesh.getConnectivity().compute(D - 1, 0);
        mesh.getConnectivity().compute(1, 0);
      }

      static Real cellSize(const Point& p)
      {
        return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
      }

      /// The enthalpy rate with lumped mass, m_i (H(T_i) - H(T_i^n))/dt, and
      /// its exact derivative m_i H'(T_i)/dt on the diagonal, m_i = (1, phi_i).
      /// Rodin's SNES residual/Jacobian are wrapped: the original assembles the
      /// form, the wrapper adds these nodal terms on the T rows (offset + i).
      void setLumpedEnthalpy(size_t offset)
      {
        m_flux = Integral(m_one, m_z);
        m_flux.assemble();
        PetscErrorCode ierr;
        PetscInt rb = 0, re = 0;
        ierr = VecGetOwnershipRange(m_flux.getVector(), &rb, &re);
        assert(ierr == PETSC_SUCCESS);
        const PetscScalar* m;
        ierr = VecGetArrayRead(m_flux.getVector(), &m);
        assert(ierr == PETSC_SUCCESS);
        m_lumped.assign(m, m + (re - rb));
        ierr = VecRestoreArrayRead(m_flux.getVector(), &m);
        assert(ierr == PETSC_SUCCESS);
        m_rows.resize(re - rb);
        for (PetscInt k = 0; k < re - rb; ++k)
          m_rows[k] = static_cast<PetscInt>(offset) + rb + k;

        ::SNES snes = m_snes.getHandle();
        ierr = SNESGetFunction(snes, PETSC_NULLPTR, &m_rodinResidual, &m_rodinResidualCtx);
        assert(ierr == PETSC_SUCCESS);
        ::Mat A = PETSC_NULLPTR, P = PETSC_NULLPTR;
        ierr = SNESGetJacobian(snes, &A, &P, &m_rodinJacobian, &m_rodinJacobianCtx);
        assert(ierr == PETSC_SUCCESS);
        ierr = SNESSetFunction(snes, PETSC_NULLPTR, &NaturalConvectionPCM::residual, this);
        assert(ierr == PETSC_SUCCESS);
        ierr = SNESSetJacobian(snes, A, P, &NaturalConvectionPCM::jacobian, this);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      /// Owned T values (state or t^n), FES-local numbering.
      template <class F>
      PetscErrorCode lumped(const ScalarGF& T, const ScalarGF& TOld, F&& f) const
      {
        const PetscScalar *t, *t0;
        PetscCall(VecGetArrayRead(T.getData(), &t));
        PetscCall(VecGetArrayRead(TOld.getData(), &t0));
        for (size_t k = 0; k < m_lumped.size(); ++k)
          f(k, t[k], t0[k]);
        PetscCall(VecRestoreArrayRead(TOld.getData(), &t0));
        PetscCall(VecRestoreArrayRead(T.getData(), &t));
        return PETSC_SUCCESS;
      }

      static PetscErrorCode residual(::SNES snes, ::Vec x, ::Vec f, void* ctx)
      {
        auto& self = *static_cast<NaturalConvectionPCM*>(ctx);
        PetscCall(self.m_rodinResidual(snes, x, f, self.m_rodinResidualCtx));  // syncs the state
        const PCM& m = self.m_cfg.pcm;
        const Real dt = self.m_cfg.dt;
        std::vector<PetscScalar> values(self.m_lumped.size());
        PetscCall(self.lumped(self.m_T.getSolution(), self.m_TOld, [&](size_t k, Real t, Real t0) {
          values[k] = self.m_lumped[k] * (m.H(t) - m.H(t0)) / dt; }));
        PetscCall(VecSetValues(f, static_cast<PetscInt>(values.size()), self.m_rows.data(),
                               values.data(), ADD_VALUES));
        PetscCall(VecAssemblyBegin(f));
        PetscCall(VecAssemblyEnd(f));
        return PETSC_SUCCESS;
      }

      static PetscErrorCode jacobian(::SNES snes, ::Vec x, ::Mat J, ::Mat P, void* ctx)
      {
        auto& self = *static_cast<NaturalConvectionPCM*>(ctx);
        PetscCall(self.m_rodinJacobian(snes, x, J, P, self.m_rodinJacobianCtx));
        // With PETSc >= 3.20 Rodin reuses the operator assembled at the same
        // iterate (state cache); the diagonal is then already in it.
#if PETSC_VERSION_GE(3, 20, 0)
        PetscObjectState state;
        PetscCall(VecGetState(x, &state));
        if (state == self.m_jacobianState)
          return PETSC_SUCCESS;
        self.m_jacobianState = state;
#endif
        const PCM& m = self.m_cfg.pcm;
        const Real dt = self.m_cfg.dt;
        PetscCall(self.lumped(self.m_T.getSolution(), self.m_TOld, [&](size_t k, Real t, Real) {
          const PetscInt row = self.m_rows[k];
          const PetscScalar v = self.m_lumped[k] * m.dH(t) / dt;
          (void)MatSetValues(P, 1, &row, 1, &row, &v, ADD_VALUES); }));
        PetscCall(MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY));
        if (J != P)
        {
          PetscCall(MatCopy(P, J, SAME_NONZERO_PATTERN));
        }
        return PETSC_SUCCESS;
      }

      /// Line-search post-check on the T block [offset, offset + size) of the
      /// monolithic iterate (u | p | T): a node whose update crosses the whole
      /// mushy interval is moved to T_m.
      void setPhaseChangeTruncation(size_t offset, size_t size)
      {
        m_truncation = { static_cast<PetscInt>(offset), static_cast<PetscInt>(offset + size),
                         m_cfg.pcm.Tlow(), m_cfg.pcm.Thigh(), m_cfg.pcm.Tm };
        ::SNES snes = m_snes.getHandle();
        PetscErrorCode ierr = SNESSetType(snes, SNESNEWTONLS);
        assert(ierr == PETSC_SUCCESS);
        ::SNESLineSearch ls;
        ierr = SNESGetLineSearch(snes, &ls);
        assert(ierr == PETSC_SUCCESS);
        ierr = SNESLineSearchSetPostCheck(ls, &NaturalConvectionPCM::truncate, &m_truncation);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      struct Truncation
      {
        PetscInt begin = 0, end = 0;
        Real Ts = 0.0, Tl = 0.0, Tm = 0.0;
      };

      static PetscErrorCode truncate(::SNESLineSearch, ::Vec X, ::Vec, ::Vec W,
                                     PetscBool*, PetscBool* changedW, void* ctx)
      {
        const auto& c = *static_cast<const Truncation*>(ctx);
        PetscInt lo = 0, hi = 0;
        PetscCall(VecGetOwnershipRange(W, &lo, &hi));
        const PetscScalar* x;
        PetscScalar* w;
        PetscCall(VecGetArrayRead(X, &x));
        PetscCall(VecGetArray(W, &w));
        PetscInt changed = 0;
        for (PetscInt i = std::max(lo, c.begin); i < std::min(hi, c.end); ++i)
        {
          const PetscScalar xi = x[i - lo], wi = w[i - lo];
          if ((xi <= c.Ts && wi >= c.Tl) || (xi >= c.Tl && wi <= c.Ts))  // [T_low, T_high]
          {
            w[i - lo] = c.Tm;
            ++changed;
          }
        }
        PetscCall(VecRestoreArray(W, &w));
        PetscCall(VecRestoreArrayRead(X, &x));
        *changedW = changed > 0 ? PETSC_TRUE : PETSC_FALSE;
        return PETSC_SUCCESS;
      }

      /// P1 interpolants of the T-dependent coefficients at the current iterate.
      void updateCoefficients()
      {
        const auto& T = m_T.getSolution();
        const PCM& m = m_cfg.pcm;
        const auto at = [&T](auto f) {
          return RealFunction([&T, f](const Point& p) { return f(T.getValue(p)); });
        };
        m_C.project(at([&m](Real t) { return m.dH(t); }));
        m_kappa.project(at([&m](Real t) { return m.kappa(t); }));
        m_dKappa.project(at([&m](Real t) { return m.dKappa(t); }));
        m_A.project(at([&m](Real t) { return m.A(t); }));
        m_dA.project(at([&m](Real t) { return m.dA(t); }));
        m_Y.project(at([&m](Real t) { return m.Y(t); }));
      }

      /// Top-face data at time t: h_a = 8.91 + 2 w, T_sky = 0.0552 T_a^1.5.
      void updateWeather(Real t)
      {
        m_sample = m_weather(t);
        const Real Ta = T0 + m_sample.Ta;
        const Real Tsky = 0.0552 * std::pow(Ta, 1.5);
        m_haValue = 8.91 + 2.0 * m_sample.w;
        m_ha = m_haValue;
        m_Ta = Ta;
        m_Tsky4 = Tsky * Tsky * Tsky * Tsky;
        m_I = m_sample.I;
      }

      /// Cellwise tau (P0) and the projections, both at t^n.
      void updateStabilization()
      {
        const PCM& m = m_cfg.pcm;
        const Real rho = m.rhoL, nu = m.mu / m.rhoL;
        const auto tau1 = [this, &m, rho, nu](const Point& p) {
          const auto u = m_uOld.getValue(p);
          const Real h = cellSize(p);
          return 1.0 / (4.0 * nu / (h * h) + 2.0 * std::sqrt(Math::dot(u, u)) / h
                        + m.A(m_TOld.getValue(p)) / rho);
        };
        m_tau1.project(RealFunction(tau1));
        m_tau2.project(RealFunction([rho, tau1](const Point& p) {
          const Real h = cellSize(p);
          return rho * h * h / (4.0 * tau1(p)); }));
        m_rhoCpTauT.project(RealFunction([this, &m](const Point& p) {
          const auto u = m_uOld.getValue(p);
          const Real T = m_TOld.getValue(p), h = cellSize(p);
          const Real rhoCp = m.rho(T) * m.cp(T), alpha = m.kappa(T) / rhoCp;
          return rhoCp / (4.0 * alpha / (h * h) + 2.0 * std::sqrt(Math::dot(u, u)) / h); }));

        const auto& U = m_u.getSolution();
        m_piConv.project(Mult(Jacobian(U), U));
        m_piGradP.project(Grad(m_p.getSolution()));
        m_piDiv.project(Div(U));
        m_piAdv.project(Dot(U, Grad(m_T.getSolution())));
      }

      void setupProblem()
      {
        const PCM& m = m_cfg.pcm;
        const Real rho = m.rhoL, dt = m_cfg.dt, mu = m.mu;
        const Real rhoBeta = m.rhoL * m.beta, Tref = m.Tl();
        const Real sigmaEps = Sigma * m_cfg.emissivity, alphaS = m_cfg.absorptivity;
        const Real hIn = m_cfg.hIn, tIn = m_cfg.tIn, epsP = m_cfg.pressurePenalty;
        const auto& L = m_cfg.labels;

        // Gravity in the roof frame: x along the slope, y across the layer.
        const Real theta = m_cfg.angle * M_PI / 180.0;
        const VectorFunction g{ -m_cfg.gravity * std::sin(theta), -m_cfg.gravity * std::cos(theta) };

        // State (current Newton iterate) and corrections (trial functions).
        const auto& U = m_u.getSolution();
        const auto& P = m_p.getSolution();
        const auto& Th = m_T.getSolution();
        const auto& du = m_u;
        const auto& dp = m_p;
        const auto& dT = m_T;

        const auto convU = Mult(Jacobian(U), U);                                // (grad U) U
        const auto dConvU = Mult(Jacobian(du), U) + Mult(Jacobian(U), du);      // its derivative
        const auto streamV = Mult(Jacobian(m_v), U);                            // (grad v) U
        const auto streamW = Dot(U, Grad(m_w));                                 // U.grad w
        const auto advT = Dot(U, Grad(Th));                                     // U.grad T
        const auto dAdvT = Dot(U, Grad(dT));                                    // its derivative,
        const auto dAdvU = Dot(du, Grad(Th));                                   // split by field

        m_problem =
            // ---- Jacobian: momentum --------------------------------------
              Integral(Dot((rho / dt) * du + 0.5 * rho * (Div(du) * U) + 0.5 * rho * (Div(U) * du), m_v))
            + Integral(dConvU, rho * m_v + rho * (m_tau1 * streamV))
            + rho * Integral(du, Mult(Transpose(Jacobian(m_v)), m_tau1 * (convU - m_piConv.get())))
            + mu * Integral(Jacobian(du), Jacobian(m_v))
            + Integral(m_tau2 * Div(du), Div(m_v))
            - Integral(dp, Div(m_v))
            + Integral(m_A * du, m_v)                                  // Carman-Kozeny
            + Integral(m_dA * dT, Dot(U, m_v))                         //   A'(T) dT u
            + Integral(rhoBeta * dT, Dot(g, m_v))                      // buoyancy

            // ---- Jacobian: continuity ------------------------------------
            + Integral(Div(du), m_q)
            + (1.0 / rho) * Integral(m_tau1 * Grad(dp), Grad(m_q))
            + epsP * Integral(dp, m_q)

            // ---- Jacobian: energy ----------------------------------------
            + Integral(m_C * dAdvT, m_w)  // M_L H'(T)/dt: added in jacobian()
            + Integral(m_rhoCpTauT * dAdvT, streamW)
            + Integral(m_kappa * Grad(dT), Grad(m_w))
            + Integral(m_dKappa * dT, Dot(Grad(Th), Grad(m_w)))
            + Integral(dAdvU, m_C * m_w + m_rhoCpTauT * streamW)
            + Integral(m_rhoCpTauT * (advT - m_piAdv.get()) * du, Grad(m_w))
            + BoundaryIntegral((m_ha + 4.0 * sigmaEps * Pow(Th, 3)) * dT, m_w).over(L.top)
            + hIn * BoundaryIntegral(dT, m_w).over(L.bottom)

            // ---- Residual: momentum --------------------------------------
            + Integral(Dot((rho / dt) * (U - m_uOld) + rho * convU + 0.5 * rho * (Div(U) * U), m_v))
            + rho * Integral(m_tau1 * (convU - m_piConv.get()), streamV)
            + mu * Integral(Jacobian(U), Jacobian(m_v))
            + Integral(m_tau2 * (Div(U) - m_piDiv.get()) - P, Div(m_v))
            + Integral(Dot(m_A * U + (rhoBeta * (Th - Tref)) * g, m_v))

            // ---- Residual: continuity ------------------------------------
            + Integral(Div(U), m_q)
            + (1.0 / rho) * Integral(m_tau1 * (Grad(P) - m_piGradP.get()), Grad(m_q))
            + epsP * Integral(P, m_q)

            // ---- Residual: energy ----------------------------------------
            + Integral(m_C * advT, m_w)   // dH/dt: lumped, added in residual()
            + Integral(m_rhoCpTauT * (advT - m_piAdv.get()), streamW)
            + Integral(m_kappa * Grad(Th), Grad(m_w))
            + BoundaryIntegral((m_ha * (Th - m_Ta) + sigmaEps * (Pow(Th, 4) - m_Tsky4)
                                - alphaS * m_I) * m_w).over(L.top)
            + hIn * BoundaryIntegral((Th - tIn) * m_w).over(L.bottom)

            // ---- No-slip on every wall: correction = 0 - state -----------
            + DirichletBC(du, -U).on(FlatSet<Attribute>{ L.bottom, L.top, L.sides });
      }

      template <class Expression>
      Real integral(const Expression& f, Attribute a)
      {
        m_flux = BoundaryIntegral(f, m_z).over(a);
        m_flux.assemble();
        return m_flux(m_one);
      }

      template <class Expression>
      Real integral(const Expression& f)
      {
        m_flux = Integral(f, m_z);
        m_flux.assemble();
        return m_flux(m_one);
      }

      /// Surface means of T, mean liquid fraction and heat flux into the room
      /// h_in (T - T_in) averaged over the bottom face (W/m^2, > 0 heats it).
      void report(int step, int steps, Real elapsed)
      {
        const auto& L = m_cfg.labels;
        const auto& u = m_u.getSolution();
        const auto& T = m_T.getSolution();

        const Real topT = integral(T, L.top) / m_topArea - T0;
        const Real bottomT = integral(T, L.bottom) / m_bottomArea - T0;
        const Real qRoom = m_cfg.hIn * (bottomT + T0 - m_cfg.tIn);
        const Real Ybar = integral(m_Y) / m_volume;
        const Real umax = std::max(std::abs(u.max()), std::abs(u.min()));
        const Real Tmax = T.max() - T0, Tmin = T.min() - T0;
        const int its = static_cast<int>(m_snes.getIterationNumber());

        if (!isRoot())
          return;
        m_csv << m_t << ',' << m_sample.hour << ',' << its << ',' << umax << ',' << Tmin << ','
              << Tmax << ',' << topT << ',' << bottomT << ',' << Ybar << ',' << qRoom << ','
              << m_sample.I << ',' << m_sample.Ta << ',' << m_haValue << '\n';
        m_csv.flush();
        if (step % m_cfg.printEvery == 0 || step == steps)
        {
          Alert::Info() << "step " << step << "/" << steps << "  t=" << m_t / 3600.0
                        << " h  newton=" << its << "  max|u|=" << umax << " m/s  T in ["
                        << Tmin << ", " << Tmax << "] C  top=" << topT << " C  bottom="
                        << bottomT << " C  Y=" << Ybar << "  q_room=" << qRoom
                        << " W/m^2  I=" << m_sample.I << "  [" << elapsed << " s]"
                        << Alert::Raise;
          std::cout.flush();
        }
      }

      bool isRoot() const { return m_mesh.getContext().getCommunicator().rank() == 0; }

      Config m_cfg;
      Weather m_weather;
      MeshType m_mesh;
      IO::XDMF m_xdmf;
      VectorFES m_vh;
      ScalarFES m_sh;
      CellFES m_ch;

      // Trial functions are Newton corrections; their solutions hold the state.
      VectorTrial m_u;
      ScalarTrial m_p;
      ScalarTrial m_T;
      VectorTest m_v;
      ScalarTest m_q;
      ScalarTest m_w;
      VectorGF m_uOld;
      ScalarGF m_TOld;

      // T-dependent coefficients (P1 interpolants of the iterate).
      ScalarGF m_C, m_kappa, m_dKappa, m_A, m_dA, m_Y;

      // Top-face weather data at t^{n+1} (constant fields).
      ScalarGF m_ha, m_Ta, m_Tsky4, m_I;

      OrthogonalProjection<VectorFES> m_piConv, m_piGradP;
      OrthogonalProjection<ScalarFES> m_piDiv, m_piAdv;

      CellGF m_tau1, m_tau2, m_rhoCpTauT;

      CoupledProblem m_problem;
      Solver::KSP m_ksp;
      Solver::SNES m_snes;

      ScalarTest m_z;
      ScalarGF m_one;
      LinearForm<ScalarFES, ::Vec> m_flux;

      Truncation m_truncation;

      // Lumped enthalpy rate: (1, phi_i), monolithic T rows, wrapped callbacks.
      std::vector<PetscScalar> m_lumped;
      std::vector<PetscInt> m_rows;
      PetscErrorCode (*m_rodinResidual)(::SNES, ::Vec, ::Vec, void*) = nullptr;
      void* m_rodinResidualCtx = nullptr;
      PetscErrorCode (*m_rodinJacobian)(::SNES, ::Vec, ::Mat, ::Mat, void*) = nullptr;
      void* m_rodinJacobianCtx = nullptr;
      PetscObjectState m_jacobianState = -1;  ///< x state of the last diagonal update

      Real m_t = 0.0;
      Weather::Sample m_sample{};
      Real m_haValue = 0.0;
      Real m_topArea = 0.0, m_bottomArea = 0.0, m_volume = 0.0;
      std::ofstream m_csv;
  };
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const auto setDefault = [](const char* name, const char* value) {
    PetscBool set = PETSC_FALSE;
    PetscErrorCode ierr = PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name, &set);
    if (ierr == PETSC_SUCCESS && !set)
      ierr = PetscOptionsSetValue(PETSC_NULLPTR, name, value);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  };
  // Newton step: direct (MUMPS). Projections: CG + Jacobi on the mass matrix.
  setDefault("-ksp_type", "preonly");
  setDefault("-pc_type", "lu");
  setDefault("-pc_factor_mat_solver_type", "mumps");
  setDefault("-mat_mumps_icntl_20", "0");
  setDefault("-mat_mumps_icntl_21", "0");
  setDefault("-snes_linesearch_type", "basic");  // the truncated iterate is the step
  setDefault("-pcm_proj_ksp_type", "cg");
  setDefault("-pcm_proj_pc_type", "jacobi");
  setDefault("-pcm_proj_ksp_rtol", "1e-10");

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  int status = 0;
  try
  {
    using Rodin::Examples::Heat::NaturalConvectionPCM;
    using Rodin::Examples::Heat::T0;
    NaturalConvectionPCM::Config cfg;

    char buffer[512];
    PetscBool got = PETSC_FALSE;
    PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_weather", buffer, sizeof(buffer), &got);
    if (got)
      cfg.weatherFile = buffer;
    PetscInt nx = static_cast<PetscInt>(cfg.nx), ny = static_cast<PetscInt>(cfg.ny);
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_nx", &nx, PETSC_NULLPTR);
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_ny", &ny, PETSC_NULLPTR);
    cfg.nx = static_cast<size_t>(nx);
    cfg.ny = static_cast<size_t>(ny);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_dt", &cfg.dt, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_tend", &cfg.tEnd, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_angle", &cfg.angle, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_mushy", &cfg.pcm.dT, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_round", &cfg.pcm.round, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_absorptivity", &cfg.absorptivity, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_hin", &cfg.hIn, PETSC_NULLPTR);
    PetscReal tIn = cfg.tIn - T0;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_tin", &tIn, PETSC_NULLPTR);
    cfg.tIn = T0 + tIn;
    PetscInt every = cfg.outputEvery, print = cfg.printEvery;
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_output_every", &every, PETSC_NULLPTR);
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-pcm_print_every", &print, PETSC_NULLPTR);
    cfg.outputEvery = static_cast<int>(every);
    cfg.printEvery = static_cast<int>(print);

    NaturalConvectionPCM simulation(context, cfg);
    status = simulation.run();
  }
  catch (const std::exception& e)
  {
    std::cerr << "Fatal error: " << e.what() << "\n";
    status = 1;
  }
  PetscFinalize();
  return status;
}
