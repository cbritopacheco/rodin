/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ArterialLesion3D_viscous_logarithmic_implicit.h
 * @brief Pulsatile sPTT/Oldroyd-B flow through an idealised stenosis or
 *        fusiform aneurysm in a 3D pipe; the formulation of
 *        ArterialLesion2D_viscous_logarithmic_implicit, term by term.
 *
 * Formulation. T = 2 eta_s eps(u) + sigma, sigma = (eta_p/lambda0)(exp(psi)
 * - I), (u, p, psi) solved by Newton, P1 spaces, BDF1 or BDF2
 * (-al_time_order; variable-step, with step rejection and control,
 * Config::adaptiveTimeStep). Momentum and continuity carry the term-by-term
 * (split) OSS with dynamic subscales: u1 for convection, u2 for the pressure
 * gradient, u3 for div sigma, and the algebraic pressure subscale (grad-div).
 * Three discretisations of the constitutive equation (Config::formulation):
 *
 *  - "psi_rb" (default): the psi equation (Fattal and Kupferman 2004)
 *      R_psi = d_t psi + (u.grad) psi + N(psi, L) = 0,
 *      N = -Dexp^{-1}[L E + E L^T - c eps(u)] + (f/lambda)(I - exp(-psi)),
 *    BDF on psi, with the OSS of its full residual: a log-conformation
 *    subscale psi~, (d_t + 1/tau) psi~ = -P^perp[R_psi - d_t psi_h] (the BDF
 *    derivative of psi_h is in the FE space; Config::subscaleTimeResidual), tau = (lambda/(2
 *    eta_p)) alpha3, entering as <psi~, -(u^n.grad) chi> (optionally + JN^*[chi],
 *    JN = dN/dpsi; Config::adjointReaction) and, through sigma~ =
 *    (eta_p/lambda0) Dexp[psi~], as <sigma~, eps(v)> in the momentum (the u-psi
 *    compatibility). Both are consistent and vanish where R_psi is resolved.
 *    The psi equation keeps every eigenvalue of psi at the scale O(1): in the exp form below the
 *    equation for an eigenvalue l of psi is scaled by e^l, so at the
 *    Weissenberg numbers of blood (lambda gd ~ 10^3) the small eigenvalues
 *    are not determined (psi_min ran to -20 and Newton stalled in the S75
 *    test at t = 0.05 s; psi_rb kept psi_min > -6.6 and ran on to
 *    t = 0.115 s without a rejected step).
 *  - "exp" (Moreno, Codina, Baiges, Castillo, CMAME 354, 2019;
 *    Moreno, Codina, Baiges, CMAME 367, 2020). With E = exp(psi), f the PTT
 *    factor and L = grad u (L_ij = d u_i/d x_j),
 *
 *      R = f (E - I)/(2 lambda0) - eps(u)
 *          + (lambda/(2 lambda0)) [d_t E + (u.grad) E - (L E + E L^T) + 2 eps(u)] = 0,
 *
 *    the sigma equation over 2 eta_p, tested with chi; BDF on E. The
 *    constitutive equation
 *    carries the OSS of its full residual, sum_K <s~, A(chi) + eps(v)>,
 *      A(chi) = f chi/(2 eta_p) - (lambda/(2 eta_p))((u.grad) chi + chi L + L^T chi),
 *      (lambda/(2 eta_p)) d_t s~ + s~/alpha3 = -P^perp[R],
 *    because its terms are large where they balance (stress layers) while R
 *    is small: projecting each term separately fails there (Moreno et al.
 *    2020, Table 1) and adds dissipation where the equation is satisfied.
 *    The subscales live at the quadrature points; P^perp[R] = R - Pi[R^n].
 *  - "psi": the former scheme, the psi equation with the split OSS of the
 *    former driver (S1 dynamic convection, pressure and div sigma; S2
 *    grad-div; S3 u-psi compatibility Pi[eps(u)] and psi advection). S3 acts
 *    on eps(u) - Pi[eps(u)] with a weight ~eta_p/2 whether or not the
 *    equation is resolved: an artificial viscosity of several eta_s for blood.
 *
 * The only changes from the planar driver are those of the dimension:
 *
 *  - psi is a symmetric 3x3 tensor, six unknowns per node. Rodin's vector H1
 *    space holds at most three components (Math::SpatialVector), so psi is
 *    carried by two P1 vector fields, the diagonal block psi_D = (xx, yy, zz)
 *    and the off-diagonal block psi_O = (xy, xz, yz):
 *
 *      psi = sum_k psi_D,k E_k + sum_k psi_O,k E_{3+k},
 *      E_0..2 = e_i e_i,  E_3 = e_x e_y + e_y e_x,  E_4 = e_x e_z + e_z e_x,
 *      E_5 = e_y e_z + e_z e_y,
 *
 *    the 2D Voigt basis extended to 3D. Every pairing with the test tensor
 *    chi = (chi_D, chi_O) is therefore split into its D and O blocks; the
 *    mass and advection operators are block-diagonal (E_i : E_{3+k} = 0),
 *    the Newton Jacobian JN and the velocity map are not.
 *  - The PTT factor is f = 1 + epsilon (lambda/lambda0)(tr exp(psi) - 3).
 *  - Daleckii-Krein on 3x3 eigenproblems.
 *
 * Domain. A pipe of diameter D whose radius follows
 *
 *   r(x)/R = 1 + a f(x),  f(x) = (1 + cos(pi x/ell))/2 for |x| <= ell,
 *
 * with a = sqrt(1 - S) - 1 (stenosis, S the AREA reduction at the throat) or
 * a = Gam - 1 (aneurysm), axis along x. The meshes are written in units of D
 * by make_lesion_mesh_3d.py and scaled by D here.
 *
 * Inlet. u = U(r, t) e_x, r = sqrt(y^2 + z^2), the fully developed pulsatile
 * pipe flow of the same fluid at the imposed flow rate Q(t) = Ubar (pi D^2/4)
 * q(t). Mode 0 is Poiseuille, 2 (1 - r^2/R^2); mode n >= 1 is the Womersley
 * pipe profile with the complex viscosity eta*(n omega) = eta_s + eta_p/(1 +
 * i n omega lambda). As in 2D, the shear stress of unidirectional flow obeys
 * the linear Maxwell equation, so the profile is exact for Oldroyd-B and the
 * linear-viscoelastic limit for sPTT. psi = 0 is imposed on the inlet.
 *
 * Outlet. Traction p_out n plus the directional do-nothing (backflow) term.
 *
 * Fluids (Config::fluid), as in the axisymmetric driver:
 *  - "sptt"  constant solvent viscosity eta_s (Newtonian: Wi -> 0);
 *  - "gsptt" one sPTT mode plus the shear-thinning solvent eta_s(gd) =
 *            eta_CY(gd) - eta_p/f(gd), whose steady shear viscosity is exactly
 *            the Carreau-Yasuda eta_CY; eta_s at u^n (u* with BDF2);
 *  - "cy"    Carreau-Yasuda: here the "gsptt" system with a vanishing polymer
 *            (eta_p = 1e-6 eta_inf, De = 1e-5, epsilon 0.5: passive, psi ~ 0),
 *            so eta_s = eta_CY to 1e-6; the
 *            psi unknowns are kept (same cost as "gsptt").
 * For "cy" and "gsptt" the inlet mean mode is the developed steady
 * Carreau-Yasuda profile and Re, Wo use eta_ref (default eta_inf).
 *
 * Dimensionless groups, with eta_0 = eta_s + eta_p:
 *   Re = rho Ubar D/eta_0,  Wo = (D/2) sqrt(omega rho/eta_0),
 *   Wi = lambda Ubar/D,     De = lambda/T.
 */
#ifndef EXAMPLES_VISCOELASTICFLUIDS_ARTERIALLESION3D_VISCOUS_LOGARITHMIC_IMPLICIT_H
#define EXAMPLES_VISCOELASTICFLUIDS_ARTERIALLESION3D_VISCOUS_LOGARITHMIC_IMPLICIT_H

#include <limits>
#include <array>
#include <complex>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#include "OrthogonalProjection.h"

namespace Rodin::Examples::ViscoelasticFluids
{
  /**
   * @brief Fully developed pulsatile flow of a Maxwell/Oldroyd-B fluid in a
   *        pipe of radius R, at a prescribed flow-rate waveform.
   *
   * @details q(t) = sum_n c_n e^{i n omega t}, c_0 = 1, c_{-n} = conj(c_n).
   *          Each mode has the profile of unit cross-sectional mean
   *
   *            phi_0(r) = 2 (1 - r^2/R^2),
   *            phi_n(r) = [1 - I_0(k_n r)/I_0(k_n R)]
   *                       / [1 - 2 I_1(k_n R)/(k_n R I_0(k_n R))],
   *
   *          k_n^2 = i n omega rho / eta*(n omega), the solution of
   *          eta* (phi'' + phi'/r) - i n omega rho phi = const, phi(R) = 0.
   *          U(r, t)/Ubar = sum_n c_n phi_n(r) e^{i n omega t}.
   */
  class PulsatilePipeInflow
  {
    public:
      using Real = Rodin::Real;
      using Complex = std::complex<Real>;

      PulsatilePipeInflow() = default;

      /// @brief q(t) = 1 + A sin(omega t).
      PulsatilePipeInflow& setSinusoidal(Real amplitude);

      /// @brief q(t) from a two-column file (t, Q), normalised by its mean;
      ///        only the shape is kept, the time column is rescaled to the
      ///        period of the run.
      PulsatilePipeInflow& load(const std::string& path, size_t harmonics);

      /// @brief Fixes the fluid and the pipe; computes every k_n.
      PulsatilePipeInflow& setFluid(Real rho, Real etaS, Real etaP, Real lambda,
        Real omega, Real radius);

      /// @brief q(t).
      Real flowRate(Real t) const;

      /// @brief U(r, t)/Ubar.
      Real velocity(Real r, Real t) const;

      /// @brief Replaces the mean mode phi_0 = 2 (1 - r^2/R^2) by a tabulated
      ///        profile of unit cross-sectional mean (r ascending, from 0 to
      ///        R): the developed steady profile of a shear-thinning fluid.
      PulsatilePipeInflow& setMeanProfile(std::vector<Real> r, std::vector<Real> phi);

      /// @brief max_t q(t), sampled.
      Real getPeak() const;

      size_t getHarmonicCount() const noexcept
      {
        return m_c.size() - 1;
      }

      /// @brief e^{-z} I_nu(z), nu = 0, 1, Re z >= 0.
      static Complex besselIScaled(int nu, Complex z);

    private:
      /// @brief phi_n(r), n >= 1.
      Complex mode(size_t n, Real r) const;

      std::vector<Complex> m_c{ Complex(1.0) };
      std::vector<Complex> m_k;
      std::vector<Real> m_meanR;
      std::vector<Real> m_meanPhi;
      Real m_omega = 0.0;
      Real m_radius = 1.0;
  };

  class ArterialLesion3DViscousLogImplicit
  {
    public:
      using Real = Rodin::Real;
      using Attribute = Rodin::Geometry::Attribute;
      using MeshType = Rodin::Geometry::Mesh<Rodin::Context::MPI>;
      using AttributeSet = Rodin::FlatSet<Attribute>;

      using VectorFESType =
        Rodin::Variational::H1<1, Rodin::Math::SpatialVector<Real>, MeshType>;
      using ScalarFESType = Rodin::Variational::H1<1, Real, MeshType>;

      using VectorGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<VectorFESType>;
      using ScalarGridFunctionType =
        Rodin::PETSc::Variational::GridFunction<ScalarFESType>;

      using VectorTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<VectorGridFunctionType, VectorFESType>;
      using ScalarTrialFunctionType =
        Rodin::PETSc::Variational::TrialFunction<ScalarGridFunctionType, ScalarFESType>;
      using VectorTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<VectorFESType>;
      using ScalarTestFunctionType =
        Rodin::PETSc::Variational::TestFunction<ScalarFESType>;

      using LinearSystemType = Rodin::PETSc::Math::LinearSystem;

      /// @brief (u, p, psi_D, psi_O; v, q, chi_D, chi_O). Each block of psi is
      ///        a three-component P1 field, so it shares the vector types.
      using FlowProblemType =
        Rodin::Variational::Problem<LinearSystemType,
          VectorTrialFunctionType, ScalarTrialFunctionType,
          VectorTrialFunctionType, VectorTrialFunctionType,
          VectorTestFunctionType, ScalarTestFunctionType,
          VectorTestFunctionType, VectorTestFunctionType>;

      using VectorProblemType = Rodin::Variational::Problem<LinearSystemType,
        VectorTrialFunctionType, VectorTestFunctionType>;

      using ScalarProblemType = Rodin::Variational::Problem<LinearSystemType,
        ScalarTrialFunctionType, ScalarTestFunctionType>;

      using FluxFormType = Rodin::Variational::LinearForm<ScalarFESType, ::Vec>;

      /// @brief Cellwise (P0) stabilisation parameters.
      using CellFESType = Rodin::Variational::P0<Real, MeshType>;
      using CellGridFunctionType = Rodin::PETSc::Variational::GridFunction<CellFESType>;

      /// @brief sPTT/Oldroyd-B blood, as in the 2D driver.
      struct OldroydB
      {
          Real etaS = 3.19e-3;  ///< solvent viscosity eta_s (Pa s)
          Real etaP = 4.0e-4;   ///< polymeric viscosity eta_p (Pa s)
          Real lambda = 0.1;   ///< relaxation time lambda (s)
          /// @brief lambda0 = max(k lambda, lambda0_min).
          Real lambda0Factor = 1.0;
          Real lambda0Min = 0.0;
          /// @brief Linear sPTT extensibility parameter; 0 recovers Oldroyd-B.
          Real pttEpsilon = 0.1;
      };

      /**
       * @brief Boundary labels written by make_lesion_mesh_3d.py.
       *
       * @details 1 inlet disk (x = -Lu D), 2 outlet disk (x = Ld D), 3
       *          parent-vessel wall, 5 lesion wall (|x| <= ell).
       */
      struct Labels
      {
          Attribute inlet = 1;
          Attribute outlet = 2;
          std::array<Attribute, 2> wall{{3, 5}};
          Attribute lesion = 5;
      };

      /// @brief Carreau-Yasuda viscosity of whole blood (Cho & Kensey 1991).
      struct CarreauYasuda
      {
          Real etaZero = 0.056;    ///< zero-shear viscosity (Pa s)
          Real etaInf = 0.00345;   ///< infinite-shear viscosity (Pa s)
          Real lambda = 1.902;     ///< time constant (s)
          Real n = 0.22;
          Real a = 1.25;
      };

      /**
       * @brief Owens (2006) blood model, JNNFM 140:57-70, single mode: red
       *        cells aggregate into rouleaux of mean size n (Hookean
       *        dumbbells, friction proportional to the size). Eqs. (23), (24),
       *        (10), (27), (28) of the paper:
       *          tau + mu tau^nabla = G mu gd,  G = N0 kB T = eta_p,inf/lambda_H,
       *          mu(n, gd) = n lambda_H / (1 + (b n (n - 1)/2 + aN0) lambda_H),
       *          D n/Dt = -(1/2) b (n - n_st)(n + n_st - 1),
       *          n_st = (eta(gd) - eta_s)/eta_p,inf (1 + (3/2) aN0 lambda_H),
       *          b = aN0/(n_st (n_st - 1)),
       *        so that G mu(n_st) = eta(gd) - eta_s exactly. eta(gd) is the
       *        Carreau-Yasuda of "cy" (Owens fits a Cross model: any steady
       *        fit will do, Eq. (29)); the plasma is the solvent, as in Fang &
       *        Owens (2006) and Moyers-Gonzalez et al. (JFM 2008), and eta_p,inf
       *        = eta_CY(inf) - eta_plasma.
       *        Aggregation rate aN0(gd): Moyers-Gonzalez, Owens & Fang (JFM
       *        2008, Part I, Appendix A, Eq. A1), the C1 cubics with zero slope
       *        at 0, gd_c = 5.78 1/s (Murata & Secomb) and gd_max = 900 1/s:
       *          aN0 = a0 + (ac - a0)(3 s^2 - 2 s^3),  s = gd/gd_c,          gd <= gd_c
       *          aN0 = ac (1 - 3 s^2 + 2 s^3),  s = (gd - gd_c)/(gd_max - gd_c), gd <= gd_max
       *          aN0 = 0,                                                     gd > gd_max,
       *        with a0 = a_{1,0} = 1 1/s and ac = 2 a0 (their coefficients
       *        a_{1,k}, a_{2,k} expand exactly to this). b -> 0 with aN0 near
       *        gd_max, which would freeze n in the jet; bFloor keeps b at its
       *        maximum above the shear rate where it peaks (steady viscosity
       *        changed by < 0.2 %, printed at start).
       */
      struct Owens
      {
          Real lambdaH = 0.005;      ///< single-cell relaxation time (s)
          Real etaPlasma = 1.2e-3;   ///< plasma (solvent) viscosity (Pa s)
          Real a0 = 1.0;             ///< aN0 at rest, a_{1,0} (1/s)
          Real aCrit = -1.0;         ///< aN0 at gd_c (1/s); < 0: 2 a0 (JFM 2008)
          Real gdCrit = 5.78;        ///< shear rate of the aggregation peak (1/s)
          Real gdMax = 900.0;        ///< aN0 = 0 above it (1/s)
          bool bFloor = true;
          /// @brief sPTT epsilon of the dumbbells (0: Hookean, Owens' model).
          Real epsilon = 0.0;
          /// @brief Initial aggregate size: < 0 the rest value n_st(0), 0 the
          ///        developed-profile value n_st(gd(r)), > 0 that value.
          Real nInitial = -1.0;
          /// @brief Inlet aggregate size: < 0 n_st(gd(r)) of the developed
          ///        profile (equilibrium); > 0 that value (blood that left the
          ///        heart ~1 s ago: n ~ 1 + aN0 t/2 ~ 1.5, aggregation being
          ///        ~10 s slower than the transit).
          Real nInlet = -1.0;
      };

      struct Config
      {
          /// @brief Discretisation of the constitutive equation and of its
          ///        stabilisation (file header). In all three, momentum and
          ///        continuity carry the split OSS with dynamic subscales.
          ///  - "psi_rb": the psi equation with the OSS of its full residual,
          ///    a log-conformation subscale psi~ seen by the psi equation
          ///    (streamline) and by the momentum (sigma~ = s Dexp[psi~]).
          ///  - "exp": Moreno, Codina, Baiges, Castillo (CMAME 2019) and Moreno,
          ///    Codina, Baiges (CMAME 2020). The constitutive equation is the
          ///    sigma equation over 2 eta_p written in E = exp(psi) and tested
          ///    with a stress-like chi; momentum and continuity carry the
          ///    term-by-term (split) OSS, with dynamic subscales u1
          ///    (convection), u2 (pressure gradient), u3 (div sigma) and the
          ///    algebraic pressure subscale; the constitutive equation carries
          ///    the OSS of its full residual, tested with the adjoint operator
          ///    A(chi, v), with a dynamic stress subscale. Newton is exact in
          ///    (u, p, psi): Dexp and D2exp by Daleckii-Krein.
          ///  - "psi": the former scheme, the psi equation d_t psi + u.grad psi
          ///    - Dexp^{-1}[L T + T L^T] + (f/lambda)(I - exp(-psi)) = 0 with the
          ///    term-by-term OSS of eps(u) and (u.grad) psi.
          std::string formulation = "psi_rb";
          /// @brief Stabilisation constants of "psi_rb" and "exp" (tau_psi =
          ///        (lambda/(2 eta_p)) alpha3), linear
          ///        elements (Moreno et al. 2019): alpha1 = [c1 eta_0/h1^2 + c2
          ///        rho |u|/h2]^-1, alpha2 = h1^2/(c1 alpha1), alpha3 = [c3 f/(2
          ///        eta_p) + c4 (lambda |u|/(2 eta_p h2) + lambda |grad u|/eta_p)]^-1,
          ///        h1 = |K|^(1/3), h2 the length of K along u.
          Real c1 = 4.0;
          Real c2 = 2.0;
          Real c3 = 4.0;
          Real c4 = 0.25;
          /// @brief Dynamic subscales (time-dependent, BDF1 in their own
          ///        equations): alpha1 -> (rho/dt + 1/alpha1)^-1, alpha3 ->
          ///        (lambda/(2 eta_p dt) + 1/alpha3)^-1, plus their memory terms.
          ///        Off: quasi-static subscales ("psi_rb" and "exp").
          bool dynamicSubscales = true;
          /// @brief "psi_rb": test the log-conformation subscale with the full
          ///        adjoint of the linearised psi operator, -(u^n.grad) chi +
          ///        JN^*[chi] (the reaction part JN^* = (dN/dpsi)^* at step n,
          ///        as in the stress adjoint of Moreno et al.); off (default):
          ///        the convective part only. JN carries the deformation terms
          ///        of N, of size |grad u| and of either sign, and tau |JN| =
          ///        O(1) in the stenosis throat: with it, Newton failed at t =
          ///        0.046 s in the S75 test (dt cut by 8), without it no step
          ///        was rejected to t = 0.115 s. In the sPTT pipe both give the
          ///        same errors.
          bool adjointReaction = false;
          /// @brief "psi_rb" and "exp": keep d_t of the FE field (BDF) in the
          ///        residual that drives the stress subscale (on, default);
          ///        off: without it, as in the OSS, where P^perp[d_t psi_h] = 0,
          ///        and as in the term-by-term form of Castillo, Moreno, Baiges,
          ///        Codina, ACME. Off was tried as the fix of the analysis
          ///        below and refuted (S75 R Re 600, xcoarse, 200 steps/cycle):
          ///        on ran the ramp cycle without a rejected step, off cut dt
          ///        to 2e-5 s with max|sigma| ~ 2.5e4 Pa from t = 0.157 s. "exp": s~ = alpha3' (memS s~^n - R + Pi), R ~
          ///        k_a (E - E_hist)/dtM, alpha3' -> 2 eta_p dt/lambda where
          ///        alpha3 >> that: the stress s (E - I) + s~ seen by the
          ///        momentum responds to E^{n+1} with s (1 - dt/dtM), -s/2 in
          ///        BDF2. "psi_rb": With it, psi~^{n+1} contains
          ///        -(tau'/dtM) psi_h^{n+1}, tau' = (1/dt + 1/tau)^-1 (BDF1 in
          ///        the subscale) and dtM = 2 dt/3 (BDF2): the stress seen by
          ///        the momentum, s Dexp[psi_h + psi~], responds to psi_h^{n+1}
          ///        with the factor 1 - tau'/dtM, which is -1/2 wherever tau >>
          ///        dt (slow flow, flow reversal, and everywhere as dt -> 0):
          ///        the u-psi coupling changes sign. The S75 Re 600 test grew
          ///        at ~8000 1/s from t = 0.135 s (deceleration) whatever dt.
          bool subscaleTimeResidual = true;
          /// @brief "sptt", "gsptt", "cy" or "owens" (file header, Owens).
          std::string fluid = "sptt";
          CarreauYasuda carreauYasuda;
          Owens owens;
          /// @brief Viscosity defining Re and Wo (Pa s); 0 means eta_s + eta_p
          ///        for "sptt" and eta_inf for "cy" and "gsptt".
          Real etaRef = 0.0;
          /// @brief Impose on the inlet the steady simple-shear conformation of
          ///        the developed mean profile instead of psi = 0 (requires
          ///        lambda0 = lambda).
          bool inletDevelopedConformation = false;
          /// @brief Scale that inlet log-conformation by the inflow ramp,
          ///        psi_in(t) = ramp(t) psi_dev, as the inlet velocity. The flow
          ///        starts from rest (psi = 0): without it psi jumps by psi_dev
          ///        across the first cell layer at t = 0, and in the "exp"
          ///        formulation the first Newton solve diverges at every dt.
          bool inletConformationRamp = true;

          /// @brief Tetrahedral MEDIT "Dimension 3" mesh in units of D.
          std::string meshPath =
            "../resources/examples/viscoelastic_fluids/S75_pipe_coarse.mesh";

          /// @brief Output paths; missing directories are created.
          std::string xdmfBasename = "results/ArterialLesion3D_viscous_logarithmic_implicit";
          std::string csvPath = "results/ArterialLesion3D_viscous_logarithmic_implicit.csv";

          Labels labels;

          /// @brief Parent-vessel diameter D (m); the mesh is scaled by it.
          Real diameter = 4.0e-3;

          Real rho = 1060.0;
          OldroydB oldroydB;
          /// @brief Pressure penalty of the equal-order pair.
          Real pressurePenalty = 1.0e-12;

          /// @brief Re = rho Ubar D/eta_0; sets the mean inlet velocity.
          Real reynolds = 300.0;
          /// @brief Wo = (D/2) sqrt(omega rho/eta_0); sets the period.
          Real womersley = 4.0;
          /// @brief q(t) = 1 + A sin(omega t); A = 0 is steady inflow.
          Real amplitude = 0.5;
          /// @brief Optional flow-rate shape (t, Q); replaces the sinusoid.
          std::string flowWaveformPath;
          /// @brief Fourier modes kept from the file.
          int harmonics = 10;

          /// @brief Wi = lambda Ubar/D. If > 0, it sets lambda.
          Real weissenberg = 0.0;
          /// @brief De = lambda/T. If > 0, it sets lambda (and wins over Wi).
          Real deborah = 0.0;

          /// @brief Inflow ramp over this many cycles, r(t) = (1 - cos(pi
          ///        t/T_r))/2. 0 disables it.
          Real rampCycles = 1.0;

          /// @brief Outlet traction level (Pa). Does not affect u.
          Real outletPressure = 0.0;
          Real outletBackflowStabilization = 1.0;

          int stepsPerCycle = 400;
          /// @brief Time step (s). If > 0, it sets stepsPerCycle = round(T/dt),
          ///        so that a cycle holds a whole number of steps, and the run
          ///        uses dt = T/stepsPerCycle.
          Real timeStep = 0.0;
          int cycles = 5;
          /// @brief Time integrator: 1 = BDF1 (backward Euler), 2 = BDF2.
          ///        BDF2 discretises d/dt by (3 x^{n+1} - 4 x^n + x^{n-1})/(2 dt)
          ///        in the momentum and psi equations, the Galerkin convection
          ///        advects with u* = 2 u^n - u^{n-1}, and Newton starts from
          ///        (u*, p*, psi*). The orthogonal projections of the split-OSS
          ///        stabilisation are explicit; see bdf2Extrapolate. The
          ///        stabilisation parameters, the streamline directions,
          ///        Pi[(grad u) u] and the dynamic subscale stay at step n. The
          ///        first step takes x^{-1} = x^0, exact for a start from rest
          ///        with the (1 - cos) ramp.
          int timeOrder = 1;
          /// @brief BDF2: the projections taken at the extrapolation (u*, p*,
          ///        psi*) instead of step n, as a bit mask: Pi[div u] 1,
          ///        Pi[grad p] 2, Pi[div sigma] 4, Pi[eps(u)] 8,
          ///        Pi[(u^n . grad) psi] 16, Pi[u^n . grad psi_t] 32
          ///        (axisymmetric). A projection lagged at step n leaves an
          ///        O(dt) error times its weight. The default 8 extrapolates
          ///        Pi[eps(u)], whose weight alpha3 ~ eta_p/2 does not vanish with
          ///        h. Extrapolating the others diverged: Pi[grad p], Pi[div u]
          ///        and Pi[div sigma] (masks 2, 3, 6) in the axisymmetric
          ///        pulsatile Oldroyd-B test after 2 to 3 periods, and the psi
          ///        advection projection (16) with the "gsptt" blood model within
          ///        the first period. Lagged, Pi[grad p] and Pi[div sigma]
          ///        (weights O(h^2)) leave an O(h^2 dt) error.
          int bdf2Extrapolate = 8;

          /// @brief XDMF output period, in steps; cycle boundaries always.
          int outputEvery = 20;

          /// @brief Per-term stabilisation scales, as in the 2D driver.
          Real vmsScale = 1.0;       ///< S1, dynamic convection
          Real gradDivScale = 1.0;   ///< S2, alpha2
          Real pressureScale = 1.0;  ///< S1, alpha1 <P^perp[grad p], grad q>
          Real stressDivScale = 1.0; ///< S1, alpha1 <P^perp[div sigma], div chi>
          Real stressScale = 1.0;    ///< S3, alpha3
          bool useVMS = true;

          /// @brief Newton on the coupled system.
          int conformationIterations = 6;
          Real conformationTolerance = 1.0e-6;
          Real newtonMaxStep = 2.0;
          /// @brief Impose psi = 0 (relaxed fluid) on the inlet.
          bool inletConformation = true;
          /// @brief Print the Newton history at every step, not only the first three.
          bool traceNewton = false;
          /// @brief Keep the D2exp part of the S1 div sigma linearisation
          ///        (quadratic Newton). Without it the converged solution is the
          ///        same, Newton converges linearly, and each iteration costs
          ///        less (four bilinear and two linear integrals fewer).
          bool newtonD2 = true;
          /// @brief Quadrature degree of the integrals that read the pointwise
          ///        log-conformation data (tetrahedra: 2 -> 4 points, 4 -> 11).
          int psiQuadratureOrder = 4;

          /// @brief Seed for non-axisymmetric instabilities: the inlet profile
          ///        is multiplied by 1 + delta y/R during the first
          ///        perturbationCycles cycles (the flow rate is unchanged, since
          ///        the integral of y over the inlet disk vanishes), then the
          ///        flow evolves freely. delta = 0 disables it.
          Real inletPerturbation = 0.0;
          Real perturbationCycles = 1.0;

          /// @brief Divergence guard: max|u| above this multiple of Ubar
          ///        rejects the step (stops the run without adaptiveTimeStep).
          Real maxVelocityFactor = 20.0;
          /// @brief Diagnostics: XDMF output at every accepted step with
          ///        max|sigma| > dumpStress (Pa; 0 off), at most dumpMax times;
          ///        the last accepted state is also written when the run aborts.
          Real dumpStress = 0.0;
          int dumpMax = 12;

          /// @brief Step control. A step is rejected, and retried from the same
          ///        x^n with dt * dtShrink, when Newton does not converge in
          ///        conformationIterations, when its undamped increment grows,
          ///        when the linear solver fails or when the divergence guard
          ///        trips. Nothing of x^n, x^{n-1} changes before a step is
          ///        accepted, so a retry starts from the same state. After
          ///        dtGrowAfter consecutive steps of at most newtonEasy
          ///        iterations dt grows by dtGrow, never above the nominal
          ///        T/stepsPerCycle; a step of newtonHard or more iterations
          ///        shrinks the next one. Steps are clipped to end on the cycle
          ///        boundaries; BDF2 is then variable-step. Iterations count
          ///        undamped ones only. Without it, a step
          ///        whose Newton does not converge is accepted with a warning and
          ///        any other failure stops the run (the former behaviour).
          bool adaptiveTimeStep = true;
          /// @brief Retries of one step; the smallest dt is dt0 dtShrink^maxStepRetries.
          int maxStepRetries = 6;
          Real dtShrink = 0.5;
          /// @brief Growth factor; BDF2 with variable steps is zero-stable for
          ///        dt_{n+1}/dt_n < 1 + sqrt(2), so it stays well below that.
          Real dtGrow = 1.5;
          int dtGrowAfter = 10;
          int newtonEasy = 2;
          int newtonHard = 5;

          /// @brief Newton acceptance by the quadratic prediction: after an
          ///        undamped increment d_k, the next one is ~ C d_k^2, with C
          ///        estimated from the steps that took two undamped iterations or
          ///        more. The iterate is accepted when C d_k^2 is below
          ///        newtonPredictSafety * conformationTolerance, which saves the
          ///        confirming iteration (one assembly and one solve). Every
          ///        newtonCalibrateEvery steps the prediction is not used, so
          ///        that C keeps being measured.
          bool newtonPredict = true;
          Real newtonPredictSafety = 0.1;
          int newtonCalibrateEvery = 20;

          /// @brief Lagged preconditioner of the flow solve: the preconditioner
          ///        (ILU blocks or LU factors) of an earlier matrix is reused for
          ///        up to pcReuseSolves solves, across Newton iterations and
          ///        steps. A solve with it that needs more than pcReuseGrowth
          ///        times the iterations of the first Newton solve after the
          ///        last build is kept, and the preconditioner is rebuilt for the
          ///        next solve; one that fails is repeated, from the same guess,
          ///        with a fresh preconditioner. A change of dt rebuilds it too.
          ///        1 rebuilds it at every solve. Only meaningful with a Krylov
          ///        method (-ksp_type gmres); with preonly it is ignored.
          int pcReuseSolves = 1;
          Real pcReuseGrowth = 2.0;

          /// @brief Weak inflow condition of the psi equation on the outlet:
          ///        where u*.n < 0 the term |u*.n|_- (psi - psi^n) : chi makes
          ///        the entering fluid carry the conformation found there at
          ///        step n. The psi equation is hyperbolic, and without it a
          ///        backflow on the outlet enters with no data.
          bool outletConformationInflow = true;
          /// @brief Skew-symmetric psi advection, + (1/2)(div u*) psi : chi: the
          ///        discrete div u is not zero for P1/P1, and this restores the
          ///        L2 stability of the advection of psi (Temam's term of the
          ///        momentum equation, applied to psi). Weakly consistent.
          bool psiSkewAdvection = false;
      };

      /// @brief Wall-clock cost of each phase of one step.
      struct Timing
      {
          Real vms = 0.0;
          Real assembly = 0.0;
          Real solve = 0.0;
          Real shear = 0.0;
          Real fluxes = 0.0;
          Real output = 0.0;
          Real total = 0.0;
          /// @brief Rejected attempts of the step, preconditioner builds and
          ///        Krylov iterations over all its attempts.
          int retries = 0;
          int pcBuilds = 0;
          int kspIterations = 0;
      };

      ArterialLesion3DViscousLogImplicit(const Rodin::Context::MPI& context, const Config& cfg);
      ~ArterialLesion3DViscousLogImplicit();

      ArterialLesion3DViscousLogImplicit(const ArterialLesion3DViscousLogImplicit&) = delete;
      ArterialLesion3DViscousLogImplicit& operator=(const ArterialLesion3DViscousLogImplicit&) = delete;

      ArterialLesion3DViscousLogImplicit& initialize();
      int run();

      Config& getConfig() noexcept
      {
        return m_cfg;
      }
      const Config& getConfig() const noexcept
      {
        return m_cfg;
      }

    private:
      static MeshType makeMesh(const Rodin::Context::MPI& context, const Config& cfg);
      static AttributeSet makeWallSet(const Config& cfg);

      bool isRoot() const;

      /// @brief Ubar, T and lambda from (Re, Wo, Wi or De); dt = T/steps.
      void deriveParameters();

      void setupSpaces();
      void setupFlow();
      /// @brief The "exp" formulation (Config::formulation): its forms.
      void setupFlowExp();
      /// @brief "exp": the alphas at step n (dynamic or quasi-static).
      void updateStabilizationExp();
      /// @brief "exp", after a step is accepted and before the history moves:
      ///        the subscales of step n+1 at the quadrature points, then
      ///        Pi[R^{n+1}], the projection of the constitutive residual read by
      ///        the next step.
      void updateSubscalesExp();
      /// @brief The "psi_rb" formulation (Config::formulation): its forms.
      void setupFlowPsiRB();
      /// @brief "psi_rb", as updateSubscalesExp: the subscales of step n+1 and
      ///        Pi[R_psi^{n+1}].
      void updateSubscalesPsiRB();
      /// @brief The length of the cell of p along the direction of u (h1 if
      ///        |u| vanishes).
      static Real streamlineLength(const Rodin::Geometry::Point& p,
        const Rodin::Math::SpatialVector<Real>& u);
      void setupWallShear();

      /// @brief L2 projection of a vector expression, reusing one mass matrix.
      template <class Expression>
      void projectVector(const Expression& expr, VectorGridFunctionType& out)
      {
        m_vectorProjection = Rodin::Variational::Integral(m_wTrial, m_wTest) -
          Rodin::Variational::Integral(expr, m_wTest);
        m_vectorProjection.solve(m_vectorProjectionKSP);
        out.setData(m_wTrial.getSolution().getData());
      }

      /// @brief y <- y + a x, with the ghost values refreshed.
      static void axpy(Real a, const ::Vec& x, ::Vec& y);

      /// @brief z <- a x + b y, with the ghost values refreshed.
      static void combine(Real a, const ::Vec& x, Real b, const ::Vec& y, ::Vec& z);

      /// @brief The time step of the BDF mass terms: dt (BDF1) or, for BDF2
      ///        with w = dt_{n+1}/dt_n, dt (1 + w)/(1 + 2 w) (2 dt/3 at w = 1).
      Real massStep() const;

      /// @brief w = dt_{n+1}/dt_n of the BDF2 step about to be taken; 1 on the
      ///        first step (x^{-1} = x^0) and for BDF1.
      Real bdfRatio() const;

      /// @brief Before each attempt of a step, from x^n, x^{n-1} and w: the BDF
      ///        history x_hist, with d/dt x ~ (x - x_hist)/massStep(),
      ///          x_hist = ((1 + w)^2 x^n - w^2 x^{n-1})/(1 + 2 w)  (BDF2),
      ///        and the extrapolation u* = (1 + w) u^n - w u^{n-1}.
      void prepareStep();

      /// @brief After a step is accepted (m_uIt, the psi iterates and the
      ///        pressure solution hold x^{n+1}): shifts x^n, x^{n-1}.
      void advanceHistory();

      /// @brief After a rejected attempt: the solution fields back to x^n.
      void rejectStep();

      /// @brief The orthogonal projections Pi[grad p], Pi[div sigma] and, if
      ///        strain, Pi[div u], Pi[eps(u)] and Pi[(u^n . grad) psi], at step n
      ///        (iterate = false) or at the Newton iterate (u^k, p^k, psi^k).
      void projectConsistent(bool iterate, unsigned mask);

      bool solveFlow();
      void computeWallShear();
      void closeCycle(Real elapsed);
      void computeFluxes();
      Real boundaryMeasure(const AttributeSet& tags);

      /// @brief r(t) = (1 - cos(pi t/T_r))/2 for t < T_r, then 1.
      Real ramp(Real t) const;

      /// @brief r(t) q(t): the ramped normalised flow rate.
      Real inflowFactor(Real t) const;

      /// @brief h_K, the local element size.
      static Real cellSize(const Rodin::Geometry::Point& p);

      /// @brief Codina tau_1 on eta_s, as in the 2D driver.
      Real tau1At(const Rodin::Geometry::Point& p) const;

      // ---- Fluids "gsptt" and "cy" -----------------------------------------
      Real carreauYasuda(Real gammaDot) const;
      /// @brief Steady-shear sPTT factor: f^3 - f^2 - 2 eps lambda^2 gd^2 = 0.
      Real steadyPTTFactor(Real gammaDot) const;
      /// @brief eta_s(gd): eta_s ("sptt"), eta_CY ("cy"), eta_CY - eta_p/f ("gsptt").
      Real solventViscosity(Real gammaDot) const;
      /// @brief gd = sqrt(2 eps:eps) of u^n (BDF1) or u* (BDF2).
      Real shearRateAt(const Rodin::Geometry::Point& p) const;
      Real solventViscosityAt(const Rodin::Geometry::Point& p) const;
      Real referenceViscosity() const;
      /// @brief The solvent viscosity does not depend on the shear rate.
      bool constantSolvent() const
      {
        return m_cfg.fluid == "sptt" || m_cfg.fluid == "owens";
      }
      bool isOwens() const
      {
        return m_cfg.fluid == "owens";
      }
      /// @brief Relaxation time of the steady simple shear at this rate:
      ///        lambda, or (eta_CY - eta_plasma)/G for Owens.
      Real steadyLambda(Real gammaDot) const;
      /// @brief Owens: aggregation rate aN0(gd), fragmentation rate b(gd),
      ///        steady aggregate size n_st(gd) and mu = lambda(n, gd).
      Real owensA(Real gammaDot) const;
      Real owensBRaw(Real gammaDot) const;
      Real owensB(Real gammaDot) const;
      Real owensNst(Real gammaDot) const;
      Real owensLambda(Real n, Real gammaDot) const;
      Real shearRateOf(const VectorGridFunctionType& u, const Rodin::Geometry::Point& p) const;
      /// @brief Owens: n_st at the shear rate of the developed inlet profile.
      Real inletAggregateSize(const Rodin::Geometry::Point& p) const;
      /// @brief Owens: lambda by cell for the psi equation, or nullptr.
      const CellGridFunctionType* lambdaField() const
      {
        return isOwens() ? &m_lambdaCell : nullptr;
      }
      void setupStructure();
      /// @brief Owens: n^{n+1} from the accepted (u^{n+1}, n^n), then
      ///        lambda(n^{n+1}, gd^{n+1}) for the next step.
      void solveStructure();
      void updateStructureFields();
      /// @brief Developed steady Carreau-Yasuda pipe flow at the mean flow rate:
      ///        m_inletR, m_inletPhi (u/Ubar) and m_inletShear (gd >= 0).
      void computeDevelopedProfile();
      /// @brief Linear interpolation of a table on m_inletR.
      Real interpolateInlet(const std::vector<Real>& table, Real r) const;
      /// @brief log c of steady simple shear at the inlet point p (radial shear
      ///        du_x/dr of the developed mean profile), D and O blocks.
      Rodin::Math::SpatialVector<Real> inletLogConformation(
        const Rodin::Geometry::Point& p, bool offDiagonal) const;
      /// @brief alpha3, with the PTT factor at psi^n.
      Real alpha3At(const Rodin::Geometry::Point& p) const;

      /// @brief The alphas (P0) and every Pi_h, all at step n.
      void updateStabilization();

      /// @brief lambda0 = max(k lambda, lambda0_min).
      Real lambda0() const;

      /// @brief PTT factor f = 1 + epsilon (lambda/lambda0)(tr exp(psi) - 3).
      Real pttFactor(Real traceExp) const;

      /// @brief sigma at the nodes, from (m_psiItD, m_psiItO).
      void updateConformation();

      void writeCSVHeader();
      void writeCSVRow(int cycle);

      Config m_cfg;

      PulsatilePipeInflow m_inflow;
      std::vector<Real> m_inletR;
      std::vector<Real> m_inletPhi;
      std::vector<Real> m_inletShear;
      /// @brief Mean inlet velocity Ubar (m/s) and period T (s).
      Real m_meanVelocity = 0.0;
      Real m_period = 0.0;
      /// @brief dt of the step being taken, its nominal value T/stepsPerCycle,
      ///        the dt of the last accepted step and the number of accepted steps.
      Real m_dt = 0.0;
      Real m_dtNominal = 0.0;
      Real m_dtPrev = 0.0;
      int m_history = 0;
      /// @brief (dt, massStep) the flow forms were built with: setupFlow bakes
      ///        both into the forms, so a change of either rebuilds them.
      Real m_formDt = -1.0;
      Real m_formMassStep = -1.0;

      MeshType m_mesh;
      Rodin::IO::XDMF m_xdmf;
      AttributeSet m_wallSet;

      VectorFESType m_vh;
      ScalarFESType m_sh;
      /// @brief Space of one block of psi: P1, three components.
      VectorFESType m_tauh;
      CellFESType m_ch;

      // ---- Flow ------------------------------------------------------------
      VectorTrialFunctionType m_u;
      ScalarTrialFunctionType m_p;
      VectorTestFunctionType m_v;
      ScalarTestFunctionType m_q;
      /// @brief u^n: the stabilisation parameters are evaluated at it.
      VectorGridFunctionType m_uOld;
      /// @brief Newton velocity u^k of the constitutive equation.
      VectorGridFunctionType m_uIt;
      /// @brief u^{n-1} (BDF2); the BDF history, u^n or (4 u^n - u^{n-1})/3, so
      ///        that d/dt u ~ (u - m_uHist)/dt_eff, dt_eff = dt or 2 dt/3; and the
      ///        advecting velocity of the Galerkin convection, u^n or u*.
      VectorGridFunctionType m_uNm1;
      VectorGridFunctionType m_uHist;
      VectorGridFunctionType m_uAdv;
      /// @brief p^n, p^{n-1} (BDF2) and the pressure iterate p^k.
      ScalarGridFunctionType m_pN;
      ScalarGridFunctionType m_pNm1;
      ScalarGridFunctionType m_pIt;
      /// @brief psi = log(tau): diagonal (xx, yy, zz) and off-diagonal
      ///        (xy, xz, yz) blocks; the test tensor chi likewise.
      VectorTrialFunctionType m_psiD;
      VectorTrialFunctionType m_psiO;
      VectorTestFunctionType m_chiD;
      VectorTestFunctionType m_chiO;
      VectorGridFunctionType m_psiOldD;
      VectorGridFunctionType m_psiOldO;
      /// @brief psi^{n-1} and the BDF history of psi, as for u.
      VectorGridFunctionType m_psiNm1D;
      VectorGridFunctionType m_psiNm1O;
      VectorGridFunctionType m_psiHistD;
      VectorGridFunctionType m_psiHistO;
      /// @brief Linearisation point psi^k of the Newton iterations.
      VectorGridFunctionType m_psiItD;
      VectorGridFunctionType m_psiItO;
      /// @brief Bumped whenever (u^k, psi^k) or (u^n, psi^n) change; keys the
      ///        cache of the pointwise log-conformation data.
      std::uint64_t m_psiRevision = 0;
      /// @brief Bumped when the step data (psi^n, u^n) change: the revision of
      ///        the pointwise data read at step n ("psi_rb" adjoint).
      std::uint64_t m_stepRevision = 0;
      /// @brief sigma = (eta_p/lambda0)(exp(psi) - I) at the nodes, by block.
      VectorGridFunctionType m_sigmaD;
      VectorGridFunctionType m_sigmaO;

      // ---- Reusable mass projection (wall-shear recovery) -------------------
      VectorTrialFunctionType m_wTrial;
      VectorTestFunctionType m_wTest;

      // ---- VMS, split OSS (updateStabilization) ----------------------------
      CellGridFunctionType m_tauK;    ///< dynamic convective tau_K
      CellGridFunctionType m_alpha1;  ///< tau_1/rho
      CellGridFunctionType m_alpha2;  ///< rho h^2/(4 tau_1)
      CellGridFunctionType m_alpha3;
      CellGridFunctionType m_alphaPsi; ///< 2 k^2 s alpha3, the psi advection weight
      Heart::OrthogonalProjection<VectorFESType> m_piConv;     ///< Pi[(grad u^n) u^n]
      Heart::OrthogonalProjection<VectorFESType> m_sub;        ///< dynamic subscale u'
      VectorGridFunctionType m_subOld;
      Heart::OrthogonalProjection<ScalarFESType> m_piDiv;      ///< Pi[div u^n]
      Heart::OrthogonalProjection<VectorFESType> m_piGradP;    ///< Pi[grad p^n]
      Heart::OrthogonalProjection<VectorFESType> m_piDivSigma; ///< Pi[div sigma^n]
      Heart::OrthogonalProjection<VectorFESType> m_piEpsD;     ///< Pi[eps(u^n)], D block
      Heart::OrthogonalProjection<VectorFESType> m_piEpsO;     ///< Pi[eps(u^n)], O block
      Heart::OrthogonalProjection<VectorFESType> m_piAdvPsiD;  ///< Pi[(u^n.grad) psi^n], D
      Heart::OrthogonalProjection<VectorFESType> m_piAdvPsiO;  ///< Pi[(u^n.grad) psi^n], O
      /// @brief "exp": Pi[R^n] of the constitutive residual (D and O blocks),
      ///        the static and dynamic alpha3 and the PTT factor f(psi^n), by cell.
      Heart::OrthogonalProjection<VectorFESType> m_piRD;
      Heart::OrthogonalProjection<VectorFESType> m_piRO;
      CellGridFunctionType m_alpha3Dyn;
      CellGridFunctionType m_pttN;
      /// @brief Owens: the aggregate size n (P1; n^n, n^{n-1} and its BDF
      ///        history) and lambda(n^n, gd(u^n)) by cell, read by the psi
      ///        equation and its stabilisation during step n + 1.
      ScalarTrialFunctionType m_n;
      ScalarTestFunctionType m_nTest;
      ScalarGridFunctionType m_nOld;
      ScalarGridFunctionType m_nNm1;
      ScalarGridFunctionType m_nHist;
      CellGridFunctionType m_lambdaCell;
      ScalarProblemType m_structure;
      Rodin::Solver::KSP m_structureKSP;
      Real m_nMin = 1.0;
      Real m_nMax = 1.0;
      Real m_lambdaMax = 0.0;
      int m_dumps = 0;
      /// @brief Owens: the maximum of b(gd) and where it is (bFloor).
      Real m_owensBPeak = 0.0;
      Real m_owensGdPeak = std::numeric_limits<Real>::infinity();


      // ---- Wall shear indices ----------------------------------------------
      /// @brief The three rows of 2 eps(u), recovered onto the nodes.
      VectorGridFunctionType m_wss;
      VectorGridFunctionType m_symRec0;
      VectorGridFunctionType m_symRec1;
      VectorGridFunctionType m_symRec2;
      VectorGridFunctionType m_netShear;
      ScalarGridFunctionType m_absShear;
      ScalarGridFunctionType m_shearMagnitude;
      ScalarGridFunctionType m_tawss;
      ScalarGridFunctionType m_osi;
      /// @brief TAWSS and OSI on the lesion wall (tag 5) only.
      ScalarGridFunctionType m_tawssLesion;
      ScalarGridFunctionType m_osiLesion;

      // ---- Fluxes ----------------------------------------------------------
      ScalarTestFunctionType m_qFlux;
      ScalarGridFunctionType m_one;
      FluxFormType m_flux;
      /// @brief "exp": the sweep that writes the subscales at the quadrature points.
      FluxFormType m_subscaleSweep;

      // ---- Problems and solvers --------------------------------------------
      FlowProblemType m_flow;
      Rodin::Solver::KSP m_flowKSP;
      VectorProblemType m_vectorProjection;
      Rodin::Solver::KSP m_vectorProjectionKSP;
      VectorTrialFunctionType m_wssTrial;
      VectorTestFunctionType m_wssTest;
      VectorProblemType m_wssProjection;
      Rodin::Solver::KSP m_wssKSP;

      // ---- Time-dependent coefficients read by the frozen forms ------------
      Real m_t = 0.0;
      /// @brief r(t) q(t) at the current step, read by the inlet profile.
      Real m_inflowNow = 0.0;
      Real m_inletMeasure = 1.0;
      Real m_outletMeasure = 1.0;
      Real m_qIn = 0.0;
      Real m_qOut = 0.0;
      Real m_inletPressure = 0.0;
      Real m_outletMeanPressure = 0.0;
      /// @brief Volume of the domain and the volume means of u_y, u_z over
      ///        Ubar: zero for any axisymmetric flow, so their magnitude
      ///        measures the azimuthal mode m = 1 (tilted or deflected jet).
      Real m_volume = 1.0;
      Real m_uyMean = 0.0;
      Real m_uzMean = 0.0;
      Real m_speed = 0.0;
      /// @brief max |sigma_ij| (Pa), reported with the flow.
      Real m_stress = 0.0;
      /// @brief Newton iterations of the last step and their last increment.
      int m_conformationIts = 0;
      Real m_psiIncrement = 0.0;
      /// @brief Undamped Newton iterations of the last step: the measure of
      ///        its difficulty for the step control (the damped ones of the
      ///        start from rest are not).
      int m_newtonUndamped = 0;
      /// @brief Why the last attempt was rejected (empty if accepted).
      std::string m_failure;
      /// @brief Newton accepted by the quadratic prediction (not by d_k < tol).
      bool m_newtonPredicted = false;
      /// @brief Recent estimates of C in d_{k+1} ~ C d_k^2 (ring buffer).
      std::vector<Real> m_contraction;
      size_t m_contractionNext = 0;
      /// @brief Lagged preconditioner: built at least once, solves since its
      ///        last build, and a forced rebuild (dt changed).
      bool m_pcBuilt = false;
      int m_pcSolves = 0;
      /// @brief Krylov iterations of the first Newton solve after the last
      ///        build, the reference of pcReuseGrowth (< 0 until measured).
      int m_pcReferenceIts = -1;
      bool m_pcStale = true;
      bool m_flowFieldSplitsSet = false;
      bool m_initialized = false;

      int m_step = 0;
      Timing m_timing;

      std::ofstream m_csv;
  };
}

#endif
