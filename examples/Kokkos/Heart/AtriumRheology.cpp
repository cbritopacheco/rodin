/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file AtriumRheology.cpp
 * @brief Pressure-driven flow of blood in the left atrium on a fixed domain,
 *        Newtonian, Carreau-Yasuda, viscoelastic and thixotropic-viscoelastic
 *        (Owens), assembled and solved on the device.
 *
 * The four fluids of examples/Heart/Atrium/Viscoelastic (model/Atrium) with
 * the openings, the lumped ventricle and the time loop of AtriumFlow.cpp. The
 * three non-Newtonian fluids share the steady shear viscosity eta_CY(gd):
 *
 *  - N    constant viscosity mu;
 *  - CY   Carreau-Yasuda, the constants of Atrium/Fluid;
 *  - VE   Owens with the aggregate size at equilibrium, n = n_st(gd): the
 *         relaxation time lambda = (eta_CY(gd) - eta_plasma)/G follows the
 *         local shear rate, elasticity without structural memory;
 *  - TVE  Owens (JNNFM 140, 2006): n is transported and relaxes,
 *           D n/Dt = -(1/2) b (n - n_st)(n + n_st - 1),
 *           lambda(n, gd) = n lambda_H / (1 + (b n (n - 1)/2 + a) lambda_H),
 *         with a(gd), b(gd), n_st(gd), G = (eta_inf - eta_plasma)/lambda_H as
 *         in the CPU driver, so that G lambda(n_st) = eta_CY - eta_plasma.
 *
 * Momentum and continuity: P1/P1, BDF2 (order=1: backward Euler), convection
 * by u* = 2 u^n - u^{n-1} in skew-symmetric form, stabilised by the algebraic
 * subgrid scale method with quasi-static subscales (stab=asgs: SUPG and PSPG
 * with the residual of the momentum equation without its time derivative),
 *   tau (rho u*.grad u + grad p - div E, rho u*.grad v + grad q)/rho,
 *   tau = [(2|u*|/h)^2 + (4 nu/h^2)^2]^(-1/2) by cell,
 * or by its two diagonal terms alone (stab=split: three matrices less, not
 * consistent). tau has no (2/dt)^2 term: with it the stabilised problem
 * changes with the time step and the time discretisation does not converge
 * on a mesh that needs the stabilisation (Bochev, Gunzburger & Shadid, CMAME
 * 193, 2004). For CY, eta = eta_CY(gd*) and the transposed part of 2 eta D is
 * explicit, E = eta (L*^T - tr(L*) I) with L* the gradient of u* at the
 * vertices.
 *
 * Elastic stress (VE, TVE). tau_p = G (c - I), and psi = log c is stored at
 * the vertices (Fattal & Kupferman, JNNFM 123, 2004), so c stays positive
 * definite. Each step, from u^{n-1}, u^n, psi^n, n^n:
 *
 *  1. L^n, the gradient of u^n at the vertices: the volume-weighted mean of
 *     the cell gradients, made trace free, because the divergence of a
 *     stabilised P1 velocity is not zero and would pump the isotropic part
 *     of c. At a vertex of the wall it is reduced to (du/dn) n^T with du/dn
 *     tangent, which is what the gradient is on a no-slip wall: a simple
 *     shear, for which c is bounded. The mean of the cells of a wall vertex
 *     has an extensional part, and that vertex is not transported: c would
 *     grow there without bound.
 *  2. Advection: psi and n are taken at the foot of the characteristic of
 *     every vertex, X = x - dt u^{n+1/2}, found by walking through the cells
 *     and interpolated linearly there (Lee & Xu, CMAME 195, 2006). The new
 *     value is a convex combination of old ones: no over- or undershoot of
 *     psi to be exponentiated, the vertices of the wall (u = 0) are not
 *     transported, and there is no linear system. A foot beyond an opening
 *     takes the inflow value there: on the veins developed blood, the steady
 *     conformation and aggregate size of the shear of the axial velocity
 *     across the section (inlet=relaxed: psi = 0, n = n_in); on the mitral
 *     plane the value found there. A foot beyond the wall takes the value on
 *     the wall. First order in space, as upwinding.
 *  3. Prediction with L^n, at every vertex: n by the exact solution of its
 *     Riccati equation with b and n_st frozen, c by the exact solution of
 *       dc/dt = L c + c L^T - (c - I)/lambda
 *     over dt (matrix exponential of L, Simpson for the source), psi = log c:
 *     no bound on dt/lambda or dt |L|.
 *  4. Momentum and continuity for (u^{n+1}, p^{n+1}), with the predicted
 *     stress in weak form and two terms that vanish for a resolved steady
 *     flow:
 *       ... + (eta_plasma grad u + tau_pred, grad v)
 *           + (eta_a (grad u - L^n), grad v) + (eta_d (grad u - grad u^n), grad v).
 *     eta_a = alpha G lambda, with L^n the vertex gradient of step 1, is
 *     DEVSS-G (Guenette & Fortin, JNNFM 60, 1995; Liu, Bornside, Armstrong &
 *     Brown, JNNFM 77, 1998; the both-sides diffusion of finite volumes). A
 *     stress built at the vertices from mean gradients does not see the
 *     velocity modes whose mean gradient is zero, which only the plasma
 *     viscosity would resist; the term gives them the polymer viscosity. At
 *     low Weissenberg number tau_pred - eta_a L^n = G lambda L^T, and the
 *     discrete problem is that of CY. alpha is the key devss (1).
 *     eta_d = beta G dt_e lambda_max(c_pred), dt_e = lambda (1 - exp(-dt/lambda)),
 *     with the cell gradient of u^n, is the response of the stress to the
 *     change of velocity within the step, taken implicitly (D'Avino & Hulsen,
 *     JNNFM 165, 2010): G dt_e (L c + c L^T) is the change of the stress for
 *     a change L of the gradient, and eta_d bounds it. Without it the
 *     stress is explicit, and stable only while G lambda_max(c) (pi dt/h)^2
 *     /rho < 1: tr(c) < 30 for cells of 0.45 mm at dt = 1 ms, and not even
 *     c = I for cells of 0.06 mm (Coarse.mesh has both). beta is the key
 *     implicit (2).
 *  5. Correction: step 3 again from the advected state, with the gradient of
 *     (u^n + u^{n+1})/2.
 *
 * The elastic part is first order in time (the splitting of steps 2 and 3)
 * and in space (step 2); it is not the monolithic Newton solve with
 * orthogonal subscales of the CPU driver. Owens' dumbbells are Hookean: at a
 * stagnation point where lambda times the extension rate exceeds 1/2 their
 * stress grows without bound (tr(c)_max in the log); ptt > 0 divides lambda
 * by 1 + ptt (tr c - 3) (linear Phan-Thien-Tanner), which bounds it and is
 * not part of Owens' model.
 *
 * The CSV column timeDissipation is the kinetic energy that the time scheme
 * removes per unit time, rho/(4 dt) |u^{n+1} - 2 u^n + u^{n-1}|^2 (order=1:
 * rho/(2 dt) |u^{n+1} - u^n|^2), to be set against the viscous and polymer
 * dissipation.
 *
 * Openings, lumped ventricle and mesh: as AtriumFlow.cpp. The vertex kernels
 * (steps 1-3, 5 and the stresses) are Kokkos kernels of this file; the flow
 * problem is a PETSc matrix of the Kokkos type solved by PETSc.
 *
 * Usage:
 *   KokkosAtriumRheologyPETSc mesh.mesh [key=value ...] [PETSc options]
 * Keys: fluid (TVE; N, CY, VE, TVE), stab (asgs; split), order (2; 1), cycles
 * (2), steps (850, per cycle), t_end (stop at that time), period (0.85), rho
 * (1060), mu (0.0049, fluid N), p_in (1200), p_in_amplitude (200), ramp
 * (0.25, cycles), z_in (1000), lv (1), p_out (0), inlet (developed; relaxed),
 * n_in (1.5, with inlet=relaxed), n_init (-1: n_st(0)), devss (1), implicit
 * (2), ptt (0), eta_plasma (0.0012), lambda_h (0.005), output (AtriumRheology),
 * output_every (17), tolerance (1e-8), max_iterations (5000), pc_lag (1:
 * steps between two set-ups of the preconditioner), restart (none).
 *
 * Checkpoints: at the end of every cycle the run writes <output>.restart
 * (u^n, u^{n-1}, psi, n, the state of the ventricle and the step). A run with
 * restart=<file> on the same mesh and time step continues from it: time, cycle
 * numbers and the ventricle go on, and cycles / t_end count from there. psi and
 * n are taken from the file only if it was written by the same fluid; otherwise
 * they start from c = I and n_init, so that every fluid can start from one
 * developed flow (e.g. a Newtonian spin-up).
 * PETSc options: the flow solver has no prefix; BiCGStab(4) with Jacobi by
 * default, and a step whose solve fails is solved again with GMRES(100)
 * before the run stops. The unknowns are split in u and p for -pc_type
 * fieldsplit; with a preconditioner that costs to set up (an algebraic
 * multigrid on the pressure Schur complement), pc_lag keeps it for that
 * number of steps.
 * Output: <output>.xdmf and <output>.csv, one row per step.
 */
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <cstdint>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Heart/CCMLC2014.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Kokkos.h>
#include <Rodin/Kokkos/PETSc.h>
#include <Rodin/Kokkos/Solver.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

#if defined(KOKKOS_ENABLE_CUDA)
using Backend = Rodin::Context::CUDA;
#elif defined(KOKKOS_ENABLE_HIP)
using Backend = Rodin::Context::HIP;
#elif defined(KOKKOS_ENABLE_SYCL)
using Backend = Rodin::Context::SYCL;
#else
using Backend = ::Kokkos::DefaultHostExecutionSpace;
#endif

// A global using directive on Rodin would make Kokkos ambiguous in the Kokkos headers.
using DeviceContext = Rodin::Context::Kokkos<Rodin::Context::Local, Backend>;
using Ventricle = Rodin::Heart::CCMLC2014T<>;

namespace
{
  using Real = Rodin::Real;
  using Exec = ::Kokkos::DefaultExecutionSpace;
  using Memory = Exec::memory_space;
  using RealView = ::Kokkos::View<Real*, Memory>;
  using IntView = ::Kokkos::View<int*, Memory>;

  /// @brief Fluid model.
  enum class Fluid
  {
    Newtonian,
    CarreauYasuda,
    Viscoelastic,
    Thixotropic
  };

  /// @brief Run parameters, set by key=value arguments.
  struct Parameters
  {
      Fluid fluid = Fluid::Thixotropic;
      std::string fluidName = "TVE";
      bool residual = true; ///< stab=asgs
      int order = 2;
      Real scale = 1e-3;
      Real period = 0.85;
      size_t cycles = 2;
      size_t steps = 850;
      Real endTime = -1;
      Real rho = 1060;
      Real mu = 0.0049;
      Real inletPressureMean = 1200;
      Real inletPressureAmplitude = 200;
      Real rampCycles = 0.25;
      Real inletImpedance = 1000;
      bool ventricle = true;
      Real outletPressure = 0;
      bool developedInlet = true;
      Real inletSize = 1.5;
      Real initialSize = -1;
      Real devss = 1;
      Real implicit = 2;
      Real ptt = 0;
      Real etaPlasma = 1.2e-3;
      Real lambdaH = 0.005;
      std::string output = "AtriumRheology";
      std::string restart;
      size_t outputEvery = 17;
      Real tolerance = 1e-8;
      size_t maxIterations = 5000;
      size_t preconditionerLag = 1;
  };

  /**
   * @brief Constants and material functions of the fluids, usable in device
   *        kernels. Carreau-Yasuda of Atrium/Fluid; Owens as in
   *        examples/Heart/Atrium/Viscoelastic (Config::Owens).
   */
  struct Rheology
  {
      Real etaZero = 0.0720;
      Real etaInf = 0.0049;
      Real cyLambda = 6.1529;
      Real cyN = 0.2081;
      Real cyA = 1.4173;
      Real lambdaH = 0.005;
      Real etaPlasma = 1.2e-3;
      Real a0 = 1.0;
      Real aCrit = 2.0;
      Real gdCrit = 5.78;
      Real gdMax = 900.0;
      Real bPeak = 0;      ///< Floor of b above gdPeak, where b(gd) peaks.
      Real gdPeak = 1e300;
      Real etaPolymer = 0; ///< eta_inf - eta_plasma
      Real modulus = 0;    ///< G = etaPolymer/lambdaH

    /// @param gd Shear rate. @return Carreau-Yasuda viscosity.
      KOKKOS_INLINE_FUNCTION Real viscosity(Real gd) const
      {
        return etaInf +
          (etaZero - etaInf) *
          ::Kokkos::pow(1 + ::Kokkos::pow(cyLambda * gd, cyA), (cyN - 1) / cyA);
      }

    /// @param gd Shear rate. @return Aggregation rate a N0 (1/s).
      KOKKOS_INLINE_FUNCTION Real aggregation(Real gd) const
      {
        if (gd <= gdCrit)
        {
          const Real s = gd / gdCrit;
          return a0 + (aCrit - a0) * s * s * (3 - 2 * s);
        }
        if (gd <= gdMax)
        {
          const Real s = (gd - gdCrit) / (gdMax - gdCrit);
          return aCrit * (1 - s * s * (3 - 2 * s));
        }
        return 0;
      }

    /// @param gd Shear rate. @return Steady aggregate size n_st.
      KOKKOS_INLINE_FUNCTION Real steadySize(Real gd) const
      {
        Real r = (viscosity(gd) - etaPlasma) / etaPolymer;
        if (r < 1)
          r = 1;
        return r * (1 + 1.5 * aggregation(gd) * lambdaH);
      }

    /// @param gd Shear rate. @return Fragmentation rate b, without its floor.
      KOKKOS_INLINE_FUNCTION Real fragmentationRaw(Real gd) const
      {
        const Real a = aggregation(gd);
        const Real n = steadySize(gd);
        return (a > 0 && n > 1) ? a / (n * (n - 1)) : Real(0);
      }

    /// @param gd Shear rate. @return Fragmentation rate b (1/s).
      KOKKOS_INLINE_FUNCTION Real fragmentation(Real gd) const
      {
        const Real b = fragmentationRaw(gd);
        return (gd >= gdPeak && b < bPeak) ? bPeak : b;
      }

    /// @param n Aggregate size. @param gd Shear rate. @return Relaxation time.
      KOKKOS_INLINE_FUNCTION Real relaxation(Real n, Real gd) const
      {
        if (n < 1)
          n = 1;
        const Real gn = 0.5 * fragmentation(gd) * n * (n - 1) + aggregation(gd);
        return n * lambdaH / (1 + gn * lambdaH);
      }

    /// @param gd Shear rate. @return Relaxation time at n = n_st(gd).
      KOKKOS_INLINE_FUNCTION Real steadyRelaxation(Real gd) const
      {
        Real e = viscosity(gd) - etaPlasma;
        if (e < 1e-12 * etaPolymer)
          e = 1e-12 * etaPolymer;
        return e / modulus;
      }
  };

  /// @brief Three by three matrix of a device kernel.
  struct Tensor
  {
      Real a[3][3];
  };

  /// @return Identity.
  KOKKOS_INLINE_FUNCTION Tensor identity()
  {
    Tensor r;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        r.a[i][j] = i == j ? Real(1) : Real(0);
    }
    return r;
  }

  /// @param x Left factor. @param y Right factor. @return x y.
  KOKKOS_INLINE_FUNCTION Tensor multiply(const Tensor& x, const Tensor& y)
  {
    Tensor r;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        r.a[i][j] = x.a[i][0] * y.a[0][j] + x.a[i][1] * y.a[1][j] + x.a[i][2] * y.a[2][j];
    }
    return r;
  }

  /// @param x Left factor. @param y Right factor. @return x y^T.
  KOKKOS_INLINE_FUNCTION Tensor multiplyTranspose(const Tensor& x, const Tensor& y)
  {
    Tensor r;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        r.a[i][j] = x.a[i][0] * y.a[j][0] + x.a[i][1] * y.a[j][1] + x.a[i][2] * y.a[j][2];
    }
    return r;
  }

  /// @param m Matrix of norm below one. @return exp(m), Taylor series of order 10.
  KOKKOS_INLINE_FUNCTION Tensor exponential(const Tensor& m)
  {
    Tensor p = identity();
    for (int k = 10; k >= 1; --k)
    {
      const Tensor mp = multiply(m, p);
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
          p.a[i][j] = (i == j ? Real(1) : Real(0)) + mp.a[i][j] / k;
      }
    }
    return p;
  }

  /**
   * @brief Eigenvalues and eigenvectors of a symmetric matrix, cyclic Jacobi.
   * @param[in] s Symmetric matrix.
   * @param[out] w Eigenvalues.
   * @param[out] r Eigenvectors by column: s = r diag(w) r^T.
   */
  KOKKOS_INLINE_FUNCTION void eigen(const Tensor& s, Real w[3], Tensor& r)
  {
    Tensor a = s;
    r = identity();
    for (int sweep = 0; sweep < 16; ++sweep)
    {
      const Real off =
        a.a[0][1] * a.a[0][1] + a.a[0][2] * a.a[0][2] + a.a[1][2] * a.a[1][2];
      const Real diag =
        a.a[0][0] * a.a[0][0] + a.a[1][1] * a.a[1][1] + a.a[2][2] * a.a[2][2];
      if (off <= 1e-30 * diag || off < 1e-290)
        break;
      for (int p = 0; p < 2; ++p)
      {
        for (int q = p + 1; q < 3; ++q)
        {
          const Real apq = a.a[p][q];
          if (::Kokkos::fabs(apq) < 1e-300)
            continue;
          const Real theta = (a.a[q][q] - a.a[p][p]) / (2 * apq);
          const Real t = (theta >= 0 ? Real(1) : Real(-1)) /
            (::Kokkos::fabs(theta) + ::Kokkos::sqrt(theta * theta + 1));
          const Real c = 1 / ::Kokkos::sqrt(t * t + 1);
          const Real sn = t * c;
          for (int k = 0; k < 3; ++k)
          {
            const Real akp = a.a[k][p], akq = a.a[k][q];
            a.a[k][p] = c * akp - sn * akq;
            a.a[k][q] = sn * akp + c * akq;
          }
          for (int k = 0; k < 3; ++k)
          {
            const Real apk = a.a[p][k], aqk = a.a[q][k];
            a.a[p][k] = c * apk - sn * aqk;
            a.a[q][k] = sn * apk + c * aqk;
          }
          for (int k = 0; k < 3; ++k)
          {
            const Real rkp = r.a[k][p], rkq = r.a[k][q];
            r.a[k][p] = c * rkp - sn * rkq;
            r.a[k][q] = sn * rkp + c * rkq;
          }
        }
      }
    }
    for (int k = 0; k < 3; ++k)
      w[k] = a.a[k][k];
  }

  /// @param w Eigenvalues. @param r Eigenvectors by column. @return r diag(w) r^T.
  KOKKOS_INLINE_FUNCTION Tensor spectral(const Real w[3], const Tensor& r)
  {
    Tensor s;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        s.a[i][j] = r.a[i][0] * w[0] * r.a[j][0] + r.a[i][1] * w[1] * r.a[j][1] +
          r.a[i][2] * w[2] * r.a[j][2];
      }
    }
    return s;
  }

  /// @param psi Components xx, yy, zz, xy, xz, yz. @param largest Largest
  /// eigenvalue of exp(psi). @return exp(psi).
  KOKKOS_INLINE_FUNCTION Tensor conformation(const Real psi[6], Real& largest)
  {
    Tensor s;
    s.a[0][0] = psi[0];
    s.a[1][1] = psi[1];
    s.a[2][2] = psi[2];
    s.a[0][1] = s.a[1][0] = psi[3];
    s.a[0][2] = s.a[2][0] = psi[4];
    s.a[1][2] = s.a[2][1] = psi[5];
    Real w[3];
    Tensor r;
    eigen(s, w, r);
    for (int k = 0; k < 3; ++k)
      w[k] = ::Kokkos::exp(w[k]);
    largest = w[0] > w[1] ? (w[0] > w[2] ? w[0] : w[2]) : (w[1] > w[2] ? w[1] : w[2]);
    return spectral(w, r);
  }

  /// @param psi Components xx, yy, zz, xy, xz, yz. @return exp(psi).
  KOKKOS_INLINE_FUNCTION Tensor conformation(const Real psi[6])
  {
    Real largest;
    return conformation(psi, largest);
  }

  /**
   * @brief Advances psi = log c over dt by the exact solution of
   *        dc/dt = L c + c L^T - (c - I)/lambda with L and lambda frozen.
   *
   * With P = exp(L h/2), E = P P and d = exp(-h/(2 lambda)), a substep h is
   *   c <- d^2 E c E^T + h/(6 lambda) (I + 4 d P P^T + d^2 E E^T),
   * Simpson's rule for the source. dt is cut so that |L| h and h/lambda stay
   * below 1/2.
   *
   * @param[in,out] psi Components xx, yy, zz, xy, xz, yz.
   * @param[in] L Velocity gradient, L_ij = d u_i/d x_j.
   * @param[in] lambda Relaxation time.
   * @param[in] dt Time step.
   * @param[in] epsilon Linear PTT factor: lambda/(1 + epsilon (tr c - 3)) with
   *            the trace at the start of the step; 0 for Hookean dumbbells.
   */
  KOKKOS_INLINE_FUNCTION void advance(
    Real psi[6], const Tensor& L, Real lambda, Real dt, Real epsilon)
  {
    Tensor c = conformation(psi);
    const Real stretch = c.a[0][0] + c.a[1][1] + c.a[2][2] - 3;
    if (epsilon > 0 && stretch > 0)
      lambda /= 1 + epsilon * stretch;
    Real norm = 0;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        norm += L.a[i][j] * L.a[i][j];
    }
    norm = ::Kokkos::sqrt(norm);
    const Real load = 2 * dt * (norm > 1 / lambda ? norm : 1 / lambda);
    const int substeps = load < 1 ? 1 : (load > 64 ? 64 : static_cast<int>(load) + 1);
    const Real h = dt / substeps;
    Tensor half;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        half.a[i][j] = 0.5 * h * L.a[i][j];
    }
    const Tensor P = exponential(half);
    const Tensor E = multiply(P, P);
    const Tensor PPt = multiplyTranspose(P, P);
    const Tensor EEt = multiplyTranspose(E, E);
    const Real d1 = ::Kokkos::exp(-h / (2 * lambda));
    const Real d2 = d1 * d1;
    const Real weight = h / (6 * lambda);
    for (int s = 0; s < substeps; ++s)
    {
      const Tensor EcEt = multiplyTranspose(multiply(E, c), E);
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          c.a[i][j] = d2 * EcEt.a[i][j] +
            weight *
              ((i == j ? Real(1) : Real(0)) + 4 * d1 * PPt.a[i][j] + d2 * EEt.a[i][j]);
        }
      }
    }
    Tensor sym;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        sym.a[i][j] = 0.5 * (c.a[i][j] + c.a[j][i]);
    }
    Real w[3];
    Tensor r;
    eigen(sym, w, r);
    for (int k = 0; k < 3; ++k)
      w[k] = ::Kokkos::log(w[k] > 1e-300 ? w[k] : Real(1e-300));
    const Tensor l = spectral(w, r);
    psi[0] = l.a[0][0];
    psi[1] = l.a[1][1];
    psi[2] = l.a[2][2];
    psi[3] = l.a[0][1];
    psi[4] = l.a[0][2];
    psi[5] = l.a[1][2];
  }

  /**
   * @brief Vertex and cell data of the rheology, in device memory. Vector
   *        fields are stored as the P1 spaces of Rodin, component r of vertex
   *        i at i + r nodes.
   */
  struct Nodal
  {
      int nodes = 0;
      int cells = 0;
      IntView cellVertices;   ///< Four vertices per cell.
      IntView cellNeighbours; ///< Cell across the face opposite to each vertex;
                              ///< -1 beyond the wall, -2 beyond an opening.
      RealView coordinates;
      IntView nodeOffsets; ///< Cells of a vertex, compressed rows.
      IntView nodeCells;
      RealView volume;     ///< Lumped volume of a vertex.
      RealView wallNormal; ///< Unit normal at the vertices of the wall, zero elsewhere.
      RealView veinNormal; ///< Unit normal at the vertices of the veins, zero elsewhere.
      RealView cellGradient; ///< By cell, L_ij at 3 i + j and the volume at 9.
      RealView velocity; ///< Velocity read by recoverGradient() and advect().
      RealView gradient; ///< By vertex, L_ij at 3 i + j.
      RealView shearRate;
      RealView viscosity; ///< Implicit viscosity of the momentum equation.
      RealView increment; ///< eta_d, the part of it that acts on u^{n+1} - u^n.
      RealView relaxation;
      RealView size; ///< Aggregate size n.
      RealView psi; ///< Components xx, yy, zz, xy, xz, yz.
      RealView sizeAdvected; ///< n and psi at the feet of the characteristics.
      RealView psiAdvected;
      RealView sizePredicted; ///< n and psi predicted with the gradient of u^n.
      RealView psiPredicted;
      RealView sizeIn; ///< Values taken by a foot beyond an opening.
      RealView psiIn;
      RealView trace;
      RealView stress; ///< Frobenius norm of tau_p.
      RealView extra0, extra1, extra2; ///< Rows of the explicit stress E.
      RealView density; ///< Four diagnostic densities per vertex.
  };

  /// @brief Volume integrals and maxima of the vertex fields.
  struct Diagnostics
  {
      Real solventDissipation = 0;
      Real polymerPower = 0;
      Real elasticEnergy = 0;
      Real meanSize = 0;
      Real maxWeissenberg = 0;
      Real maxTrace = 0;
  };

  /// @brief Uploads host values. @return Device copy.
  template <class T>
  ::Kokkos::View<T*, Memory> upload(const std::string& name, const T* data, size_t n)
  {
    ::Kokkos::View<T*, Memory> device(
      ::Kokkos::view_alloc(name, ::Kokkos::WithoutInitializing), n);
    ::Kokkos::View<const T*, ::Kokkos::HostSpace, ::Kokkos::MemoryUnmanaged> host(
      data, n);
    ::Kokkos::deep_copy(device, host);
    return device;
  }

  /// @brief Copies host values to a device view of the same size.
  void toDevice(const Real* data, const RealView& device)
  {
    ::Kokkos::View<const Real*, ::Kokkos::HostSpace, ::Kokkos::MemoryUnmanaged> host(
      data, device.extent(0));
    ::Kokkos::deep_copy(device, host);
  }

  /// @brief Copies a device view to the host. @return Host copy.
  Rodin::Math::Vector<Real> toHost(const RealView& device)
  {
    Rodin::Math::Vector<Real> res(device.extent(0));
    ::Kokkos::View<Real*, ::Kokkos::HostSpace, ::Kokkos::MemoryUnmanaged> host(
      res.data(), res.size());
    ::Kokkos::deep_copy(host, device);
    return res;
  }

  /**
   * @brief Step 1: velocity gradient at the vertices, from Nodal::velocity.
   *
   * The volume-weighted mean of the cell gradients, made trace free; at the
   * vertices of the wall, (du/dn) n^T with du/dn tangent.
   *
   * @param d Vertex data; cellGradient, gradient and shearRate are written.
   */
  void recoverGradient(const Nodal& d)
  {
    const int nodes = d.nodes;
    const auto cv = d.cellVertices;
    const auto x = d.coordinates;
    const auto u = d.velocity;
    const auto cg = d.cellGradient;
    ::Kokkos::parallel_for(
      "AtriumRheology::cellGradient", ::Kokkos::RangePolicy<Exec>(0, d.cells),
      KOKKOS_LAMBDA(const int c) {
        Real X[4][3], U[4][3];
        for (int a = 0; a < 4; ++a)
        {
          const int i = cv(4 * c + a);
          for (int k = 0; k < 3; ++k)
          {
            X[a][k] = x(i + k * nodes);
            U[a][k] = u(i + k * nodes);
          }
        }
        // Columns of J are the edges from vertex 0; rows of its inverse are
        // the gradients of the barycentric coordinates 1, 2, 3.
        Real J[3][3];
        for (int k = 0; k < 3; ++k)
        {
          for (int a = 1; a < 4; ++a)
            J[k][a - 1] = X[a][k] - X[0][k];
        }
        const Real det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
          J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
          J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
        Real g[4][3];
        g[1][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / det;
        g[1][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / det;
        g[1][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / det;
        g[2][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / det;
        g[2][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / det;
        g[2][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / det;
        g[3][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / det;
        g[3][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / det;
        g[3][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / det;
        for (int k = 0; k < 3; ++k)
          g[0][k] = -(g[1][k] + g[2][k] + g[3][k]);
        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            cg(10 * c + 3 * i + j) = U[0][i] * g[0][j] + U[1][i] * g[1][j] +
              U[2][i] * g[2][j] + U[3][i] * g[3][j];
          }
        }
        cg(10 * c + 9) = ::Kokkos::fabs(det) / 6;
      });
    const auto offsets = d.nodeOffsets;
    const auto incident = d.nodeCells;
    const auto gradient = d.gradient;
    const auto shearRate = d.shearRate;
    const auto wallNormal = d.wallNormal;
    ::Kokkos::parallel_for(
      "AtriumRheology::vertexGradient", ::Kokkos::RangePolicy<Exec>(0, nodes),
      KOKKOS_LAMBDA(const int i) {
        Real L[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        Real weight = 0;
        for (int k = offsets(i); k < offsets(i + 1); ++k)
        {
          const int c = incident(k);
          const Real volume = cg(10 * c + 9);
          for (int m = 0; m < 9; ++m)
            L[m] += volume * cg(10 * c + m);
          weight += volume;
        }
        for (int m = 0; m < 9; ++m)
          L[m] /= weight;
        const Real n[3] = {
          wallNormal(i), wallNormal(i + nodes), wallNormal(i + 2 * nodes)};
        if (n[0] * n[0] + n[1] * n[1] + n[2] * n[2] > 0.5)
        {
          // No-slip wall: grad u = (du/dn) n^T with du/dn tangent to the wall.
          Real a[3];
          Real an = 0;
          for (int r = 0; r < 3; ++r)
          {
            a[r] = L[3 * r] * n[0] + L[3 * r + 1] * n[1] + L[3 * r + 2] * n[2];
            an += a[r] * n[r];
          }
          for (int r = 0; r < 3; ++r)
          {
            for (int c = 0; c < 3; ++c)
              L[3 * r + c] = (a[r] - an * n[r]) * n[c];
          }
        }
        else
        {
          // div u = 0.
          const Real third = (L[0] + L[4] + L[8]) / 3;
          L[0] -= third;
          L[4] -= third;
          L[8] -= third;
        }
        for (int m = 0; m < 9; ++m)
          gradient(9 * i + m) = L[m];
        Real s2 = 0;
        for (int a = 0; a < 3; ++a)
        {
          for (int b = 0; b < 3; ++b)
          {
            const Real e = 0.5 * (L[3 * a + b] + L[3 * b + a]);
            s2 += e * e;
          }
        }
        shearRate(i) = ::Kokkos::sqrt(2 * s2);
      });
  }

  /**
   * @brief Exact solution of dn/dt = -(b/2)(n - n_st)(n + n_st - 1) with b and
   *        n_st frozen: (n - n_st)/(n + n_st - 1) decays as
   *        exp(-(b/2)(2 n_st - 1) t).
   * @param n Aggregate size at the start, >= 1. @param nst Steady size, >= 1.
   * @param b Fragmentation rate. @param t Time. @return Aggregate size at t.
   */
  KOKKOS_INLINE_FUNCTION Real aggregateSize(Real n, Real nst, Real b, Real t)
  {
    const Real q =
      (n - nst) / (n + nst - 1) * ::Kokkos::exp(-0.5 * b * (2 * nst - 1) * t);
    return (nst + q * (nst - 1)) / (1 - q);
  }

  /**
   * @brief Aggregate size and psi of every vertex over dt from the advected
   *        state, with the velocity gradient of recoverGradient() frozen.
   * @param d Vertex data; psiAdvected and sizeAdvected are read, relaxation
   *        is written.
   * @param psiTo Receives psi. @param sizeTo Receives the aggregate size.
   * @param model Fluid constants.
   * @param thixotropic Whether n is transported (TVE) or n = n_st (VE).
   * @param dt Time step.
   * @param epsilon Linear PTT factor of advance().
   */
  void advanceStructure(const Nodal& d, const RealView& psiTo, const RealView& sizeTo,
    const Rheology& model, bool thixotropic, Real dt, Real epsilon)
  {
    const int nodes = d.nodes;
    const auto gradient = d.gradient;
    const auto shearRate = d.shearRate;
    const auto relaxation = d.relaxation;
    const auto psiFrom = d.psiAdvected;
    const auto sizeFrom = d.sizeAdvected;
    ::Kokkos::parallel_for(
      "AtriumRheology::advanceStructure", ::Kokkos::RangePolicy<Exec>(0, nodes),
      KOKKOS_LAMBDA(const int i) {
        const Real gd = shearRate(i);
        Real lambda;
        if (thixotropic)
        {
          const Real nst = model.steadySize(gd);
          const Real b = model.fragmentation(gd);
          const Real n0 = sizeFrom(i) < 1 ? Real(1) : sizeFrom(i);
          lambda = model.relaxation(aggregateSize(n0, nst, b, 0.5 * dt), gd);
          const Real n = aggregateSize(n0, nst, b, dt);
          sizeTo(i) = n < 1 ? Real(1) : n;
        }
        else
        {
          sizeTo(i) = model.steadySize(gd);
          lambda = model.steadyRelaxation(gd);
        }
        relaxation(i) = lambda;
        Tensor L;
        for (int a = 0; a < 3; ++a)
        {
          for (int b = 0; b < 3; ++b)
            L.a[a][b] = gradient(9 * i + 3 * a + b);
        }
        Real s[6];
        for (int k = 0; k < 6; ++k)
          s[k] = psiFrom(i + k * nodes);
        advance(s, L, lambda, dt, epsilon);
        for (int k = 0; k < 6; ++k)
          psiTo(i + k * nodes) = s[k];
      });
  }

  /**
   * @brief Values taken by the characteristics that come from beyond an
   *        opening: developed blood on the veins.
   *
   * At a vertex of a vein with unit normal n, the shear of the axial velocity
   * across the section, S = n g^T with g = (I - n n^T) L^T n, has the shear
   * rate |g|, and for it the steady conformation of the model is
   * c = I + lambda (S + S^T) + 2 lambda^2 S S^T, with lambda and the aggregate
   * size n_st at that shear rate. Relaxed blood entering at the ring where a
   * vein meets the wall is sheared at once by the wall, a stress that blood
   * coming from a vein has built over centimetres. At the other vertices the
   * values are the current ones, so that a backflow on the mitral plane
   * carries the psi found there (as the CPU driver does).
   *
   * @param d Vertex data, with the gradient of recoverGradient(); psiIn and
   *        sizeIn are written.
   * @param model Fluid constants.
   * @param thixotropic Whether lambda is that of n_st (TVE) or of eta_CY (VE).
   * @param developed Whether the veins carry developed blood; otherwise
   *        psi = 0 and the aggregate size @p relaxedSize.
   * @param relaxedSize Aggregate size on the veins without @p developed.
   */
  void inflowValues(const Nodal& d, const Rheology& model, bool thixotropic,
    bool developed, Real relaxedSize)
  {
    const int nodes = d.nodes;
    const auto gradient = d.gradient;
    const auto normal = d.veinNormal;
    const auto psi = d.psi;
    const auto size = d.size;
    const auto psiIn = d.psiIn;
    const auto sizeIn = d.sizeIn;
    ::Kokkos::parallel_for(
      "AtriumRheology::inflowValues", ::Kokkos::RangePolicy<Exec>(0, nodes),
      KOKKOS_LAMBDA(const int i) {
        const Real n[3] = {normal(i), normal(i + nodes), normal(i + 2 * nodes)};
        if (n[0] * n[0] + n[1] * n[1] + n[2] * n[2] < 0.5)
        {
          for (int k = 0; k < 6; ++k)
            psiIn(i + k * nodes) = psi(i + k * nodes);
          sizeIn(i) = size(i);
          return;
        }
        if (!developed)
        {
          for (int k = 0; k < 6; ++k)
            psiIn(i + k * nodes) = 0;
          sizeIn(i) = thixotropic ? relaxedSize : size(i);
          return;
        }
        Real q[3], qn = 0;
        for (int j = 0; j < 3; ++j)
        {
          q[j] = n[0] * gradient(9 * i + j) + n[1] * gradient(9 * i + 3 + j) +
            n[2] * gradient(9 * i + 6 + j);
          qn += q[j] * n[j];
        }
        Real g[3], gd2 = 0;
        for (int j = 0; j < 3; ++j)
        {
          g[j] = q[j] - qn * n[j];
          gd2 += g[j] * g[j];
        }
        const Real gd = ::Kokkos::sqrt(gd2);
        const Real nst = model.steadySize(gd);
        const Real lambda =
          thixotropic ? model.relaxation(nst, gd) : model.steadyRelaxation(gd);
        Tensor c;
        for (int a = 0; a < 3; ++a)
        {
          for (int b = 0; b < 3; ++b)
          {
            c.a[a][b] = (a == b ? Real(1) : Real(0)) +
              lambda * (n[a] * g[b] + g[a] * n[b]) +
              2 * lambda * lambda * gd2 * n[a] * n[b];
          }
        }
        Real w[3];
        Tensor r;
        eigen(c, w, r);
        for (int k = 0; k < 3; ++k)
          w[k] = ::Kokkos::log(w[k] > 1e-300 ? w[k] : Real(1e-300));
        const Tensor l = spectral(w, r);
        psiIn(i) = l.a[0][0];
        psiIn(i + nodes) = l.a[1][1];
        psiIn(i + 2 * nodes) = l.a[2][2];
        psiIn(i + 3 * nodes) = l.a[0][1];
        psiIn(i + 4 * nodes) = l.a[0][2];
        psiIn(i + 5 * nodes) = l.a[1][2];
        sizeIn(i) = thixotropic ? nst : size(i);
      });
  }

  /**
   * @brief Step 2: psi and n at the feet of the characteristics.
   *
   * For every vertex the foot X = x - dt u, u the velocity of Nodal::velocity
   * at the vertex, is located by walking from a cell of the vertex across the
   * face of the most negative barycentric coordinate. Inside the mesh the
   * values of the four vertices of the cell are interpolated; when the walk
   * meets the boundary the foot is brought back to that face, and the values
   * are those of the inflow (an opening) or the current ones (the wall).
   *
   * @param d Vertex data; psi, size, psiIn, sizeIn are read, psiAdvected and
   *        sizeAdvected are written.
   * @param dt Time step.
   */
  void advect(const Nodal& d, Real dt)
  {
    const int nodes = d.nodes;
    const auto cv = d.cellVertices;
    const auto neighbours = d.cellNeighbours;
    const auto x = d.coordinates;
    const auto offsets = d.nodeOffsets;
    const auto incident = d.nodeCells;
    const auto u = d.velocity;
    const auto psi = d.psi;
    const auto size = d.size;
    const auto psiIn = d.psiIn;
    const auto sizeIn = d.sizeIn;
    const auto psiTo = d.psiAdvected;
    const auto sizeTo = d.sizeAdvected;
    ::Kokkos::parallel_for(
      "AtriumRheology::advect", ::Kokkos::RangePolicy<Exec>(0, nodes),
      KOKKOS_LAMBDA(const int i) {
        Real X[3], move2 = 0;
        for (int k = 0; k < 3; ++k)
        {
          const Real step = dt * u(i + k * nodes);
          X[k] = x(i + k * nodes) - step;
          move2 += step * step;
        }
        if (move2 == 0)
        {
          for (int k = 0; k < 6; ++k)
            psiTo(i + k * nodes) = psi(i + k * nodes);
          sizeTo(i) = size(i);
          return;
        }
        int c = incident(offsets(i));
        Real b[4];
        bool opening = false;
        for (int walk = 0; walk < 512; ++walk)
        {
          // Barycentric coordinates of X in the cell c.
          Real P[4][3];
          for (int a = 0; a < 4; ++a)
          {
            for (int k = 0; k < 3; ++k)
              P[a][k] = x(cv(4 * c + a) + k * nodes);
          }
          Real J[3][3], r[3];
          for (int k = 0; k < 3; ++k)
          {
            for (int a = 1; a < 4; ++a)
              J[k][a - 1] = P[a][k] - P[0][k];
            r[k] = X[k] - P[0][k];
          }
          const Real det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
            J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
            J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
          b[1] = (r[0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                   J[0][1] * (r[1] * J[2][2] - J[1][2] * r[2]) +
                   J[0][2] * (r[1] * J[2][1] - J[1][1] * r[2])) /
            det;
          b[2] = (J[0][0] * (r[1] * J[2][2] - J[1][2] * r[2]) -
                   r[0] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                   J[0][2] * (J[1][0] * r[2] - r[1] * J[2][0])) /
            det;
          b[3] = (J[0][0] * (J[1][1] * r[2] - r[1] * J[2][1]) -
                   J[0][1] * (J[1][0] * r[2] - r[1] * J[2][0]) +
                   r[0] * (J[1][0] * J[2][1] - J[1][1] * J[2][0])) /
            det;
          b[0] = 1 - b[1] - b[2] - b[3];
          int out = 0;
          for (int a = 1; a < 4; ++a)
          {
            if (b[a] < b[out])
              out = a;
          }
          if (b[out] >= -1e-12)
            break;
          const int next = neighbours(4 * c + out);
          if (next >= 0 && walk < 511)
          {
            c = next;
            continue;
          }
          // Beyond the boundary (or a walk that does not end): back to the cell.
          opening = next == -2;
          Real sum = 0;
          for (int a = 0; a < 4; ++a)
          {
            if (b[a] < 0)
              b[a] = 0;
            sum += b[a];
          }
          for (int a = 0; a < 4; ++a)
            b[a] /= sum;
          break;
        }
        int v[4];
        for (int a = 0; a < 4; ++a)
          v[a] = cv(4 * c + a);
        for (int k = 0; k < 6; ++k)
        {
          Real value = 0;
          for (int a = 0; a < 4; ++a)
          {
            value += b[a] * (opening ? psiIn(v[a] + k * nodes) : psi(v[a] + k * nodes));
          }
          psiTo(i + k * nodes) = value;
        }
        Real value = 0;
        for (int a = 0; a < 4; ++a)
          value += b[a] * (opening ? sizeIn(v[a]) : size(v[a]));
        sizeTo(i) = value;
      });
  }

  /**
   * @brief Viscosity and explicit stress E of the momentum equation, from
   *        the gradient of recoverGradient() and psi.
   *
   * N: the viscosity mu. CY: eta_CY(gd) and E = eta (L^T - tr(L) I). VE and
   * TVE: E = G (c - I) - eta_a L and the viscosity eta_plasma + eta_a + eta_d,
   * with eta_a = alpha G lambda and eta_d = beta G dt_e lambda_max(c),
   * dt_e = lambda (1 - exp(-dt/lambda)).
   *
   * @param d Vertex data; viscosity, increment, extra0..2 are written.
   * @param psi psi of the stress (elastic fluids).
   * @param model Fluid constants.
   * @param fluid Fluid model.
   * @param mu Viscosity of the Newtonian fluid.
   * @param dt Time step.
   * @param alpha Factor of eta_a.
   * @param beta Factor of eta_d.
   */
  void momentumStress(const Nodal& d, const RealView& psi, const Rheology& model,
    Fluid fluid, Real mu, Real dt, Real alpha, Real beta)
  {
    const int nodes = d.nodes;
    const auto gradient = d.gradient;
    const auto shearRate = d.shearRate;
    const auto viscosity = d.viscosity;
    const auto increment = d.increment;
    const auto relaxation = d.relaxation;
    const auto e0 = d.extra0;
    const auto e1 = d.extra1;
    const auto e2 = d.extra2;
    const bool elastic = fluid == Fluid::Viscoelastic || fluid == Fluid::Thixotropic;
    const bool shearThinning = fluid == Fluid::CarreauYasuda;
    ::Kokkos::parallel_for(
      "AtriumRheology::momentumStress", ::Kokkos::RangePolicy<Exec>(0, nodes),
      KOKKOS_LAMBDA(const int i) {
        Real E[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        Real eta = mu, etaD = 0;
        if (shearThinning)
        {
          Real L[3][3];
          for (int a = 0; a < 3; ++a)
          {
            for (int b = 0; b < 3; ++b)
              L[a][b] = gradient(9 * i + 3 * a + b);
          }
          eta = model.viscosity(shearRate(i));
          const Real divergence = L[0][0] + L[1][1] + L[2][2];
          for (int a = 0; a < 3; ++a)
          {
            for (int b = 0; b < 3; ++b)
              E[a][b] = eta * (L[b][a] - (a == b ? divergence : Real(0)));
          }
        }
        else if (elastic)
        {
          Real s[6];
          for (int k = 0; k < 6; ++k)
            s[k] = psi(i + k * nodes);
          Real largest;
          const Tensor c = conformation(s, largest);
          const Real lambda = relaxation(i);
          const Real etaA = alpha * model.modulus * lambda;
          etaD =
            beta * model.modulus * lambda * (1 - ::Kokkos::exp(-dt / lambda)) * largest;
          eta = model.etaPlasma + etaA + etaD;
          for (int a = 0; a < 3; ++a)
          {
            for (int b = 0; b < 3; ++b)
            {
              E[a][b] = model.modulus * (c.a[a][b] - (a == b ? Real(1) : Real(0))) -
                etaA * gradient(9 * i + 3 * a + b);
            }
          }
        }
        viscosity(i) = eta;
        increment(i) = etaD;
        for (int b = 0; b < 3; ++b)
        {
          e0(i + b * nodes) = E[0][b];
          e1(i + b * nodes) = E[1][b];
          e2(i + b * nodes) = E[2][b];
        }
      });
  }

  /**
   * @brief Volume integrals and maxima of the vertex fields, with the
   *        gradient of recoverGradient().
   * @param d Vertex data; trace, stress are written.
   * @param model Fluid constants.
   * @param elastic Whether the fluid has an elastic stress (psi, relaxation).
   * @return Volume integrals and maxima.
   */
  Diagnostics diagnose(const Nodal& d, const Rheology& model, bool elastic)
  {
    const int nodes = d.nodes;
    const auto gradient = d.gradient;
    const auto shearRate = d.shearRate;
    const auto viscosity = d.viscosity;
    const auto relaxation = d.relaxation;
    const auto size = d.size;
    const auto psi = d.psi;
    const auto trace = d.trace;
    const auto stress = d.stress;
    const auto density = d.density;
    const auto volume = d.volume;
    ::Kokkos::parallel_for(
      "AtriumRheology::diagnose", ::Kokkos::RangePolicy<Exec>(0, nodes),
      KOKKOS_LAMBDA(const int i) {
        const Real gd = shearRate(i);
        Real power = 0, energy = 0, tr = 3, norm = 0;
        if (elastic)
        {
          Real s[6];
          for (int k = 0; k < 6; ++k)
            s[k] = psi(i + k * nodes);
          const Tensor c = conformation(s);
          for (int a = 0; a < 3; ++a)
          {
            for (int b = 0; b < 3; ++b)
            {
              const Real t = model.modulus * (c.a[a][b] - (a == b ? Real(1) : Real(0)));
              power +=
                t * 0.5 * (gradient(9 * i + 3 * a + b) + gradient(9 * i + 3 * b + a));
              norm += t * t;
            }
          }
          tr = c.a[0][0] + c.a[1][1] + c.a[2][2];
          energy = 0.5 * model.modulus * (tr - 3 - (s[0] + s[1] + s[2]));
        }
        trace(i) = tr;
        stress(i) = ::Kokkos::sqrt(norm);
        // 2 eta D:D, eta the viscosity of the solvent for the elastic fluids.
        density(4 * i) = (elastic ? model.etaPlasma : viscosity(i)) * gd * gd;
        density(4 * i + 1) = power;
        density(4 * i + 2) = energy;
        density(4 * i + 3) = size(i);
      });
    Diagnostics res;
    Real sums[4] = {0, 0, 0, 0};
    for (int k = 0; k < 4; ++k)
    {
      ::Kokkos::parallel_reduce(
        "AtriumRheology::integrate", ::Kokkos::RangePolicy<Exec>(0, nodes),
        KOKKOS_LAMBDA(const int i, Real& sum) { sum += volume(i) * density(4 * i + k); },
        sums[k]);
    }
    Real total = 0;
    ::Kokkos::parallel_reduce(
      "AtriumRheology::volume", ::Kokkos::RangePolicy<Exec>(0, nodes),
      KOKKOS_LAMBDA(const int i, Real& sum) { sum += volume(i); }, total);
    res.solventDissipation = sums[0];
    res.polymerPower = sums[1];
    res.elasticEnergy = sums[2];
    res.meanSize = sums[3] / total;
    if (elastic)
    {
      ::Kokkos::parallel_reduce(
        "AtriumRheology::maxWeissenberg", ::Kokkos::RangePolicy<Exec>(0, nodes),
        KOKKOS_LAMBDA(const int i, Real& m) {
          const Real wi = relaxation(i) * shearRate(i);
          if (wi > m)
            m = wi;
        },
        ::Kokkos::Max<Real>(res.maxWeissenberg));
      ::Kokkos::parallel_reduce(
        "AtriumRheology::maxTrace", ::Kokkos::RangePolicy<Exec>(0, nodes),
        KOKKOS_LAMBDA(const int i, Real& m) {
          if (trace(i) > m)
            m = trace(i);
        },
        ::Kokkos::Max<Real>(res.maxTrace));
    }
    return res;
  }

  /// @brief Smooth step on [0, 1].
  Real smooth(Real s)
  {
    s = s < 0 ? 0 : (s > 1 ? 1 : s);
    return s * s * (3 - 2 * s);
  }

  /// @brief Ventricular activation of Atrium/Fluid, periodic.
  Real activation(Real period, Real t)
  {
    constexpr Real rampStart = 0.15, rampEnd = 0.21, plateauEnd = 0.36, relaxEnd = 0.45,
                   negativeEnd = 0.6, positive = 35, negative = -20;
    const Real tau = t - period * std::floor(t / period);
    if (tau < rampStart)
      return 0;
    if (tau < rampEnd)
      return positive * smooth((tau - rampStart) / (rampEnd - rampStart));
    if (tau < plateauEnd)
      return positive;
    if (tau < relaxEnd)
    {
      return positive +
        (negative - positive) * smooth((tau - plateauEnd) / (relaxEnd - plateauEnd));
    }
    if (tau < negativeEnd)
      return negative;
    return negative * (1 - smooth((tau - negativeEnd) / (period - negativeEnd)));
  }

  /**
   * @brief Lumped ventricle of Atrium/Fluid, its parameters unchanged.
   * @param period Cardiac period.
   * @param atrialPressure Mean pressure of the mitral plane, read at every step.
   */
  Ventricle::Input ventricleInput(Real period, const Real* atrialPressure)
  {
    Ventricle::Input input;
    input.rho = 1.0e3;
    input.R0 = 2.45e-2;
    input.d0 = 1.4e-2;
    input.Es = 3.0e6;
    input.mu = 70.0;
    input.eta = 70.0;
    input.alpha = 1.5;
    input.alphaR = 0.12;
    input.k0 = 1.0e5;
    input.sigma0 = 1.2e5;
    input.Rp = 5.0e7;
    input.Cp = 5.0e-9;
    input.Rd = 1.0e8;
    input.Cd = 5.0e-10;
    input.proximalRadius = 0.0125;
    input.proximalLength = 0.4;
    input.distalRadius = 0.00175;
    input.distalLength = 0.2;
    input.mu_0 = 0.301;
    input.mu_Inf = 0.0055;
    input.lambda = 16.152;
    input.n = 0.21;
    input.yasuda = 0.77;
    input.Kat = 1.0e-6;
    input.Kp = 5.0e-10;
    input.Kar = 2.0e-7;
    input.cavityCapacity = 5.0e-12;
    input.localTolerance = 1.0e-12;
    input.localMaxIterations = 50;
    input.localDamping = 1.0;
    input.absRegularization = 1.0e-14;
    input.pSv = [](Real) { return Real(1.0e3); };
    input.pAt = [atrialPressure](Real) { return *atrialPressure; };
    input.u = [period](Real t) { return activation(period, t); };
    input.m0 = [](Real ec) {
      constexpr Real low = 1.6, high = 1.0, lowEc = 0.0, highEc = 2.0;
      if (ec >= highEc)
        return high;
      if (ec <= lowEc)
        return low;
      const Real s = (ec - lowEc) / (highEc - lowEc);
      return (1 - s) * low + s * high;
    };
    input.dm0 = [](Real ec) {
      constexpr Real low = 1.6, high = 1.0, lowEc = 0.0, highEc = 2.0;
      if (ec <= lowEc || ec >= highEc)
        return Real(0);
      return (high - low) / (highEc - lowEc);
    };
    using PassiveEnergy = std::decay_t<decltype(input.passiveEnergy)>;
    typename PassiveEnergy::Parameters hp;
    hp.C0 = 1.9e3;
    hp.C1 = 1.1e-1;
    hp.C2 = 1.9e3;
    hp.C3 = 1.1e-1;
    input.passiveEnergy = PassiveEnergy(hp);
    return input;
  }

  /// @brief Reads key=value arguments into the parameters.
  Parameters parse(int argc, char** argv)
  {
    Parameters P;
    std::map<std::string, std::string> kv;
    for (int i = 2; i < argc; ++i)
    {
      const std::string a = argv[i];
      const auto eq = a.find('=');
      if (a.rfind("-", 0) == 0 || eq == std::string::npos)
        continue; // Kokkos and PETSc arguments.
      kv[a.substr(0, eq)] = a.substr(eq + 1);
    }
    const auto real = [&](const char* key, Real& out) {
      if (kv.count(key))
        out = std::stod(kv[key]);
    };
    const auto count = [&](const char* key, size_t& out) {
      if (kv.count(key))
        out = std::stoul(kv[key]);
    };
    const auto choice = [&](const char* key, const char* yes, const char* no, bool& out) {
      if (!kv.count(key))
        return;
      if (kv[key] != yes && kv[key] != no)
      {
        throw std::invalid_argument(std::string(key) + " is " + yes + " or " + no + ".");
      }
      out = kv[key] == yes;
    };
    real("scale", P.scale);
    real("period", P.period);
    count("cycles", P.cycles);
    count("steps", P.steps);
    real("t_end", P.endTime);
    real("rho", P.rho);
    real("mu", P.mu);
    real("p_in", P.inletPressureMean);
    real("p_in_amplitude", P.inletPressureAmplitude);
    real("ramp", P.rampCycles);
    real("z_in", P.inletImpedance);
    real("p_out", P.outletPressure);
    real("n_in", P.inletSize);
    real("n_init", P.initialSize);
    real("devss", P.devss);
    real("implicit", P.implicit);
    real("ptt", P.ptt);
    real("eta_plasma", P.etaPlasma);
    real("lambda_h", P.lambdaH);
    real("tolerance", P.tolerance);
    count("max_iterations", P.maxIterations);
    count("pc_lag", P.preconditionerLag);
    count("output_every", P.outputEvery);
    choice("stab", "asgs", "split", P.residual);
    choice("inlet", "developed", "relaxed", P.developedInlet);
    if (kv.count("lv"))
      P.ventricle = kv["lv"] != "0";
    if (kv.count("output"))
      P.output = kv["output"];
    if (kv.count("restart"))
      P.restart = kv["restart"];
    if (kv.count("order"))
    {
      P.order = std::stoi(kv["order"]);
      if (P.order != 1 && P.order != 2)
        throw std::invalid_argument("order is 1 or 2.");
    }
    if (kv.count("fluid"))
    {
      std::string f = kv["fluid"];
      std::transform(f.begin(), f.end(), f.begin(),
        [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
      if (f == "N")
        P.fluid = Fluid::Newtonian;
      else if (f == "CY")
        P.fluid = Fluid::CarreauYasuda;
      else if (f == "VE")
        P.fluid = Fluid::Viscoelastic;
      else if (f == "TVE")
        P.fluid = Fluid::Thixotropic;
      else
        throw std::invalid_argument("fluid is N, CY, VE or TVE.");
      P.fluidName = f;
    }
    return P;
  }

  /**
   * @brief Copies a sequential PETSc vector into a host vector.
   * @param[in] x PETSc vector.
   * @return Host copy of @p x.
   */
  Rodin::Math::Vector<Real> toHost(::Vec x)
  {
    PetscInt n;
    PetscCallAbort(PETSC_COMM_SELF, VecGetLocalSize(x, &n));
    const PetscScalar* data;
    PetscCallAbort(PETSC_COMM_SELF, VecGetArrayRead(x, &data));
    Rodin::Math::Vector<Real> res = Eigen::Map<const Rodin::Math::Vector<Real>>(data, n);
    PetscCallAbort(PETSC_COMM_SELF, VecRestoreArrayRead(x, &data));
    return res;
  }
}

/**
   * @brief Entry (i, k) of E - eta_d grad u^n at a point: E from the vertex
   *        values of its row i, grad u^n that of the cell.
   */
template <class VectorField, class ScalarField>
struct StressEntry
{
    VectorField row; ///< Row i of E.
    ScalarField etaD; ///< eta_d.
    VectorField old; ///< u^n.
    ::Rodin::Index i, k;

    KOKKOS_INLINE_FUNCTION Real operator()(
      const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& point) const
    {
      return row.getValue(point)(k) -
        etaD.getValue(point) * old.u.interpolateJacobian(point)(i, k);
    }
};

/// @brief Writes and reads the checkpoint of a run (see the file header).
struct Checkpoint
{
    static constexpr char magic[8] = {'A', 'T', 'R', 'C', 'K', 'P', 'T', '1'};

    template <class T>
    static void put(std::ofstream& f, const T& x)
    {
      static_assert(std::is_trivially_copyable_v<T>);
      f.write(reinterpret_cast<const char*>(&x), sizeof(T));
    }

    template <class T>
    static void get(std::ifstream& f, T& x)
    {
      static_assert(std::is_trivially_copyable_v<T>);
      f.read(reinterpret_cast<char*>(&x), sizeof(T));
    }

    static void putVector(std::ofstream& f, const Real* x, std::uint64_t n)
    {
      put(f, n);
      f.write(
        reinterpret_cast<const char*>(x), static_cast<std::streamsize>(n * sizeof(Real)));
    }

    static bool getVector(std::ifstream& f, Real* x, std::uint64_t n)
    {
      std::uint64_t m = 0;
      get(f, m);
      if (m != n)
        return false;
      f.read(reinterpret_cast<char*>(x), static_cast<std::streamsize>(n * sizeof(Real)));
      return static_cast<bool>(f);
    }
};

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;
  using PointView = ::Rodin::Kokkos::View<::Rodin::Geometry::Point>;

  if (argc < 2)
  {
    std::cerr << "Usage: KokkosAtriumRheologyPETSc mesh.mesh [key=value ...]"
              << std::endl;
    return EXIT_FAILURE;
  }
  ::Kokkos::ScopeGuard runtime(argc, argv);
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);
  {
    const Parameters P = parse(argc, argv);
    const Fluid fluid = P.fluid;
    const bool elastic = fluid == Fluid::Viscoelastic || fluid == Fluid::Thixotropic;
    const bool thixotropic = fluid == Fluid::Thixotropic;

    Rheology model;
    model.lambdaH = P.lambdaH;
    model.etaPlasma = P.etaPlasma;
    model.etaPolymer = model.etaInf - model.etaPlasma;
    model.modulus = model.etaPolymer / model.lambdaH;
    // b(gd) falls to zero with a(gd) near gd_max, which would freeze n in the
    // jets: b keeps its maximum above the shear rate where it peaks.
    {
      Real peak = 0, at = 1e300;
      for (Real gd = 1e-2; gd <= model.gdMax; gd *= 1.02)
      {
        const Real b = model.fragmentationRaw(gd);
        if (b > peak)
        {
          peak = b;
          at = gd;
        }
      }
      model.bPeak = peak;
      model.gdPeak = at;
    }
    const Real inletSize = P.inletSize > 0 ? P.inletSize : model.steadySize(0);
    const Real initialSize = P.initialSize < 0
      ? model.steadySize(0)
      : (P.initialSize > 0 ? P.initialSize : inletSize);
    std::cout << "[fluid] " << P.fluidName << "  stabilisation "
              << (P.residual ? "asgs" : "split");
    if (fluid == Fluid::Newtonian)
      std::cout << "  mu = " << P.mu << " Pa s";
    else
    {
      std::cout << "  eta_CY(0) = " << model.viscosity(0)
                << "  eta_CY(inf) = " << model.etaInf << " Pa s";
    }
    if (elastic)
    {
      std::cout << "  G = " << model.modulus
                << " Pa  lambda(0) = " << model.steadyRelaxation(0)
                << " s  lambda_H = " << model.lambdaH
                << " s  n_st(0) = " << model.steadySize(0);
    }
    if (thixotropic)
    {
      std::cout << "  n_in = " << inletSize << "  n_init = " << initialSize
                << "  b floor = " << model.bPeak << " 1/s above " << model.gdPeak
                << " 1/s";
    }
    std::cout << std::endl;

    // Labels of the Atrium meshes.
    constexpr Attribute wall = 4;
    constexpr Attribute mitral = 6;
    constexpr std::array<Attribute, 4> veins{1, 2, 3, 5};

    Mesh host;
    host.load(argv[1], IO::FileFormat::MEDIT);
    if (host.getDimension() != 3 || host.getCellCount() == 0 ||
      host.getPolytopeCount(Polytope::Type::Tetrahedron) != host.getCellCount())
    {
      std::cerr << argv[1]
                << " is not a tetrahedral mesh of the fluid (a Git LFS pointer, "
                << "or a combined fluid-solid mesh?)." << std::endl;
      return EXIT_FAILURE;
    }
    host.scale(P.scale);
    // Face to cell incidence for the boundary; the remaining incidences for the
    // host space of the XDMF output.
    for (size_t d = 1; d <= 3; ++d)
    {
      for (size_t e = 0; e < d; ++e)
      {
        host.getConnectivity().compute(d, e);
        host.getConnectivity().compute(e, d);
      }
    }
    const Mesh<DeviceContext> mesh(DeviceContext{}, host);

    H1 vh(std::integral_constant<size_t, 1>{}, mesh, 3);
    H1 qh(std::integral_constant<size_t, 1>{}, mesh);
    PETSc::Variational::TrialFunction u(vh);
    PETSc::Variational::TrialFunction p(qh);
    PETSc::Variational::TestFunction v(vh);
    PETSc::Variational::TestFunction q(qh);
    u.setName("u");
    p.setName("p");
    // Device test functions of the diagnostic forms, kept on the host.
    TestFunction vk(vh);
    TestFunction qk(qh);

    const size_t nodes = host.getVertexCount();
    const size_t cells = host.getCellCount();
    if (vh.getSize() != 3 * nodes || qh.getSize() != nodes)
    {
      std::cerr << "The P1 spaces are not numbered by vertex." << std::endl;
      return EXIT_FAILURE;
    }
    using HostVector = Math::Vector<Real>;
    GridFunction velocity(vh);
    GridFunction pressure(qh);
    GridFunction viscosity(qh);
    GridFunction shearRate(qh);
    GridFunction relaxation(qh);
    GridFunction aggregate(qh);
    GridFunction trace(qh);
    GridFunction stress(qh);
    GridFunction increment(qh);
    GridFunction extra0(vh);
    GridFunction extra1(vh);
    GridFunction extra2(vh);
    GridFunction psiDiagonal(vh);
    GridFunction psiOffDiagonal(vh);
    // Convecting velocity u* and history of the time derivative.
    GridFunction advecting(vh);
    GridFunction past(vh);
    velocity = Zero(3);
    pressure = 0.0;
    increment.setData(HostVector(HostVector::Zero(nodes)));
    viscosity.setData(HostVector(HostVector::Constant(
      nodes, fluid == Fluid::Newtonian ? P.mu : model.viscosity(0))));
    shearRate.setData(HostVector(HostVector::Zero(nodes)));
    relaxation.setData(
      HostVector(HostVector::Constant(nodes, model.steadyRelaxation(0))));
    aggregate.setData(HostVector(HostVector::Constant(nodes, initialSize)));
    trace.setData(HostVector(HostVector::Constant(nodes, 3)));
    stress.setData(HostVector(HostVector::Zero(nodes)));
    extra0.setData(HostVector(HostVector::Zero(3 * nodes)));
    extra1.setData(HostVector(HostVector::Zero(3 * nodes)));
    extra2.setData(HostVector(HostVector::Zero(3 * nodes)));
    psiDiagonal.setData(HostVector(HostVector::Zero(3 * nodes)));
    psiOffDiagonal.setData(HostVector(HostVector::Zero(3 * nodes)));
    std::cout << "[mesh] " << argv[1] << "  cells=" << cells << " vertices=" << nodes
              << " unknowns=" << (vh.getSize() + qh.getSize()) << std::endl;

    // Vertex data of the rheology kernels, indexed as the data of the grid
    // functions: by the degree of freedom of the vertex in the scalar space,
    // which is numbered by cell and is not the index of the vertex.
    Nodal nodal;
    HostVector lumpedVolume = HostVector::Zero(nodes);
    {
      std::vector<int> dof(nodes);
      {
        std::vector<char> taken(nodes, 0);
        for (size_t i = 0; i < nodes; ++i)
        {
          const auto& scalar = qh.getDOFs(0, i);
          const auto& vector = vh.getDOFs(0, i);
          const size_t d = static_cast<size_t>(scalar(0));
          bool valid = scalar.size() == 1 && vector.size() == 3 && d < nodes && !taken[d];
          for (size_t r = 0; valid && r < 3; ++r)
            valid = static_cast<size_t>(vector(r)) == d + r * nodes;
          if (!valid)
          {
            std::cerr << "Unexpected numbering of the P1 spaces at vertex " << i
                      << std::endl;
            return EXIT_FAILURE;
          }
          taken[d] = 1;
          dof[i] = static_cast<int>(d);
        }
      }
      std::vector<int> cv(4 * cells);
      std::vector<Real> x(3 * nodes);
      std::vector<int> offsets(nodes + 1, 0);
      for (size_t c = 0; c < cells; ++c)
      {
        const auto& key = host.getConnectivity().getPolytope(3, c);
        for (std::uint8_t a = 0; a < 4; ++a)
        {
          cv[4 * c + a] = dof[key(a)];
          ++offsets[dof[key(a)] + 1];
        }
      }
      for (size_t i = 0; i < nodes; ++i)
      {
        const auto X = host.getVertexCoordinates(i);
        for (size_t k = 0; k < 3; ++k)
          x[dof[i] + k * nodes] = X(k);
        offsets[i + 1] += offsets[i];
      }
      std::vector<int> incident(4 * cells);
      std::vector<int> next(offsets.begin(), offsets.end() - 1);
      for (size_t c = 0; c < cells; ++c)
      {
        Real e[3][3];
        for (int a = 1; a < 4; ++a)
        {
          for (int k = 0; k < 3; ++k)
            e[a - 1][k] = x[cv[4 * c + a] + k * nodes] - x[cv[4 * c] + k * nodes];
        }
        const Real det = e[0][0] * (e[1][1] * e[2][2] - e[1][2] * e[2][1]) -
          e[0][1] * (e[1][0] * e[2][2] - e[1][2] * e[2][0]) +
          e[0][2] * (e[1][0] * e[2][1] - e[1][1] * e[2][0]);
        for (int a = 0; a < 4; ++a)
        {
          incident[next[cv[4 * c + a]]++] = static_cast<int>(c);
          lumpedVolume(cv[4 * c + a]) += std::abs(det) / 24;
        }
      }
      nodal.nodes = static_cast<int>(nodes);
      nodal.cells = static_cast<int>(cells);
      nodal.cellVertices = upload("AtriumRheology::cellVertices", cv.data(), cv.size());
      nodal.coordinates = upload("AtriumRheology::coordinates", x.data(), x.size());
      nodal.nodeOffsets =
        upload("AtriumRheology::nodeOffsets", offsets.data(), offsets.size());
      nodal.nodeCells =
        upload("AtriumRheology::nodeCells", incident.data(), incident.size());
      nodal.volume = upload("AtriumRheology::volume", lumpedVolume.data(), nodes);
      // Cell across each face, for the walk of advect(): faces are matched by
      // their sorted vertices; the unmatched ones are on the boundary.
      if (nodes >= (size_t(1) << 21))
      {
        std::cerr << "More than 2^21 vertices." << std::endl;
        return EXIT_FAILURE;
      }
      const auto faceKey = [](std::uint64_t a, std::uint64_t b, std::uint64_t c) {
        if (a > b)
          std::swap(a, b);
        if (b > c)
          std::swap(b, c);
        if (a > b)
          std::swap(a, b);
        return (a << 42) | (b << 21) | c;
      };
      std::vector<int> neighbours(4 * cells, -1);
      {
        std::unordered_map<std::uint64_t, int> open;
        open.reserve(2 * cells);
        for (size_t c = 0; c < cells; ++c)
        {
          for (int a = 0; a < 4; ++a)
          {
            const std::uint64_t key = faceKey(
              cv[4 * c + (a + 1) % 4], cv[4 * c + (a + 2) % 4], cv[4 * c + (a + 3) % 4]);
            const auto found = open.find(key);
            if (found == open.end())
              open.emplace(key, static_cast<int>(4 * c + a));
            else
            {
              neighbours[4 * c + a] = found->second / 4;
              neighbours[found->second] = static_cast<int>(c);
              open.erase(found);
            }
          }
        }
        for (Index f = 0; f < host.getFaceCount(); ++f)
        {
          if (!host.isBoundary(f))
            continue;
          const auto attribute = host.getAttribute(2, f);
          if (!attribute || *attribute == wall)
            continue;
          const auto& key = host.getConnectivity().getPolytope(2, f);
          const auto found = open.find(faceKey(dof[key(0)], dof[key(1)], dof[key(2)]));
          if (found != open.end())
            neighbours[found->second] = -2;
        }
      }
      nodal.cellNeighbours =
        upload("AtriumRheology::cellNeighbours", neighbours.data(), neighbours.size());
      // Normal of a boundary at its vertices: the dominant direction of the
      // sum of |f| n n^T over its faces at the vertex, which needs no
      // orientation of the faces.
      const auto vertexNormals = [&](auto&& selected, size_t& count) {
        std::vector<Real> spread(6 * nodes, 0);
        for (Index f = 0; f < host.getFaceCount(); ++f)
        {
          if (!host.isBoundary(f))
            continue;
          const auto attribute = host.getAttribute(2, f);
          if (!attribute || !selected(*attribute))
            continue;
          const auto& key = host.getConnectivity().getPolytope(2, f);
          Real a[3], b[3];
          for (int k = 0; k < 3; ++k)
          {
            a[k] = x[dof[key(1)] + k * nodes] - x[dof[key(0)] + k * nodes];
            b[k] = x[dof[key(2)] + k * nodes] - x[dof[key(0)] + k * nodes];
          }
          const Real n[3] = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
          const Real twice = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
          if (!(twice > 0))
            continue;
          for (std::uint8_t c = 0; c < 3; ++c)
          {
            Real* t = &spread[6 * dof[key(c)]];
            t[0] += n[0] * n[0] / twice;
            t[1] += n[1] * n[1] / twice;
            t[2] += n[2] * n[2] / twice;
            t[3] += n[0] * n[1] / twice;
            t[4] += n[0] * n[2] / twice;
            t[5] += n[1] * n[2] / twice;
          }
        }
        std::vector<Real> normals(3 * nodes, 0);
        count = 0;
        for (size_t i = 0; i < nodes; ++i)
        {
          const Real* t = &spread[6 * i];
          if (!(t[0] + t[1] + t[2] > 0))
            continue;
          Tensor m, vectors;
          m.a[0][0] = t[0];
          m.a[1][1] = t[1];
          m.a[2][2] = t[2];
          m.a[0][1] = m.a[1][0] = t[3];
          m.a[0][2] = m.a[2][0] = t[4];
          m.a[1][2] = m.a[2][1] = t[5];
          Real values[3];
          eigen(m, values, vectors);
          const int k = values[0] >= values[1] && values[0] >= values[2]
            ? 0
            : (values[1] >= values[2] ? 1 : 2);
          for (int c = 0; c < 3; ++c)
            normals[i + c * nodes] = vectors.a[c][k];
          ++count;
        }
        return normals;
      };
      size_t wallNodes = 0, veinNodes = 0;
      const auto wallNormal =
        vertexNormals([&](Attribute attribute) { return attribute == wall; }, wallNodes);
      const auto veinNormal = vertexNormals(
        [&](Attribute attribute) {
          return std::find(veins.begin(), veins.end(), attribute) != veins.end();
        },
        veinNodes);
      nodal.wallNormal =
        upload("AtriumRheology::wallNormal", wallNormal.data(), wallNormal.size());
      nodal.veinNormal =
        upload("AtriumRheology::veinNormal", veinNormal.data(), veinNormal.size());
      nodal.psiIn = RealView("AtriumRheology::psiIn", 6 * nodes);
      nodal.sizeIn = RealView("AtriumRheology::sizeIn", nodes);
      std::cout << "[mesh] wall vertices = " << wallNodes
                << "  vein vertices = " << veinNodes << std::endl;
      nodal.cellGradient = RealView("AtriumRheology::cellGradient", 10 * cells);
      nodal.velocity = RealView("AtriumRheology::velocity", 3 * nodes);
      nodal.gradient = RealView("AtriumRheology::gradient", 9 * nodes);
      nodal.shearRate = RealView("AtriumRheology::shearRate", nodes);
      nodal.viscosity = RealView("AtriumRheology::viscosity", nodes);
      nodal.increment = RealView("AtriumRheology::increment", nodes);
      nodal.sizeAdvected = RealView("AtriumRheology::sizeAdvected", nodes);
      nodal.psiAdvected = RealView("AtriumRheology::psiAdvected", 6 * nodes);
      nodal.sizePredicted = RealView("AtriumRheology::sizePredicted", nodes);
      nodal.psiPredicted = RealView("AtriumRheology::psiPredicted", 6 * nodes);
      nodal.relaxation = RealView("AtriumRheology::relaxation", nodes);
      nodal.size = RealView("AtriumRheology::size", nodes);
      nodal.psi = RealView("AtriumRheology::psi", 6 * nodes);
      nodal.trace = RealView("AtriumRheology::trace", nodes);
      nodal.stress = RealView("AtriumRheology::stress", nodes);
      nodal.extra0 = RealView("AtriumRheology::extra0", 3 * nodes);
      nodal.extra1 = RealView("AtriumRheology::extra1", 3 * nodes);
      nodal.extra2 = RealView("AtriumRheology::extra2", 3 * nodes);
      nodal.density = RealView("AtriumRheology::density", 4 * nodes);
      ::Kokkos::deep_copy(nodal.size, initialSize);
      ::Kokkos::deep_copy(nodal.relaxation, model.steadyRelaxation(0));
    }

    // Areas of the openings, and the functionals that give the flow rates and
    // the mean pressure of the mitral plane: int_G phi n and int_G psi, assembled
    // once on the device and kept on the host.
    const BoundaryNormal normal(mesh);
    LinearForm veinFlux(vk);
    LinearForm mitralFlux(vk);
    LinearForm veinMean(qk);
    LinearForm mitralMean(qk);
    veinFlux = BoundaryIntegral(normal, vk).over(veins[0], veins[1], veins[2], veins[3]);
    mitralFlux = BoundaryIntegral(normal, vk).over(mitral);
    veinMean =
      BoundaryIntegral(RealFunction(1), qk).over(veins[0], veins[1], veins[2], veins[3]);
    mitralMean = BoundaryIntegral(RealFunction(1), qk).over(mitral);
    veinFlux.assemble();
    mitralFlux.assemble();
    veinMean.assemble();
    mitralMean.assemble();
    const auto veinFluxHost = veinFlux.getVector().toHost();
    const auto mitralFluxHost = mitralFlux.getVector().toHost();
    const auto veinMeanHost = veinMean.getVector().toHost();
    const auto mitralMeanHost = mitralMean.getVector().toHost();
    const Real veinArea = veinMeanHost.sum();
    const Real mitralArea = mitralMeanHost.sum();
    const Real volume = host.getVolume();
    const Real h = std::cbrt(volume / static_cast<Real>(cells));
    std::cout << "[boundary] venous area = " << veinArea
              << " m^2  mitral area = " << mitralArea << " m^2  volume = " << 1e6 * volume
              << " mL  mean cell size = " << 1e3 * h << " mm" << std::endl;
    if (!(veinArea > 0) || !(mitralArea > 0))
    {
      std::cerr << "No venous or mitral faces: labels 1, 2, 3, 5 and 6 are expected."
                << std::endl;
      return EXIT_FAILURE;
    }

    // Lumped ventricle, on the host. It reads the mean pressure of the mitral
    // plane at the last step as its atrial pressure.
    Real atrialPressure = P.rampCycles > 0 ? Real(0) : P.inletPressureMean;
    Ventricle::Input input = ventricleInput(P.period, &atrialPressure);
    Ventricle ventricle(input);
    {
      Ventricle::State s0;
      s0.t = 0.0;
      s0.y = 0.0;
      s0.v = 0.0;
      s0.pv = atrialPressure - 100.0;
      s0.par = 11000.0;
      s0.pd = 10000.0;
      s0.w = input.m0(s0.ec);
      ventricle.setMaxIterations(200)
        .setAbsoluteTolerance(1.0e-8)
        .setRelativeTolerance(1.0e-8)
        .setStepTolerance(1.0e-10)
        .setDampingFactor(1.0);
      ventricle.initialize(s0);
    }

    // Restart: flow, history, ventricle and, for the same fluid, psi and n.
    const Real dt = P.period / static_cast<Real>(P.steps);
    size_t stepOffset = 0;
    HostVector previousVelocity = HostVector::Zero(3 * nodes); // u^{n-1}
    if (!P.restart.empty())
    {
      std::ifstream f(P.restart, std::ios::binary);
      char magic[8] = {};
      f.read(magic, 8);
      std::uint64_t savedNodes = 0, savedStep = 0;
      std::int32_t savedFluid = -1;
      Real savedDt = 0, savedPressure = 0;
      Checkpoint::get(f, savedNodes);
      Checkpoint::get(f, savedStep);
      Checkpoint::get(f, savedFluid);
      Checkpoint::get(f, savedDt);
      Checkpoint::get(f, savedPressure);
      if (!f || std::string(magic, 8) != std::string(Checkpoint::magic, 8))
      {
        std::cerr << "Cannot read the checkpoint " << P.restart << std::endl;
        return EXIT_FAILURE;
      }
      if (savedNodes != nodes || std::abs(savedDt - dt) > 1e-12 * dt)
      {
        std::cerr << "The checkpoint " << P.restart
                  << " is of another mesh or time step (" << savedNodes
                  << " vertices, dt = " << savedDt << " s)" << std::endl;
        return EXIT_FAILURE;
      }
      HostVector un(3 * nodes), psi(6 * nodes), size(nodes);
      bool ok = Checkpoint::getVector(f, un.data(), 3 * nodes) &&
        Checkpoint::getVector(f, previousVelocity.data(), 3 * nodes) &&
        Checkpoint::getVector(f, psi.data(), 6 * nodes) &&
        Checkpoint::getVector(f, size.data(), nodes);
      Ventricle::State state;
      Ventricle::History history;
      std::uint8_t hasOlder = 0;
      Checkpoint::get(f, state);
      Checkpoint::get(f, history.n);
      Checkpoint::get(f, history.nm1);
      Checkpoint::get(f, hasOlder);
      if (hasOlder)
      {
        Ventricle::State older;
        Checkpoint::get(f, older);
        history.nm2 = older;
      }
      std::uint64_t unknowns = 0;
      Checkpoint::get(f, unknowns);
      Ventricle::DenseVector x(static_cast<Eigen::Index>(unknowns));
      f.read(reinterpret_cast<char*>(x.data()),
        static_cast<std::streamsize>(unknowns * sizeof(Real)));
      if (!ok || !f)
      {
        std::cerr << "The checkpoint " << P.restart << " is incomplete" << std::endl;
        return EXIT_FAILURE;
      }
      velocity.setData(un);
      ventricle.restore(state, history, x, ventricle.getReport());
      atrialPressure = savedPressure;
      stepOffset = savedStep;
      const bool sameFluid = savedFluid == static_cast<std::int32_t>(fluid);
      if (elastic && sameFluid)
      {
        toDevice(psi.data(), nodal.psi);
        toDevice(size.data(), nodal.size);
        aggregate.setData(size);
      }
      std::cout << "[restart] " << P.restart << ": step " << savedStep
                << ", t = " << savedStep * dt << " s"
                << (elastic ? (sameFluid ? ", psi and n from the file"
                                         : ", psi = 0 and n = n_init (other fluid)")
                            : "")
                << std::endl;
    }
    const auto writeCheckpoint = [&](size_t step, const HostVector& un) {
      const std::string name = P.output + ".restart";
      {
        std::ofstream f(name + ".tmp", std::ios::binary);
        f.write(Checkpoint::magic, 8);
        Checkpoint::put(f, static_cast<std::uint64_t>(nodes));
        Checkpoint::put(f, static_cast<std::uint64_t>(step));
        Checkpoint::put(f, static_cast<std::int32_t>(fluid));
        Checkpoint::put(f, dt);
        Checkpoint::put(f, atrialPressure);
        const HostVector psi = toHost(nodal.psi), size = toHost(nodal.size);
        Checkpoint::putVector(f, un.data(), 3 * nodes);
        Checkpoint::putVector(f, previousVelocity.data(), 3 * nodes);
        Checkpoint::putVector(f, psi.data(), 6 * nodes);
        Checkpoint::putVector(f, size.data(), nodes);
        const auto& history = ventricle.getHistory();
        Checkpoint::put(f, ventricle.getState());
        Checkpoint::put(f, history.n);
        Checkpoint::put(f, history.nm1);
        const std::uint8_t hasOlder = history.nm2.has_value() ? 1 : 0;
        Checkpoint::put(f, hasOlder);
        if (hasOlder)
          Checkpoint::put(f, *history.nm2);
        const auto& x = ventricle.getUnknowns();
        Checkpoint::put(f, static_cast<std::uint64_t>(x.size()));
        f.write(reinterpret_cast<const char*>(x.data()),
          static_cast<std::streamsize>(x.size() * sizeof(Real)));
      }
      std::rename((name + ".tmp").c_str(), name.c_str());
    };

    IO::XDMF xdmf(P.output);
    {
      auto grid = xdmf.grid();
      grid.setMesh(mesh).add("velocity", velocity).add("pressure", pressure);
      grid.add("shear_rate", shearRate).add("viscosity", viscosity);
      if (elastic)
      {
        grid.add("relaxation_time", relaxation)
          .add("aggregate_size", aggregate)
          .add("trace_conformation", trace)
          .add("elastic_stress", stress)
          .add("psi_diagonal", psiDiagonal)
          .add("psi_off_diagonal", psiOffDiagonal);
      }
    }
    xdmf.write(stepOffset * dt).flush();
    std::ofstream csv(P.output + ".csv");
    csv << "t,cycle,pVenous,pOut,valveOpen,qIn,qOut,pInletMean,pOutletMean,"
           "maxSpeedComponent,kineticEnergy,timeDissipation,viscousDissipation,"
           "polymerPower,"
           "elasticEnergy,maxWeissenberg,maxTraceConformation,meanAggregateSize,"
           "rheology,assembly,solve,iterations\n";

    // Flow problem and solver; the sparsity pattern is the same at every step.
    Problem flow(u, p, v, q);
    PETSc::Solver::KSP solver(flow);
    solver.setType(KSPBCGSL);
    // No divergence tolerance: PETSc tests the intermediate residuals of the
    // BiCG part of BiCGStab(l), which spike by many orders of magnitude in
    // solves that converge.
    solver.setTolerances(
      P.tolerance, Real(1e-50), Real(1e30), static_cast<PetscInt>(P.maxIterations));
    {
      ::PC pc;
      // BiCGStab(4) unless the options choose a method or its degree.
      char type[64] = "";
      PetscBool typed = PETSC_FALSE, given = PETSC_FALSE;
      PetscCallAbort(PETSC_COMM_SELF,
        PetscOptionsGetString(
          PETSC_NULLPTR, PETSC_NULLPTR, "-ksp_type", type, sizeof(type), &typed));
      PetscCallAbort(PETSC_COMM_SELF,
        PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, "-ksp_bcgsl_ell", &given));
      if (!given && (!typed || std::string(type) == KSPBCGSL))
      {
        PetscCallAbort(
          PETSC_COMM_SELF, PetscOptionsSetValue(PETSC_NULLPTR, "-ksp_bcgsl_ell", "4"));
      }
      PetscCallAbort(PETSC_COMM_SELF, KSPGetPC(solver.getHandle(), &pc));
      PetscCallAbort(PETSC_COMM_SELF, PCSetType(pc, PCJACOBI));
      // Each step starts from the solution of the previous one.
      PetscCallAbort(
        PETSC_COMM_SELF, KSPSetInitialGuessNonzero(solver.getHandle(), PETSC_TRUE));
    }
    bool fieldSplits = false;
    size_t fallbacks = 0;

    const Real rho = P.rho;
    const Real pi = Math::Constants::pi();
    size_t total = P.cycles * P.steps;
    if (P.endTime > 0)
      total = std::min<size_t>(total, static_cast<size_t>(std::llround(P.endTime / dt)));
    const bool second = P.order == 2;
    const Real dtEff = second ? 2 * dt / 3 : dt; // BDF2: (3/2)(u - past)/dt
    const VectorFunction axis0{Real(1), Real(0), Real(0)};
    const VectorFunction axis1{Real(0), Real(1), Real(0)};
    const VectorFunction axis2{Real(0), Real(0), Real(1)};
    const size_t last = stepOffset + total;
    std::chrono::duration<double> rheologyTime{0}, assembly{0}, solve{0};
    bool valveWasOpen = true;
    for (size_t step = stepOffset + 1; step <= last; ++step)
    {
      const Real time = step * dt;
      const Real rampTime = P.rampCycles * P.period;
      const Real ramp = rampTime > 0 && time < rampTime
        ? 0.5 * (1 - std::cos(pi * time / rampTime))
        : Real(1);
      const Real pIn = ramp *
        (P.inletPressureMean +
          P.inletPressureAmplitude * std::sin(2 * pi * time / P.period));

      // Mitral plane: ventricular pressure and valve law, or a prescribed pressure.
      Real pOut = ramp * P.outletPressure;
      Real mitralImpedance = 0;
      bool valveOpen = true;
      if (P.ventricle)
      {
        if (!ventricle.step(dt).converged)
        {
          std::cerr << "The lumped ventricle did not converge at t = " << time
                    << std::endl;
          return EXIT_FAILURE;
        }
        pOut = ventricle.getState().pv;
        valveOpen = atrialPressure > pOut;
        mitralImpedance = mitralArea / (valveOpen ? input.Kat : input.Kp);
      }

      // Time levels of the velocity: u* = 2 u^n - u^{n-1} convects,
      // u^{n+1/2} = (3 u^n - u^{n-1})/2 carries psi and n, (4 u^n - u^{n-1})/3
      // is the history of BDF2; all u^n for order 1.
      const HostVector un = velocity.getData().toHost();
      const HostVector uStar = second ? HostVector(2 * un - previousVelocity) : un;
      advecting.setData(uStar);
      past.setData(second ? HostVector((4 * un - previousVelocity) / 3) : un);

      const auto r0 = std::chrono::steady_clock::now();
      Diagnostics diagnostics;
      if (elastic)
      {
        // Steps 1 to 3: gradient of u^n, advection, prediction.
        toDevice(un.data(), nodal.velocity);
        recoverGradient(nodal);
        inflowValues(nodal, model, thixotropic, P.developedInlet, inletSize);
        if (second)
        {
          const HostVector uHalf = 1.5 * un - 0.5 * previousVelocity;
          toDevice(uHalf.data(), nodal.velocity);
        }
        advect(nodal, dt);
        advanceStructure(
          nodal, nodal.psiPredicted, nodal.sizePredicted, model, thixotropic, dt, P.ptt);
        momentumStress(
          nodal, nodal.psiPredicted, model, fluid, P.mu, dt, P.devss, P.implicit);
      }
      else
      {
        // The viscosity and the stress of CY, from u*.
        toDevice(uStar.data(), nodal.velocity);
        recoverGradient(nodal);
        momentumStress(nodal, nodal.psi, model, fluid, P.mu, dt, P.devss, P.implicit);
        diagnostics = diagnose(nodal, model, false);
      }
      viscosity.setData(toHost(nodal.viscosity));
      increment.setData(toHost(nodal.increment));
      extra0.setData(toHost(nodal.extra0));
      extra1.setData(toHost(nodal.extra1));
      extra2.setData(toHost(nodal.extra2));
      const auto r1 = std::chrono::steady_clock::now();

      // Step 4. tau by cell, from u*, the viscosity and the cell volume
      // |det J|/6 of a tetrahedron.
      using VectorView =
        ::Rodin::Kokkos::Variational::View<std::decay_t<decltype(velocity)>>;
      using ScalarView =
        ::Rodin::Kokkos::Variational::View<std::decay_t<decltype(viscosity)>>;
      const VectorView w(advecting);
      const VectorView old(velocity);
      const VectorView row0(extra0);
      const VectorView row1(extra1);
      const VectorView row2(extra2);
      const ScalarView eta(viscosity);
      const ScalarView etaD(increment);
      RealFunction tau(KOKKOS_LAMBDA(const PointView& point)->Real {
        const auto wp = w.getValue(point);
        Real speed2 = 0;
        for (Index i = 0; i < 3; ++i)
          speed2 += wp(i) * wp(i);
        const Real hk = ::Kokkos::cbrt(point.getDistortion() / 6);
        const Real b2 = 4 * speed2 / (hk * hk);
        const Real c = 4 * eta.getValue(point) / (rho * hk * hk);
        return 1 / ::Kokkos::sqrt(b2 + c * c);
      });
      RealFunction pressureTau(KOKKOS_LAMBDA(const PointView& point)->Real {
        const auto wp = w.getValue(point);
        Real speed2 = 0;
        for (Index i = 0; i < 3; ++i)
          speed2 += wp(i) * wp(i);
        const Real hk = ::Kokkos::cbrt(point.getDistortion() / 6);
        const Real b2 = 4 * speed2 / (hk * hk);
        const Real c = 4 * eta.getValue(point) / (rho * hk * hk);
        return 1 / (rho * ::Kokkos::sqrt(b2 + c * c));
      });
      // div E, constant by cell: the force of the residual of the subscales.
      VectorFunction force(
        KOKKOS_LAMBDA(const PointView& point)->Real {
          return row0.u.interpolateDivergence(point);
        },
        KOKKOS_LAMBDA(const PointView& point)->Real {
          return row1.u.interpolateDivergence(point);
        },
        KOKKOS_LAMBDA(const PointView& point)->Real {
          return row2.u.interpolateDivergence(point);
        });
      // Column k of E - eta_d grad u^n, with the cell gradient of u^n.
      using Entry = StressEntry<VectorView, ScalarView>;
      VectorFunction column0(Entry{row0, etaD, old, 0, 0}, Entry{row1, etaD, old, 1, 0},
        Entry{row2, etaD, old, 2, 0});
      VectorFunction column1(Entry{row0, etaD, old, 0, 1}, Entry{row1, etaD, old, 1, 1},
        Entry{row2, etaD, old, 2, 1});
      VectorFunction column2(Entry{row0, etaD, old, 0, 2}, Entry{row1, etaD, old, 1, 2},
        Entry{row2, etaD, old, 2, 2});

      // Galerkin terms, openings and history; then the two diagonal terms of
      // the stabilisation and the two that complete the residual; then the
      // terms of the explicit stress, which each fluid assembles or not.
      const auto galerkin = (rho / dtEff) * Integral(u, v) +
        (0.5 * rho) * Integral(Dot(Mult(Jacobian(u), advecting), v)) -
        (0.5 * rho) * Integral(Dot(u, Mult(Jacobian(v), advecting))) +
        Integral(Mult(viscosity, Jacobian(u)), Jacobian(v)) - Integral(p, Div(v)) +
        Integral(Div(u), q) +
        pIn * BoundaryIntegral(normal, v).over(veins[0], veins[1], veins[2], veins[3]) +
        pOut * BoundaryIntegral(normal, v).over(mitral) +
        P.inletImpedance *
          BoundaryIntegral(u, v).over(veins[0], veins[1], veins[2], veins[3]) +
        mitralImpedance * BoundaryIntegral(Dot(normal, u), Dot(normal, v)).over(mitral) -
        (rho / dtEff) * Integral(past, v);
      const auto diagonal = rho *
          Integral(
            Dot(Mult(tau, Mult(Jacobian(u), advecting)), Mult(Jacobian(v), advecting))) +
        Integral(Dot(Mult(pressureTau, Grad(p)), Grad(q)));
      const auto cross = Integral(Dot(Mult(tau, Grad(p)), Mult(Jacobian(v), advecting))) +
        Integral(Dot(Mult(tau, Mult(Jacobian(u), advecting)), Grad(q)));
      const auto stressResidual =
        Integral(Dot(Mult(tau, force), Mult(Jacobian(v), advecting))) +
        Integral(Dot(Mult(pressureTau, force), Grad(q)));
      const auto stressWeak = Integral(Dot(column0, Mult(Jacobian(v), axis0))) +
        Integral(Dot(column1, Mult(Jacobian(v), axis1))) +
        Integral(Dot(column2, Mult(Jacobian(v), axis2)));
      const auto noSlip = DirichletBC(u, Zero(3)).on(wall);
      if (fluid == Fluid::Newtonian)
      {
        if (P.residual)
          flow = galerkin + diagonal + cross + noSlip;
        else
          flow = galerkin + diagonal + noSlip;
      }
      else if (elastic)
      {
        // E enters as int E : grad v, the adjoint of the vertex gradient.
        if (P.residual)
          flow = galerkin + diagonal + cross + stressWeak - stressResidual + noSlip;
        else
          flow = galerkin + diagonal + stressWeak + noSlip;
      }
      else
      {
        // CY: E enters as int div E . v.
        if (P.residual)
          flow =
            galerkin + diagonal + cross - Integral(force, v) - stressResidual + noSlip;
        else
          flow = galerkin + diagonal - Integral(force, v) + noSlip;
      }

      const auto t0 = std::chrono::steady_clock::now();
      flow.assemble();
      if (!fieldSplits)
      {
        // Named splits u and p for -pc_type fieldsplit; harmless otherwise.
        flow.setFieldSplits();
        fieldSplits = true;
      }
      const auto t1 = std::chrono::steady_clock::now();
      if (P.preconditionerLag > 1)
      {
        // Kept between two set-ups, except when the valve changes the matrix.
        const bool keep =
          (step - 1) % P.preconditionerLag != 0 && valveOpen == valveWasOpen;
        PetscCallAbort(PETSC_COMM_SELF,
          KSPSetReusePreconditioner(solver.getHandle(), keep ? PETSC_TRUE : PETSC_FALSE));
        valveWasOpen = valveOpen;
      }
      flow.solve(solver);
      PetscInt iterations = 0;
      KSPConvergedReason reason;
      PetscCallAbort(PETSC_COMM_SELF, KSPGetConvergedReason(solver.getHandle(), &reason));
      if (reason < 0)
      {
        // BiCGStab(l) can break down: one more attempt, GMRES from a zero guess.
        char type[64] = "";
        PetscBool given = PETSC_FALSE;
        PetscCallAbort(PETSC_COMM_SELF,
          PetscOptionsGetString(
            PETSC_NULLPTR, PETSC_NULLPTR, "-ksp_type", type, sizeof(type), &given));
        std::cerr << "Step " << step << ": " << KSPConvergedReasons[reason]
                  << ", solving again with GMRES" << std::endl;
        PetscCallAbort(
          PETSC_COMM_SELF, PetscOptionsSetValue(PETSC_NULLPTR, "-ksp_type", KSPGMRES));
        PetscBool restart = PETSC_FALSE;
        PetscCallAbort(PETSC_COMM_SELF,
          PetscOptionsHasName(
            PETSC_NULLPTR, PETSC_NULLPTR, "-ksp_gmres_restart", &restart));
        if (!restart)
        {
          PetscCallAbort(PETSC_COMM_SELF,
            PetscOptionsSetValue(PETSC_NULLPTR, "-ksp_gmres_restart", "100"));
        }
        PetscCallAbort(PETSC_COMM_SELF, VecSet(flow.getLinearSystem().getSolution(), 0));
        solver.setType(KSPGMRES);
        flow.solve(solver);
        solver.setType(KSPBCGSL);
        if (given)
        {
          PetscCallAbort(
            PETSC_COMM_SELF, PetscOptionsSetValue(PETSC_NULLPTR, "-ksp_type", type));
        }
        else
        {
          PetscCallAbort(
            PETSC_COMM_SELF, PetscOptionsClearValue(PETSC_NULLPTR, "-ksp_type"));
        }
        PetscCallAbort(
          PETSC_COMM_SELF, KSPGetConvergedReason(solver.getHandle(), &reason));
        ++fallbacks;
      }
      const auto t2 = std::chrono::steady_clock::now();
      PetscCallAbort(
        PETSC_COMM_SELF, KSPGetIterationNumber(solver.getHandle(), &iterations));
      if (reason < 0)
      {
        std::cerr << "The linear solver failed at step " << step << ": "
                  << KSPConvergedReasons[reason] << " after " << iterations
                  << " iterations" << std::endl;
        return EXIT_FAILURE;
      }
      velocity.setData(toHost(u.getSolution().getData()));
      pressure.setData(toHost(p.getSolution().getData()));
      const HostVector uh = velocity.getData().toHost();
      const auto ph = pressure.getData().toHost();

      // Step 5: correction of psi and n with the gradient of the midpoint.
      const auto r2 = std::chrono::steady_clock::now();
      if (elastic)
      {
        const HostVector uMid = second ? HostVector(0.5 * (un + uh)) : uh;
        toDevice(uMid.data(), nodal.velocity);
        recoverGradient(nodal);
        advanceStructure(nodal, nodal.psi, nodal.size, model, thixotropic, dt, P.ptt);
        diagnostics = diagnose(nodal, model, true);
      }
      const auto r3 = std::chrono::steady_clock::now();
      const double rheologySeconds =
        std::chrono::duration<double>((r1 - r0) + (r3 - r2)).count();
      const double assemblySeconds = std::chrono::duration<double>(t1 - t0).count();
      const double solveSeconds = std::chrono::duration<double>(t2 - t1).count();
      rheologyTime += (r1 - r0) + (r3 - r2);
      assembly += t1 - t0;
      solve += t2 - t1;

      // Diagnostics on the host copies of the solution. n is outward: the
      // venous flux is negated to read as the inflow.
      const Real qIn = -veinFluxHost.dot(uh);
      const Real qOut = mitralFluxHost.dot(uh);
      const Real pVein = veinMeanHost.dot(ph) / veinArea;
      const Real planePressure = mitralMeanHost.dot(ph) / mitralArea;
      // With the skew-symmetric convection, the natural condition of the mitral
      // plane balances p + rho (u.n)^2/2, not p, against p_out + Z u.n. The
      // ventricle reads the mean of that balance, p_out + Q/K: its inflow
      // K (p_at - p_v) is then the flow rate of the plane, and the valve
      // closes when the flow reverses rather than when the jet lowers p.
      atrialPressure =
        P.ventricle ? pOut + mitralImpedance * qOut / mitralArea : planePressure;
      const Real speed = uh.cwiseAbs().maxCoeff();
      if (!std::isfinite(speed) || speed > 50)
      {
        std::cerr << "Diverged at step " << step << ": max |u_i| = " << speed << " m/s"
                  << std::endl;
        return EXIT_FAILURE;
      }
      // Kinetic energy, and the energy removed by the time scheme per unit
      // time: rho/(2 dt) |u^{n+1} - u^n|^2 for backward Euler, rho/(4 dt)
      // |u^{n+1} - 2 u^n + u^{n-1}|^2 for BDF2 (lumped mass).
      Real kinetic = 0, timeDissipation = 0;
      for (size_t i = 0; i < nodes; ++i)
      {
        for (size_t c = 0; c < 3; ++c)
        {
          const size_t k = i + c * nodes;
          const Real jump =
            second ? uh(k) - 2 * un(k) + previousVelocity(k) : uh(k) - un(k);
          kinetic += 0.5 * rho * lumpedVolume(i) * uh(k) * uh(k);
          timeDissipation += lumpedVolume(i) * jump * jump;
        }
      }
      timeDissipation *= rho / ((second ? 4 : 2) * dt);
      previousVelocity = un;
      const size_t cycle = (step - 1) / P.steps + 1;

      csv << time << ',' << cycle << ',' << pIn << ',' << pOut << ','
          << (valveOpen ? 1 : 0) << ',' << qIn << ',' << qOut << ',' << pVein << ','
          << planePressure << ',' << speed << ',' << kinetic << ',' << timeDissipation
          << ',' << diagnostics.solventDissipation << ',' << diagnostics.polymerPower
          << ',' << diagnostics.elasticEnergy << ',' << diagnostics.maxWeissenberg << ','
          << diagnostics.maxTrace << ',' << diagnostics.meanSize << ',' << rheologySeconds
          << ',' << assemblySeconds << ',' << solveSeconds << ',' << iterations << '\n';
      if (step % P.steps == 0 || step == last)
        writeCheckpoint(step, uh);
      if (step % P.outputEvery == 0 || step == last)
      {
        shearRate.setData(toHost(nodal.shearRate));
        if (elastic)
        {
          relaxation.setData(toHost(nodal.relaxation));
          aggregate.setData(toHost(nodal.size));
          trace.setData(toHost(nodal.trace));
          stress.setData(toHost(nodal.stress));
          const HostVector psi = toHost(nodal.psi);
          psiDiagonal.setData(HostVector(psi.head(3 * nodes)));
          psiOffDiagonal.setData(HostVector(psi.tail(3 * nodes)));
        }
        xdmf.write(time).flush();
        csv.flush();
        std::cout << "Step " << step << " / " << last << "  t=" << time
                  << " s  p_pv=" << pIn << "  p_out=" << pOut << " Pa  valve "
                  << (valveOpen ? "open" : "closed") << "  qIn=" << qIn
                  << "  qOut=" << qOut << " m^3/s  max|u_i|=" << speed << " m/s";
        if (elastic)
        {
          std::cout << "  Wi_max=" << diagnostics.maxWeissenberg
                    << "  tr(c)_max=" << diagnostics.maxTrace
                    << "  n_mean=" << diagnostics.meanSize;
        }
        std::cout << "  rheology " << rheologySeconds << " s  assembly "
                  << assemblySeconds << " s  solve " << solveSeconds
                  << " s  its=" << iterations << std::endl;
      }
    }
    std::cout << "Rheology " << rheologyTime.count() << " s, assembly "
              << assembly.count() << " s, solve " << solve.count() << " s in " << total
              << " steps, " << fallbacks << " solved again with GMRES" << std::endl;
  }
  PetscFinalize();
  return 0;
}
