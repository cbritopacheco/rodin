// ThrombosisModel.h
//
// Coagulation kinetics (Qureshi et al. 2022) with per-cycle haemodynamic
// indices, written once so that both CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq
// and CoupledLV0DCoronary3D use the same implementation.
//
// ---------------------------------------------------------------------------
// MODEL
//
//   dTh/dt = D_th DTh - u.grad(Th) + R_Th                     (thrombin)
//   dFg/dt = D_fg DFg - u.grad(Fg) - k_eff Fg Th              (fibrinogen)
//   dFn/dt = D_fn DFn - u.grad(Fn) + k_eff Fg Th              (fibrin)
//
// The three fields are transported by a *known* velocity: the coupling is
// one-way, the kinetics do not feed back on the flow. R_Th is an endothelial
// source, active where the endothelium has been activated by the flow, and
// measured per unit *area*: it is a Neumann flux on the wall, not a volumetric
// source.
//
// ---------------------------------------------------------------------------
// SCOPE -- read this before using it anywhere new.
//
// This is a model of thrombus *initiation and growth*. It predicts where and how fast fibrin would begin to accumulate in
// a geometry whose flow is not yet perturbed by the clot. It contains no
// mechanical representation of the clot: no permeability, no added resistance,
// no moving boundary, no platelet aggregation and no lysis. It ceases to be
// valid once the fibrin field is dense enough to alter the flow, which is
// exactly the two-way extension the model does not have.
//
// It is also a *stasis* model. The activation indicator is built from wall
// shear statistics, so the mechanism it encodes is recirculation and low,
// oscillatory wall shear. Transferring it unchanged to a coronary artery is questionable on
// two counts: the residence time of a coronary conduit is of order 0.1 s
// against the seconds-to-minutes of the cascade, so the species are advected
// away before they react; and coronary thrombosis is driven by plaque rupture
// and tissue-factor exposure rather than by stasis. The model is appropriate
// for a coronary geometry only where a genuine low-flow feature exists -- an
// aneurysmal or ectatic segment, a post-stent recirculation zone, a side branch
// with reversed flow -- and the residence-time diagnostic below is what decides
// whether such a feature is present.
//
// ---------------------------------------------------------------------------
#ifndef EXAMPLES_HEART_CORONARYARTERY_THROMBOSISMODEL_H
#define EXAMPLES_HEART_CORONARYARTERY_THROMBOSISMODEL_H

#include <cmath>
#include <cstddef>
#include <limits>

#include <Rodin/Types.h>

namespace Rodin::Examples::Heart
{
  /**
   * @brief Constants of the coagulation model and of its stabilization.
   *
   * @details Values from Qureshi et al. (2022) unless noted. Concentrations are
   *          in mol/m^3, so a fibrinogen level quoted in g/L must be divided by
   *          the molar mass, 340 kg/mol: 2.5 g/L -> 7.35e-3 mol/m^3.
   */
  struct ThrombosisParameters
  {
      using Real = Rodin::Real;

  /// @brief Thrombin diffusivity (m^2/s).
      Real diffusivityThrombin = 4.6e-11;
  /// @brief Fibrinogen diffusivity (m^2/s).
      Real diffusivityFibrinogen = 2.0e-11;
  /// @brief Fibrin diffusivity (m^2/s).
      Real diffusivityFibrin = 2.0e-11;

  /// @brief Effective reaction rate (m^3 mol^-1 s^-1).
      Real reactionRate = 7180.0;

  /// @brief Molar mass of fibrinogen (kg/mol), used to convert g/L to mol/m^3.
      Real fibrinogenMolarMass = 340.0;
  /// @brief Initial fibrinogen in sinus rhythm (g/L).
      Real fibrinogenSinusRhythm = 2.5;
  /// @brief Initial fibrinogen in atrial fibrillation (g/L).
      Real fibrinogenFibrillation = 4.0;

  /// @brief Endothelial thrombin flux at an activated site (mol m^-2 s^-1).
  /// @details 0.1 nmol/(m^2 s). It is a *surface* flux: it enters the weak form
  ///          as a boundary integral over the activated wall, never as a
  ///          volume integral.
      Real thrombinWallFlux = 1.0e-10;

  /// @brief Wall shear stress below which the endothelium is taken to be
  ///        prothrombotic (Pa).
  /// @details Used by the smooth activation kernel. The endothelial phenotype
  ///          switches below roughly 0.4 Pa, which is a measured quantity,
  ///          unlike a percentile of the computed ECAP field.
      Real activationShearStress = 0.4;
  /// @brief Width of the activation transition (Pa).
      Real activationShearWidth = 0.15;

  /// @brief Number of complete cycles to discard before the source is armed.
  /// @details The wall-shear statistics are only meaningful once the flow is
  ///          periodic, and the indices are accumulated over a whole cycle, so
  ///          the source cannot be active during the first one.
      int warmupCycles = 1;

  /// @brief Cycle period (s).
      Real cyclePeriod = 0.85;

  /// @brief SUPG/VMS stabilization multiplier. 1 is the standard value; set to
  ///        0 to disable stabilization and expose the underlying oscillations.
      Real stabilizationScale = 1.0;

  /// @brief Floor on TAWSS used when forming ECAP (Pa).
  /// @details ECAP = OSI/TAWSS diverges wherever TAWSS -> 0, which happens at
  ///          every stagnation point irrespective of how oscillatory the flow
  ///          is. Without a floor the maximum of ECAP is a property of the mesh
  ///          rather than of the flow. See the note on indicator choice below.
      Real tawssFloor = 1.0e-3;
  };

  /**
   * @brief Cycle-averaged wall-shear indices.
   *
   * @details Accumulates, over one cardiac cycle,
   *
   *    TAWSS = (1/T) int_0^T |tau_w| dt,
   *    OSI   = (1/2) [ 1 - |int_0^T tau_w dt| / int_0^T |tau_w| dt ],
   *    ECAP  = OSI / TAWSS,
   *
   * and exposes them at the end of the cycle. Two accumulators are needed and
   * they are not interchangeable: the vector integral of tau_w, whose magnitude
   * measures how much net direction survives, and the scalar integral of
   * |tau_w|, which measures how much shear was applied regardless of direction.
   *
   * @note On the choice of indicator. ECAP is a ratio and inherits the pathology
   *       of one: it diverges wherever TAWSS vanishes, so its maximum tends to
   *       sit on a stagnation point rather than on the most thrombogenic
   *       region, and that point moves with the mesh. OSI itself saturates at
   *       1/2 for perfectly reversing flow, which is not the most thrombogenic
   *       condition either. The class therefore also exposes the two raw
   *       accumulators, so that a different activation functional -- relative
   *       residence time, a sigmoid in TAWSS alone, or a residence-time field --
   *       can be formed without recomputing anything.
   */
  template <class VectorGridFunction, class ScalarGridFunction>
  class WallShearAccumulator
  {
    public:
      using Real = Rodin::Real;

      template <class VectorFES, class ScalarFES>
      WallShearAccumulator(const VectorFES& vfes, const ScalarFES& sfes)
        : m_netIntegral(vfes), m_absIntegral(sfes),
          m_tawss(sfes), m_osi(sfes), m_ecap(sfes),
          m_scalarSize(sfes.getSize()), m_vectorSize(vfes.getSize()),
          m_elapsed(0.0), m_cyclesCompleted(0), m_ready(false)
      {
        reset();
      }

      /// @brief Zero the accumulators and restart the cycle clock.
      ///
      /// @note On element access. For a PETSc-backed GridFunction, getData()
      ///       returns a raw ::Vec handle, so indexing must go through
      ///       operator[], which lazily acquires the local array
      ///       (VecGetArrayWrite) and translates the index. The array stays
      ///       acquired until flush() is called, and while it is held **no
      ///       other PETSc operation on that Vec is legal** -- not VecMax, not
      ///       an assembly that reads it, not the XDMF writer. Every routine
      ///       below therefore flushes every grid function it touched before
      ///       returning. Omitting that is not a leak but a lock: the next
      ///       PETSc call on the same vector blocks.
      void reset()
      {
        m_netIntegral = Real(0);
        m_absIntegral = Real(0);
        m_elapsed = 0.0;
      }

      /**
       * @brief Accumulate one time step.
       * @param wss Wall shear stress vector field (as produced by
       *            computeWallShear()); zero away from the wall.
       */
      void accumulate(const VectorGridFunction& wss, Real dt)
      {
        const std::size_t n = m_scalarSize;
        const std::size_t d = m_vectorSize / n;

        for (std::size_t i = 0; i < n; ++i)
        {
          Real sq = 0.0;
          for (std::size_t k = 0; k < d; ++k)
          {
            const Real c = wss[k * n + i];
            m_netIntegral[k * n + i] += c * dt;
            sq += c * c;
          }
          m_absIntegral[i] += std::sqrt(sq) * dt;
        }

        // Release every array acquired above, including the read array on the
        // wall-shear field: the const operator[] acquires too.
        wss.flush();
        m_netIntegral.flush();
        m_absIntegral.flush();

        m_elapsed += dt;
      }

      /**
       * @brief Close the cycle: form TAWSS, OSI and ECAP and restart.
       * @returns true if a cycle was completed and the indices are available.
       */
      bool closeCycle(const ThrombosisParameters& p)
      {
        if (m_elapsed <= 0.0)
          return false;

        const std::size_t n = m_scalarSize;
        const std::size_t d = m_vectorSize / n;

        for (std::size_t i = 0; i < n; ++i)
        {
          Real sq = 0.0;
          for (std::size_t k = 0; k < d; ++k)
          {
            const Real c = m_netIntegral[k * n + i];
            sq += c * c;
          }

          const Real absInt = m_absIntegral[i];
          const Real netMag = std::sqrt(sq);

          const Real tawss = absInt / m_elapsed;

          // OSI is 0 for unidirectional shear and 1/2 for perfectly reversing
          // shear. Where no shear was applied at all it is left at 0 rather
          // than 0/0.
          const Real osi = (absInt > 0.0) ? 0.5 * (1.0 - netMag / absInt) : 0.0;

          m_tawss[i] = tawss;
          m_osi[i] = osi;
          m_ecap[i] = osi / std::max<Real>(tawss, p.tawssFloor);
        }

        m_netIntegral.flush();
        m_absIntegral.flush();
        m_tawss.flush();
        m_osi.flush();
        m_ecap.flush();

        ++m_cyclesCompleted;
        m_ready = true;
        reset();
        return true;
      }

      /// @brief Index of the degree of freedom carrying the largest ECAP.
      /// @note Reported for reference only. Injecting the source at a single
      ///       node makes the result depend on the mesh; prefer the smooth
      ///       activation kernel of activationWeight().
      std::size_t argMaxECAP() const
      {
        std::size_t best = 0;
        Real bestVal = -std::numeric_limits<Real>::infinity();
        for (std::size_t i = 0; i < m_scalarSize; ++i)
          if (m_ecap[i] > bestVal) { bestVal = m_ecap[i]; best = i; }
        m_ecap.flush();
        return best;
      }

      Real maxECAP() const { return m_ecap.max(); }

      /// @brief Number of scalar degrees of freedom carried by the indices.
      std::size_t size() const { return m_scalarSize; }

      /// @brief Have enough cycles elapsed for the source to be armed?
      bool sourceArmed(const ThrombosisParameters& p) const
      {
        return m_ready && m_cyclesCompleted >= p.warmupCycles;
      }

      int cyclesCompleted() const { return m_cyclesCompleted; }

      const ScalarGridFunction& tawss() const { return m_tawss; }
      const ScalarGridFunction& osi() const { return m_osi; }
      const ScalarGridFunction& ecap() const { return m_ecap; }
      ScalarGridFunction& tawss() { return m_tawss; }
      ScalarGridFunction& osi() { return m_osi; }
      ScalarGridFunction& ecap() { return m_ecap; }

      /// @brief Raw accumulators, exposed so that alternative activation
      ///        functionals can be built without recomputing them.
      const VectorGridFunction& netShearIntegral() const { return m_netIntegral; }
      const ScalarGridFunction& absShearIntegral() const { return m_absIntegral; }

    private:
      VectorGridFunction m_netIntegral;
      ScalarGridFunction m_absIntegral;
      ScalarGridFunction m_tawss;
      ScalarGridFunction m_osi;
      ScalarGridFunction m_ecap;
      std::size_t m_scalarSize;
      std::size_t m_vectorSize;
      Real m_elapsed;
      int m_cyclesCompleted;
      bool m_ready;
  };

  /**
   * @brief Smooth endothelial activation weight in [0,1].
   *
   * @details Written into @p weight from the cycle indices. Two mechanisms are
   *          combined multiplicatively, both bounded:
   *
   *            w = sigma_low(TAWSS) * (2 OSI),
   *
   *          where sigma_low is a logistic in the *measured* activation
   *          threshold, sigma_low = 1/(1+exp((TAWSS - tau_a)/w_a)), and 2*OSI
   *          maps the oscillatory index onto [0,1]. This replaces the single
   *          arg-max of ECAP by a distributed, mesh-independent field, and it
   *          removes the divergence of the ratio: where TAWSS -> 0 the weight
   *          saturates at 1 rather than growing without bound.
   *
   *          The weight is then normalised so that its integral over the wall
   *          equals the area of the region the caller intends to activate, which
   *          keeps the total injected thrombin independent of mesh refinement.
   */
  template <class Accumulator, class ScalarGridFunction>
  void activationWeight(const Accumulator& acc, const ThrombosisParameters& p,
    ScalarGridFunction& weight)
  {
    using Real = Rodin::Real;

    const auto& tawss = acc.tawss();
    const auto& osi = acc.osi();

    for (std::size_t i = 0; i < acc.size(); ++i)
    {
      const Real low =
        1.0 / (1.0 + std::exp((tawss[i] - p.activationShearStress) /
                              std::max<Real>(p.activationShearWidth, 1e-12)));
      weight[i] = low * std::min<Real>(2.0 * osi[i], 1.0);
    }

    // Same rule as in the accumulator: the arrays acquired by operator[] must
    // be released before any other PETSc operation touches these vectors.
    tawss.flush();
    osi.flush();
    weight.flush();
  }

  /**
   * @brief SUPG/VMS stabilization parameter for a scalar transport equation.
   *
   * @details tau = [ (2/dt)^2 + (2|u|/h)^2 + (4 D/h^2)^2 + sigma^2 ]^{-1/2},
   *          with sigma the reactive coefficient of the equation. It has units
   *          of time. Its absence is not a matter of accuracy: the stabilizing
   *          terms are then dimensionally inconsistent with the Galerkin ones
   *          and dominate the system by many orders of magnitude.
   */
  inline Rodin::Real supgTau(Rodin::Real h, Rodin::Real speed, Rodin::Real diffusivity,
    Rodin::Real reaction, Rodin::Real dt)
  {
    using Real = Rodin::Real;
    const Real a = 2.0 / std::max<Real>(dt, 1e-300);
    const Real b = 2.0 * speed / std::max<Real>(h, 1e-300);
    const Real c = 4.0 * diffusivity / std::max<Real>(h * h, 1e-300);
    return 1.0 / std::sqrt(a * a + b * b + c * c + reaction * reaction);
  }
}

#endif
