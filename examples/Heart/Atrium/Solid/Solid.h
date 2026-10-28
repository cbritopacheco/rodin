/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_SOLID_SOLID_H
#define RODIN_ATRIUM_SOLID_SOLID_H

#include <memory>
#include <Rodin/Solid.h>
#include <Rodin/PETSc/Solver/SNES.h>
#include "ProblemData.h"

namespace Rodin::Examples::Heart::Atrium
{
  class Output;

  /**
   * @brief Transient active atrial wall on a Scotch-partitioned MPI wedge mesh.
   *
   * ## Architecture
   *
   * The component retains its reference mesh, reference fibers, active fields, finite-element
   * spaces, active material, nonlinear problem, and solver between time steps.
   * Newmark updates advance displacement, velocity, and acceleration. Output
   * observes the retained fields. Distributed fields use PETSc storage and KSP/SNES solvers. The solid
   * and fluid retain their own partitioned meshes and assembly state.
   */
  class Solid
  {
    public:
      using Parameters = Options;
      using Diagnostics = DiagnosticState;

      /// @brief Returns cached scalar diagnostics, or empty before a mechanical state is written.
      const Optional<Diagnostics>& getDiagnostics() const noexcept
      {
        return m_diagnostics;
      }

      /**
       * @brief Constructs the wall component and prescribed material inputs.
       * @param parent Combined parent mesh, already distributed by the coordinator.
       * @param material Solid volume label to extract.
       * @param parameters Validated wall, material, wave, and time configuration.
       */
      Solid(const MeshType& parent, const Parameters& parameters,
        Geometry::Attribute material = 4);
      ~Solid();
      Solid(const Solid&) = delete;
      Solid& operator=(const Solid&) = delete;

      /**
       * @brief Solves passive equilibrium and writes the relaxed initial state.
       * @returns This component.
       */
      Solid& initialize();

      /**
       * @brief Advances one Newmark time step and writes the state.
       * @note A failed step cannot be retried on the same component.
       */
      void step();

      /// @brief Finalizes output; no further time steps may be taken.
      void finish();

      /// @brief Advances the configured number of intervals and finalizes output.
      void run();

      /// @brief Current physical time in seconds.
      Real getTime() const noexcept;

      /// @brief Configured time increment in seconds.
      Real getTimeStep() const noexcept;

      /// @brief Number of successfully completed time steps.
      size_t getStepCount() const noexcept;

      /// @brief Returns the wall configuration.
      const Parameters& getParameters() const noexcept;

      /// @brief Returns the reference wall mesh, with coordinates in mm.
      const SubMeshType& getMesh() const noexcept;

      /// @brief Returns the current displacement in mm.
      const VectorField& getDisplacement() const noexcept;

      /// @brief Returns the current velocity in mm/s.
      const VectorField& getVelocity() const noexcept;

      /// @brief Returns the current acceleration in mm/s².
      const VectorField& getAcceleration() const noexcept;

      /// @brief Returns the solved discontinuous active extension.
      const CellScalarField& getActiveExtension() const noexcept;

      /// @brief Returns the square root of contractile stiffness.
      const CellScalarField& getActiveGamma() const noexcept;

      /// @brief Returns contractile stress divided by the square root of stiffness.
      const CellScalarField& getActiveBeta() const noexcept;

    private:
      using PassiveMaterial =
        Rodin::Solid::KelvinVoigt<Rodin::Solid::DemirayHolzapfelGasserOgden>;
      using Material =
        Rodin::Solid::ActiveContraction<PassiveMaterial, Rodin::Solid::ActiveFiberLaw>;
      using Trial = PETSc::Variational::TrialFunction<VectorField, VectorSpace>;
      using Test = PETSc::Variational::TestFunction<VectorSpace>;
      using ActiveTrial =
        PETSc::Variational::TrialFunction<CellScalarField, CellScalarSpace>;
      using ActiveTest = PETSc::Variational::TestFunction<CellScalarSpace>;
      using LinearSystem = PETSc::Math::LinearSystem;
      using MechanicsProblem = Variational::Problem<LinearSystem, Trial, ActiveTrial,
        ActiveTrial, ActiveTrial, Test, ActiveTest, ActiveTest, ActiveTest>;
      using EquilibriumProblem = Variational::Problem<LinearSystem, Trial, Test>;
      using LinearSolver = Solver::KSP;
      using Newton = Solver::SNES;
      using InternalWork = Rodin::Solid::InternalVirtualWork<Material, VectorField>;
      using PassiveWork = Rodin::Solid::InternalVirtualWork<PassiveMaterial, VectorField>;

      static constexpr Real NewmarkBeta = 0.25;
      static constexpr Real NewmarkGamma = 0.5;
      // 1 kPa s²/mm² = 10⁹ kg/m³.
      static constexpr Real DensityConversion = 1e-9;
      // Quadrature policy for nonlinear constitutive terms on nonaffine wedges.
      static constexpr size_t QuadratureOrder = 4;
      // Heuristic work caps; reaching either limit fails the time step.
      static constexpr size_t MaximumNewtonIterations = 80;

      /// Function-evaluation budget allowing line-search backtracking at each Newton step.
      static constexpr size_t MaximumFunctionEvaluations = 10000;
      // Absolute assembled residual tolerance in the example's mm/kPa/s units.
      // Relative stopping is disabled because the initial active residual can
      // be large while a mechanically unresolved residual remains significant.
      static constexpr Real AbsoluteResidualTolerance = 1e-8;

      /**
       * @brief Collectively extracts the solid cells from the common parent shard.
       * @param parent Distributed mesh owned by the atrium coordinator.
       * @param material Volume material label selected for the solid.
       * @returns Solid submesh retaining its parent ancestry maps.
       */
      static SubMeshType makeMesh(const MeshType& parent, Geometry::Attribute material);

      /**
       * @brief Constructs passive elasticity, tissue viscosity, and active kinetics.
       * @param cfg Solid material and active coefficients.
       * @returns Configured constitutive law.
       */
      static Material makeMaterial(const Parameters& cfg);

      /**
       * @brief Populates rate, fiber, and independent active-field inputs.
       * @param point Constitutive point to populate.
       */
      void setInputs(Rodin::Solid::ConstitutivePoint& point) const;

      /**
       * @brief Active residual, tangent, and momentum couplings at one material point.
       *
       * Owns only the blocks used by the mixed coefficient expressions. The
       * reference mesh is fixed; a synchronization epoch identifies the current
       * fields, history, predictors, and activation time. Entries are refreshed
       * after each Newton state scatter and before each time-step assembly.
       */
      struct Prepared
      {
        size_t epoch = 0;
        Rodin::Solid::ActiveFiberLaw::Residual residual;
        Math::SpatialMatrix<Real> momentum{3, 3}, extension{3, 3};
      };

      /**
       * @brief Prepares the coupled material at a supplied reference point.
       * @param point Reference evaluation point.
       * @returns Owned active blocks valid until the next preparation or state update.
       */
      const Prepared& prepare(const Geometry::Point& point) const;

      /// @brief State epoch for constitutive coefficients; advanced on every scatter.
      size_t m_constitutiveEpoch = 1;

      /// @brief Material-point keys use cell identity and exact reference coordinates.
      using ConstitutiveKey = std::pair<Index, std::array<Real, 3>>;

      /// @brief Rank-local reuse across coefficient expressions on the fixed mesh.
      mutable UnorderedMap<ConstitutiveKey, Prepared> m_constitutiveCache;

      /// @brief Owned scratch storage when constitutive reuse is disabled.
      mutable Prepared m_uncachedPrepared;

      /// @brief Rank-local preparation requests and actual constitutive evaluations.
      mutable size_t m_prepareRequests = 0, m_prepareEvaluations = 0;

      /**
       * @brief Initializes active fields with the residual's volume quadrature.
       *
       * The cellwise stiffness predictor averages the pointwise squared gamma
       * updates before taking the square root. Inactive quadrature samples thus
       * retain the constitutive floor even when a cell vertex is activated.
       * Beta uses the volume mean of its pointwise update. These are initial
       * guesses only; all three active fields remain in the global Newton solve.
       */
      void predictActiveState();

      /// @brief Packs mechanical and independent active fields into the Newton vector.
      void pack();

      /**
       * @brief Synchronizes the four retained fields from a global SNES iterate.
       * @param values Distributed block iterate with field offsets defined by the spaces.
       */
      void scatter(const PETSc::Math::Vector& values);

      /// @brief Returns whether the current iterate has positive stretches and wedges.
      bool isAdmissible() const;

      Parameters m_parameters;
      SubMeshType m_mesh;
      ProblemData m_data;
      VectorSpace m_space;
      CellScalarSpace m_activeSpace;
      VectorField m_displacement, m_velocity, m_acceleration;
      VectorField m_predictedDisplacement, m_predictedVelocity, m_previousDisplacement;
      CellScalarField m_extension, m_gamma, m_beta;
      CellScalarField m_previousExtension, m_previousGamma, m_previousBeta;
      Real m_timeStep, m_density, m_accelerationCoefficient, m_massCoefficient;
      Real m_velocityCoefficient;
      Real m_time = 0;
      size_t m_step = 0;
      Material m_law;
      Trial m_increment;
      ActiveTrial m_extensionIncrement, m_gammaIncrement, m_betaIncrement;
      Test m_test;
      ActiveTest m_extensionTest, m_gammaTest, m_betaTest;
      InternalWork m_work;
      MechanicsProblem m_mechanics;
      LinearSolver m_linearSolver;
      Optional<Newton> m_newton;
      PassiveWork m_passiveWork;
      EquilibriumProblem m_equilibrium;
      LinearSolver m_equilibriumSolver;
      Optional<Newton> m_equilibriumNewton;
      Optional<Diagnostics> m_diagnostics;
      std::unique_ptr<Output> m_output;
      bool m_initialized = false, m_finished = false, m_failed = false;
  };
}

#endif
