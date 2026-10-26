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
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include "ProblemData.h"

namespace Rodin::Examples::Heart::AtriumContraction
{
  class Output;

  /**
   * @brief Transient active atrial wall on a local wedge mesh.
   *
   * ## Architecture
   *
   * The component retains its reference mesh, reference fibers, active fields, finite-element
   * spaces, active material, nonlinear problem, and solver between time steps.
   * Newmark updates advance displacement, velocity, and acceleration. Output
   * observes the retained fields. This component uses the local Eigen backend;
   * its mesh and interface are independent of the MPI fluid component.
   */
  class Solid
  {
    public:
      using Parameters = Options;

      /**
       * @brief Constructs the wall component and prescribed material inputs.
       * @param parameters Validated wall, material, wave, and time configuration.
       */
      explicit Solid(const Parameters& parameters);
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
      const Geometry::LocalMesh& getMesh() const noexcept;

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
      using PassiveMaterial = Rodin::Solid::KelvinVoigt<Rodin::Solid::DemirayGasser>;
      using Material =
        Rodin::Solid::ActiveContraction<PassiveMaterial, Rodin::Solid::ActiveFiberLaw>;
      using Trial = Variational::TrialFunction<VectorField, VectorSpace>;
      using Test = Variational::TestFunction<VectorSpace>;
      using ActiveTrial = Variational::TrialFunction<CellScalarField, CellScalarSpace>;
      using ActiveTest = Variational::TestFunction<CellScalarSpace>;
      using LinearSystem =
        Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
      using MechanicsProblem = Variational::Problem<LinearSystem, Trial, ActiveTrial,
        ActiveTrial, ActiveTrial, Test, ActiveTest, ActiveTest, ActiveTest>;
      using LinearSolver = Solver::SparseLU<LinearSystem>;
      using Newton = Solver::NewtonSolver<LinearSolver>;
      using InternalWork = Rodin::Solid::InternalVirtualWork<Material, VectorField>;

      static constexpr Real NewmarkBeta = 0.25;
      static constexpr Real NewmarkGamma = 0.5;
      // 1 kPa s²/mm² = 10⁹ kg/m³.
      static constexpr Real DensityConversion = 1e-9;
      // Quadrature policy for nonlinear constitutive terms on nonaffine wedges.
      static constexpr size_t QuadratureOrder = 4;
      static constexpr size_t MaximumNewtonIterations = 80;
      static constexpr size_t MaximumBacktracks = 60;
      static constexpr Real AbsoluteResidualTolerance = 1e-8;

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
       * @brief Prepares the coupled material at a supplied reference point.
       * @param point Reference evaluation point.
       * @returns Material with owned kinematics and active residuals.
       */
      Material prepare(const Geometry::Point& point) const;

      /// @brief Packs mechanical and independent active fields into the Newton vector.
      void pack();

      /// @brief Scatters the current Newton vector to its finite-element fields.
      void scatter();

      /// @brief Returns whether the current iterate has positive stretches and wedges.
      bool isAdmissible() const;

      Parameters m_parameters;
      Geometry::LocalMesh m_mesh;
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
      Newton m_newton;
      Math::Vector<Real> m_iterate;
      std::unique_ptr<Output> m_output;
      bool m_initialized = false, m_finished = false, m_failed = false;
  };
}

#endif
