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
   * The component retains its reference mesh, prescribed fields, finite-element
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
       * @brief Computes the initial acceleration and writes the initial state.
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

    private:
      using Material = Rodin::Solid::ActiveContraction<Rodin::Solid::NeoHookean,
        Rodin::Solid::ActiveFiberLaw>;
      using Trial = Variational::TrialFunction<VectorField, VectorSpace>;
      using Test = Variational::TestFunction<VectorSpace>;
      using LinearSystem =
        Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
      using MechanicsProblem = Variational::Problem<LinearSystem, Trial, Test>;
      using LinearSolver = Solver::SparseLU<LinearSystem>;
      using Newton = Solver::NewtonSolver<LinearSolver>;
      using InternalWork = Rodin::Solid::InternalVirtualWork<Material, VectorField>;

      static constexpr Real NewmarkBeta = 0.25;
      static constexpr Real NewmarkGamma = 0.5;
      // 1 kPa s²/mm² = 10⁹ kg/m³.
      static constexpr Real DensityConversion = 1e-9;
      static constexpr size_t QuadratureOrder = 4;

      /**
       * @brief Constructs the passive law and prescribed active fiber branch.
       * @param cfg Solid material parameters.
       * @returns Configured constitutive law.
       */
      static Material makeMaterial(const Parameters& cfg);

      Parameters m_parameters;
      Geometry::LocalMesh m_mesh;
      ProblemData m_data;
      VectorSpace m_space;
      VectorField m_displacement, m_velocity, m_acceleration;
      VectorField m_predictedDisplacement, m_predictedVelocity;
      Real m_timeStep, m_density, m_accelerationCoefficient, m_massCoefficient;
      Real m_time = 0;
      size_t m_step = 0;
      Material m_law;
      Trial m_increment;
      Test m_test;
      InternalWork m_work;
      MechanicsProblem m_mechanics;
      LinearSolver m_linearSolver;
      Newton m_newton;
      std::unique_ptr<Output> m_output;
      bool m_initialized = false, m_finished = false, m_failed = false;
  };
}

#endif
