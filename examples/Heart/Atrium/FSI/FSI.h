/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_FSI_FSI_H
#define RODIN_ATRIUM_FSI_FSI_H

#include "../Solid/Solid.h"
#include "../Fluid/Fluid.h"

namespace Rodin::Examples::Heart::Atrium
{
  /**
   * @brief Partitioned Robin--Robin coupling on a conforming MPI interface.
   *
   * The solid condition is @f$ P_sN_s+\alpha_0 v_\Gamma=T+\alpha_0u_f @f$; the
   * fluid uses the same reference-area impedance and the lagged traction. The
   * update @f$ T^{k+1}=T^k+\alpha_0(u_f^{k+1}-v_\Gamma^{k+1}) @f$ is the fluid
   * variational reaction, without relaxation. Transmission uses the BDF1 wall rate
   * @f$ (d_s^{n+1}-d_s^n)/\Delta t @f$ on both sides. Solid dynamics retain
   * Newmark integration. Relative solid displacement change determines early
   * termination; budget exhaustion accepts the last completed sweep.
   *
   * ## Architecture
   *
   * Parent vertex identities define a retained PETSc trace scatter plan. Solid
   * displacement is harmonically extended to the fluid mesh. Matching P1 traces
   * permit the algebraic traction update. Previous physical histories remain fixed
   * throughout the sweeps. The geometry is updated for every accepted pair.
   */
  class FSI
  {
    public:
      struct Parameters
      {
        size_t maximumIterations =
          5; ///< Sweep budget; the last completed sweep is accepted.
        Real relativeTolerance =
          1e-6; ///< Relative global solid displacement-change tolerance.
        Real robinCoefficient = 0; ///< Pa s/m; zero selects sqrt(rho_s E_matrix).
      };

      struct Diagnostics
      {
        size_t iterations = 0;
        Real tractionDataUpdate = 0; ///< Pa, maximum Robin-data or traction update.
        Real minimumMeshJacobian = 1;
        Real wallFlux = 0; ///< m^3/s, outward from the fluid.
        Real displacementChange =
          0; ///< Relative L2 change, or absolute L2 change in m at zero state.
        bool converged =
          false; ///< Displacement-change criterion satisfied after at least two sweeps.
        bool budgetExhausted =
          false; ///< Last sweep accepted without satisfying the criterion.
        Real velocityMismatch = 0; ///< m/s, fluid/BDF1 interface velocity mismatch.
        Real newmarkVelocityMismatch = 0; ///< m/s, fluid/Newmark endpoint mismatch.
        Real robinCoefficient = 0; ///< Pa s/m, reference-area impedance.
        Real displacementMismatch = 0; ///< m, maximum interface geometry mismatch.
      };

      /**
       * @brief Binds components and constructs their shared reference interface.
       * @param solid Active solid component, not yet initialized.
       * @param fluid Fluid component, not yet initialized.
       * @param parameters Coupling tolerances and Robin parameter.
       */
      FSI(Solid& solid, Fluid& fluid, const Parameters& parameters);
      ~FSI();
      FSI(const FSI&) = delete;
      FSI& operator=(const FSI&) = delete;

      /// @brief Matches the initial fluid geometry to passive solid equilibrium.
      void initialize();

      /// @brief Advances and commits the last completed coupling sweep.
      void step();

      /**
       * @brief Returns diagnostics for the latest coupled step.
       * @returns Cached coupling diagnostics.
       */
      const Diagnostics& getDiagnostics() const noexcept
      {
        return m_diagnostics;
      }

    private:
      /**
       * @brief Builds a nodal scatter using shared distributed parent identities.
       * @param sourceMesh Mesh of the source field.
       * @param source Source field defining the source layout.
       * @param targetMesh Mesh of the target field.
       * @param target Target field defining the target layout.
       * @returns Owned scatter plan, destroyed by this coupling object.
       */
      template <class SourceField, class TargetField>
      static ::VecScatter makeScatter(const Fluid::MeshType& sourceMesh,
        const SourceField& source, const Fluid::MeshType& targetMesh,
        const TargetField& target);

      /**
       * @brief Transfers an interface field and refreshes its ghosts.
       * @param plan Retained interface scatter.
       * @param source Source field.
       * @param target Target field.
       * @param scale Unit conversion applied after the scatter.
       * @param direction Forward or reverse direction of the plan.
       */
      template <class SourceField, class TargetField>
      static void transfer(::VecScatter plan, const SourceField& source,
        TargetField& target, Real scale, ::ScatterMode direction = SCATTER_FORWARD);

      /// @brief Restores the reference fluid geometry before the harmonic solve.
      void restoreReference();

      /// @brief Extends solid displacement and updates the ALE geometry and velocity.
      void moveFluidMesh();

      /// @brief Updates the reference Robin traction from the interface velocity jump.
      void updateTraction();

      Solid& m_solid;
      Fluid& m_fluid;
      Parameters m_parameters;
      Diagnostics m_diagnostics;
      Fluid::VectorTrialFunctionType m_motionTrial;
      Fluid::VectorTestFunctionType m_motionTest;
      Fluid::VectorProblemType m_motion;
      Solver::KSP m_motionSolver;
      Fluid::VectorGridFunctionType m_previousMeshDisplacement;
      VectorField m_importedTraction, m_trialTraction, m_importedVelocity;
      VectorField m_interfaceVelocity, m_previousRobinData;
      VectorField m_previousDisplacementIteration, m_residualDifference;
      VectorField m_tractionResidual;
      std::vector<Math::SpatialPoint> m_referenceVertices;
      std::vector<Real> m_referenceCellDeterminants, m_previousCellMeasures;
      std::vector<Real> m_referenceFaceMeasures;
      ::VecScatter m_solidToFluid = nullptr;
      Real m_robinCoefficient;
  };
}

#endif
