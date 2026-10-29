/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_FSI_H
#define RODIN_ATRIUM_FSI_H

#include "Solid/Solid.h"
#include "Fluid.h"

namespace Rodin::Examples::Heart::Atrium
{
  /**
   * @brief Partitioned Robin--Dirichlet coupling on a conforming MPI interface.
   *
   * The reference-area Robin equation is
   * @f$ P_sN_s+\alpha v_s=g @f$, with the fluid update
   * @f$ g=J_\Gamma(-\sigma_f n_f)+\alpha u_f @f$.
   * Velocity continuity cancels the Robin terms at convergence. The load is
   * obtained from an L2 projection of pressure and recovered-gradient stress.
   *
   * ## Architecture
   *
   * Parent vertex identities define retained PETSc scatter plans. Solid motion
   * is harmonically extended to the fluid mesh; an interface mass problem
   * projects fluid traction back to the reference wall. Coupling trials retain
   * the previous physical state and are committed only after convergence.
   */
  class FSI
  {
    public:
      struct Parameters
      {
        size_t maximumIterations =
          100; ///< Coupling trial budget; failure does not commit.
        Real relativeTolerance = 1e-5;
        Real absoluteTolerance = 1e-4; ///< Robin residual tolerance in Pa.
        Real robinCoefficient = 0; ///< Pa s/m; zero selects rho V^(1/3)/dt.
        Real initialRelaxation = .5; ///< Initial Aitken relaxation factor.
        Real minimumRelaxation = .05; ///< Heuristic lower bound on Aitken relaxation.
        Real maximumRelaxation = 1; ///< Upper bound preventing over-relaxation.
      };

      struct Diagnostics
      {
        size_t iterations = 0;
        Real robinResidual = 0; ///< Pa, maximum nodal component difference.
        Real minimumMeshJacobian = 1;
        Real wallFlux = 0; ///< m^3/s, outward from the fluid.
        Real robinTolerance = 0; ///< Pa, stopping threshold for the accepted trial.
        Real velocityMismatch = 0; ///< m/s, maximum interface velocity mismatch.
        Real displacementMismatch = 0; ///< m, maximum interface geometry mismatch.
        Real relaxation = 1;
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

      /// @brief Advances and commits one converged coupled timestep.
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
       * @brief Builds the material-4 reference interface.
       * @param fluid Fluid submesh with boundary connectivity.
       * @returns Reference interface submesh with parent ancestry.
       */
      static Fluid::MeshType makeInterface(const Fluid& fluid);

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

      /// @brief Projects the fluid reaction per reference interface area.
      void updateTraction();

      Solid& m_solid;
      Fluid& m_fluid;
      Parameters m_parameters;
      Diagnostics m_diagnostics;
      Fluid::MeshType m_interface;
      Fluid::VectorFESType m_interfaceSpace;
      Fluid::VectorTrialFunctionType m_motionTrial, m_tractionTrial;
      Fluid::VectorTestFunctionType m_motionTest, m_tractionTest;
      Fluid::VectorProblemType m_motion, m_traction;
      Solver::KSP m_motionSolver, m_tractionSolver;
      Fluid::VectorGridFunctionType m_previousMeshDisplacement;
      VectorField m_importedTraction, m_importedVelocity, m_previousRobinData;
      VectorField m_previousResidual, m_residualDifference;
      std::vector<Math::SpatialPoint> m_referenceVertices;
      std::vector<Real> m_referenceCellDeterminants, m_previousCellMeasures;
      ::VecScatter m_solidToFluid = nullptr;
      ::VecScatter m_tractionToSolid = nullptr;
      Real m_robinCoefficient;
  };
}

#endif
