/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Lagrangian.h
 * @brief Lagrangian variational advection for scalar fields.
 *
 * This file provides the Lagrangian class, which implements semi-Lagrangian
 * advection schemes for scalar fields in variational form.
 */
#ifndef RODIN_MODELS_ADVECTION_LAGRANGIAN_H
#define RODIN_MODELS_ADVECTION_LAGRANGIAN_H

#include "Rodin/Variational/Flow.h"

#include "Rodin/Math/RungeKutta/RK4.h"
#include "Rodin/Solver/CG.h"

#include "Rodin/Variational/BoundaryNormal.h"
#include "Rodin/Variational/BoundaryIntegral.h"
#include "Rodin/Variational/ProblemBody.h"
#include <functional>

namespace Rodin::Advection
{
  /**
   * @brief Boundary policy for characteristic tracing: stop-at-boundary with inward offset.
   *
   * This policy is intended to be used with semi-Lagrangian / characteristic-based
   * transport where the footpoint of a characteristic must remain inside the computational
   * domain (no sampling outside).
   *
   * When the tracing algorithm detects that the characteristic hits a boundary face
   * before consuming the full remaining time @p hit.tau, this policy:
   *  - sanitizes the reference coordinates to lie on the reported face and inside the cell,
   *  - computes an outward unit normal in *physical* coordinates at the hit point,
   *  - nudges the point slightly inside the domain along the inward normal direction,
   *  - maps back to reference coordinates and clamps into the cell,
   *  - terminates tracing by setting @p hit.tau = 0.
   *
   * In the literature this behavior is commonly described as:
   *  - **truncated characteristics**,
   *  - **stop-at-boundary**,
   *  - **clamped / clipped footpoint** (with an interior offset),
   * depending on the application.
   *
   * @tparam Velocity Type of the velocity field functor used by the tracer.
   *
   * @warning This policy does **not** enforce a physical inflow boundary condition and does
   *          not reflect characteristics. It is a *geometric* policy whose primary goal is to
   *          prevent out-of-domain sampling. If your PDE requires prescribed inflow data,
   *          you need a different policy/modeling choice.
   *
   * @par Typical use cases
   * - **Level set transport** (most common in your current setting):
   *   you want a stable, non-invasive handling of boundary hits without introducing
   *   tangential sliding artifacts or inflow modeling. Truncation + inward nudge is a good default.
   *
   * - **Semi-Lagrangian advection of a scalar** when you deliberately choose:
   *   “if the departure point would be outside, sample the closest interior point”
   *   (often used as a robust fallback in graphics/engineering codes).
   *
   * - **Geometric backtracing in ALE / mesh motion** when boundary points must remain
   *   strictly inside a cell to avoid degeneracies in interpolation/shape functions.
   *
   * @par Non-goals / when NOT to use
   * - Problems with **true inflow boundaries** where values outside must be provided from
   *   prescribed data (Dirichlet inflow), ghost-cell extension, or an extension PDE.
   * - Cases where you explicitly want **reflection / bounce** or **sliding along the boundary**
   *   (tangential continuation).
   *
   * @par Assumptions
   * - The geometric mapping is sufficiently regular for using @f$ J^{-T} n_{\mathrm{ref}} @f$
   *   as a physical normal direction on the face. (Best behavior for affine/P1 geometry.)
   *
   * @par Robustness notes
   * - The inward nudge is done in *physical space* so that the offset represents a true
   *   physical distance, independent of the reference parametrization.
   * - The policy orients the normal using the cell centroid as an interior probe, which avoids
   *   relying on a pre-defined orientation convention of face normals.
   *
   * @par Parameters
   * - @p eps_in_phys: desired inward physical offset when the boundary is hit. The effective
   *   offset is @c max(eps_in_phys, 50*sqrt(machine_epsilon)).
   *
   * @par Complexity
   * Constant time per boundary hit: one projection to cell + face, a few transforms, and a clamp.
   */
  template <class Velocity>
  class StopInsideBoundaryPolicy
  {
    public:
      /**
       * @brief Constructs the boundary policy.
       *
       * @param[in] dt          Signed time step used by the tracer (kept for API symmetry).
       * @param[in] mesh        Mesh on which tracing occurs.
       * @param[in] velocity    Velocity field functor (kept for API symmetry; not used by this policy).
       * @param[in] eps_in_phys Inward physical offset used to place the footpoint strictly inside.
       */
      StopInsideBoundaryPolicy(
          Real dt,
          const Geometry::Mesh<Context::Local>& mesh,
          const Velocity& velocity,
          Real eps_in_phys = Real(1e-12))
        : m_dt(dt),
          m_mesh(mesh),
          m_vel(velocity),
          m_eps_in_phys(eps_in_phys)
      {}

      /**
       * @brief Handles a boundary hit during characteristic tracing.
       *
       * The tracer calls this when it detects that the characteristic would hit a boundary face
       * before consuming all remaining time @p hit.tau.
       *
       * On return:
       * - @p hit.rref is updated to a point strictly inside the domain (same cell),
       * - @p hit.tau is set to 0 to stop further tracing for this step.
       *
       * @param[in,out] hit Boundary hit record (references owned by the tracer).
       * @return Always returns @c true to indicate boundary handling succeeded.
       *
       * @note If the reported face index is invalid or the physical normal degenerates,
       *       the policy conservatively stops tracing and leaves a clamped point.
       */
      bool operator()(const Variational::BoundaryHit& hit) const
      {
        if (!(hit.tau > Real(0)))
          return true;

        const auto& mesh = m_mesh.get();
        const size_t cd  = mesh.getDimension();
        const Index  c   = hit.cell;

        const auto g   = mesh.getGeometry(cd, c);
        const auto& hs = Geometry::Polytope::Traits(g).getHalfSpace();

        const size_t j = hit.face;
        if (j >= static_cast<size_t>(hs.vector.size()))
        {
          hit.tau = 0;
          return true;
        }

        // -------- (1) sanitize reference point: inside + on face --------
        Math::SpatialPoint r = hit.rref;
        Math::SpatialPoint rtmp;
        Geometry::Polytope::Project(g).cell(rtmp, r);
        Geometry::Polytope::Project(g).face(j, r, rtmp);

        // -------- (2) physical boundary point --------
        Math::SpatialPoint x;
        mesh.getPolytopeTransformation(cd, c).transform(x, r);

        // -------- (3) compute outward unit normal in physical space --------
        const auto nref = hs.matrix.row(j).transpose();

        const auto itc  = mesh.getPolytope(cd, c);
        const auto& cell = *itc;

        Geometry::Point qface(cell, r, x);
        const auto JinvT = qface.getJacobianInverse().transpose();

        Math::SpatialVector<Real> nu = JinvT * nref; // unnormalized physical normal
        Real nn = nu.dot(nu);
        if (!(nn > Real(0)) || !std::isfinite(nn))
        {
          hit.rref = r;
          hit.tau  = 0;
          return true;
        }

        Math::SpatialVector<Real> nhat = nu / std::sqrt(nn);

        // Orient nhat so the cell centroid is on the "interior" side:
        // interior condition: (nhat·x_face - nhat·x_centroid) >= 0
        const auto rcent = Geometry::Polytope::Traits(g).getCentroid();
        Math::SpatialPoint xc;
        mesh.getPolytopeTransformation(cd, c).transform(xc, rcent);

        const Real cplane = nhat.dot(x);
        const Real sd_in  = cplane - nhat.dot(xc);
        if (sd_in < Real(0))
          nhat = -nhat;

        // -------- (4) nudge inside in physical space (inside = -nhat) --------
        const Real eps_machine = std::numeric_limits<Real>::epsilon();
        const Real eps_floor   = Real(50) * std::sqrt(eps_machine);
        const Real eps_in      = std::max<Real>(m_eps_in_phys, eps_floor);

        x -= eps_in * nhat;

        // -------- (5) map back + clamp --------
        mesh.getPolytopeTransformation(cd, c).inverse(rtmp, x);
        Geometry::Polytope::Project(g).cell(r, rtmp);

        hit.rref = r;
        hit.tau  = 0;
        return true;
      }

    private:
      const Real m_dt;
      std::reference_wrapper<const Geometry::Mesh<Context::Local>> m_mesh;
      std::reference_wrapper<const Velocity> m_vel; // kept for symmetry; not used here
      const Real m_eps_in_phys;
  };

  /**
   * @brief Lagrangian variational advection for scalar fields.
   *
   * This class implements semi-Lagrangian time-stepping for the advection
   * equation:
   * @f[
   *   \frac{\partial u}{\partial t} + \mathbf{v} \cdot \nabla u = 0
   * @f]
   * where @f$ u @f$ is the scalar field and @f$ \mathbf{v} @f$ is the velocity
   * field.
   *
   * ## Semi-Lagrangian Method
   * At each time step, the method:
   * 1. Traces characteristics backwards: @f$ \mathbf{X}(t) = \mathbf{x} - \int_0^{\Delta t} \mathbf{v} \, ds @f$
   * 2. Interpolates the solution at departure points
   * 3. Solves a projection problem to obtain the new solution
   *
   * The variational formulation at each step is:
   * @f[
   *   \int_\Omega u^{n+1} v \, dx = \int_\Omega u^n(\mathbf{X}) v \, dx
   * @f]
   *
   * ## Advantages
   * - Unconditionally stable for large CFL numbers
   * - Naturally handles complex geometries
   * - Mass-conservative in variational form
   *
   * @tparam Params Parameter pack for class specialization
   */
  template <class ... Params>
  class Lagrangian;

  /**
   * @brief Lagrangian advection specialization for trial/test function formulation.
   *
   * @tparam FES Finite element space type
   * @tparam Data Data storage type for grid function
   * @tparam Initial Type of initial condition
   * @tparam VectorField Type of velocity field
   * @tparam Step Type of time-stepping scheme (default: RK4)
   */
  template <class FES, class Data, class Initial, class VectorField, class Step>
  class Lagrangian<
    Variational::TrialFunction<Variational::GridFunction<FES, Data>, FES>,
    Variational::TestFunction<FES>, Initial, VectorField, Step>
  {
    public:
      /// Finite element space type
      using FESType =
        FES;

      /// Data storage type for grid function
      using DataType =
        Data;

      /// Type of the initial condition
      using InitialType =
        Initial;

      /// Type of the velocity field
      using VectorFieldType =
        VectorField;

      /// Type of the time-stepping scheme
      using StepType =
        Step;

      /// Solution type: grid function holding the advected field
      using SolutionType =
        Variational::GridFunction<FES, Data>;

      /// Trial function type for variational formulation
      using TrialFunctionType =
        Variational::TrialFunction<Variational::GridFunction<FES, Data>, FES>;

      /// Test function type for variational formulation
      using TestFunctionType =
        Variational::TestFunction<FES>;

      /**
       * @brief Constructs a Lagrangian advection solver.
       *
       * @tparam U0 Type of initial condition (deduced)
       * @tparam VVel Type of velocity field (deduced)
       * @tparam S Type of time-stepping scheme (deduced, default: RK4)
       * @param[in,out] u Trial function holding the solution
       * @param[in,out] v Test function for variational formulation
       * @param[in] u0 Initial condition
       * @param[in] vel Velocity field
       * @param[in] st Time-stepping scheme (default: RK4)
       */
      template <class U0, class VVel, class S = StepType>
      Lagrangian(TrialFunctionType& u, TestFunctionType& v, U0&& u0, VVel&& vel, S&& st = S{})
        : m_t(0),
          m_u(u), m_v(v),
          m_initial(std::forward<U0>(u0)),
          m_velocity(std::forward<VVel>(vel)),
          m_step(std::forward<S>(st))
      {}

      /**
       * @brief Advances the solution by one time step.
       *
       * Performs a semi-Lagrangian time step:
       * 1. Computes characteristics backwards by @f$ \Delta t @f$
       * 2. Interpolates solution at departure points
       * 3. Projects onto finite element space
       *
       * @param[in] dt Time step size
       *
       * @note For the first step (@f$ t = 0 @f$), uses the initial condition.
       * For subsequent steps, uses the current solution as departure values.
       */
      void step(const Real& dt)
      {
        using namespace Variational;

        auto& u = m_u.get();
        auto& v = m_v.get();

        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();

        const StopInsideBoundaryPolicy bp(-dt, mesh, m_velocity);
        const Math::RungeKutta::RK4 step;
        const DefaultTangentPolicy tp;

        Problem pb(u, v);
        if (m_t > 0)
        {
          pb = Integral(u, v)
             - Integral(Flow(-dt, u.getSolution(), m_velocity, step, bp, tp), v);
        }
        else
        {
          pb = Integral(u, v)
             - Integral(Flow(-dt, m_initial, m_velocity, step, bp, tp), v);
        }

        Solver::CG(pb).solve();

        m_t += dt;
      }

    private:
      Real m_t;  ///< Current simulation time

      std::reference_wrapper<TrialFunctionType> m_u;  ///< Reference to trial function
      std::reference_wrapper<TestFunctionType> m_v;   ///< Reference to test function

      InitialType m_initial;        ///< Initial condition
      VectorFieldType m_velocity;   ///< Velocity field
      StepType m_step;              ///< Time-stepping scheme
  };

  /**
   * @brief Deduction guide for Lagrangian with default RK4 stepper.
   *
   * Allows construction without explicitly specifying the Step template parameter.
   */
  template <class FES, class Data, class Initial, class VVel>
  Lagrangian(Variational::TrialFunction<Variational::GridFunction<FES, Data>, FES>&,
             Variational::TestFunction<FES>&,
             Initial&&,
             VVel&&)
  -> Lagrangian<
       Variational::TrialFunction<Variational::GridFunction<FES,Data>,FES>,
       Variational::TestFunction<FES>,
       Initial,
       VVel,
       Math::RungeKutta::RK4>;

  /**
   * @brief Deduction guide for Lagrangian with custom stepper.
   *
   * Allows construction with explicit time-stepping scheme specification.
   */
  template <class FES, class Data, class Initial, class VVel, class SStep>
  Lagrangian(Variational::TrialFunction<Variational::GridFunction<FES, Data>,FES>&,
             Variational::TestFunction<FES>&,
             Initial&&,
             VVel&&,
             SStep&&)
  -> Lagrangian<
       Variational::TrialFunction<Variational::GridFunction<FES, Data>, FES>,
       Variational::TestFunction<FES>,
       Initial,
       VVel,
       SStep>;
}

#endif
