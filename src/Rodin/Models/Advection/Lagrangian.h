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
 * This file provides the Lagrangian class and FirstOrderBoundaryPolicy,
 * which implement semi-Lagrangian advection schemes for scalar fields
 * in variational form.
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

namespace Rodin::Models::Advection
{
  /**
   * @brief First-order boundary policy for Lagrangian advection.
   *
   * This class implements a boundary handling policy for semi-Lagrangian
   * advection, ensuring that characteristics (particle trajectories) properly
   * interact with domain boundaries.
   *
   * @tparam Velocity Type of the velocity field
   *
   * ## Mathematical Background
   * In semi-Lagrangian methods, we trace characteristics backwards in time:
   * @f[
   *   \frac{d\mathbf{X}}{dt} = -\mathbf{v}(\mathbf{X}, t)
   * @f]
   * This policy handles cases where characteristics intersect boundaries.
   */
  template <class Velocity>
  class FirstOrderBoundaryPolicy
  {
    public:
      /**
       * @brief Constructs a boundary policy.
       *
       * @param[in] dt Time step size
       * @param[in] mesh Computational mesh
       * @param[in] velocity Velocity field
       * @param[in] eps Tolerance parameter for boundary treatment (default: 1e-12)
       */
      FirstOrderBoundaryPolicy(
          Real dt,
          const Geometry::Mesh<Context::Local>& mesh,
          const Velocity& velocity,
          Real eps = 1e-12)
        : m_dt(dt),
          m_mesh(mesh),
          m_vel(velocity),
          m_eps(eps)
      {}

      /**
       * @brief Applies boundary policy at a characteristic footpoint.
       *
       * Adjusts the characteristic footpoint when it would exit the domain,
       * projecting it back onto the boundary with appropriate handling.
       *
       * @param[in,out] tau Remaining time to trace
       * @param[in,out] cell Current cell index
       * @param[in,out] rref Reference coordinates in the cell
       * @return True if boundary treatment was successful, false otherwise
       */
      bool operator()(Real& tau, Index& cell, Math::SpatialPoint& rref) const
      {
        static thread_local Math::Vector<Real> s_nphys_u;
        static thread_local Math::SpatialPoint s_rtmp;
        static thread_local Math::SpatialPoint s_xint;

        const Real trace_sign = std::abs(m_dt);

        if (tau <= 0)
          return true;

        const auto& mesh = m_mesh.get();
        const size_t cd = mesh.getDimension();
        const auto g = mesh.getGeometry(cd, cell);
        const auto& hs = Geometry::Polytope::Traits(g).getHalfSpace();

        size_t jbest = 0;
        Real gbest = std::numeric_limits<Real>::infinity();
        for (size_t j = 0; j < static_cast<size_t>(hs.matrix.rows()); ++j)
        {
          const auto n = hs.matrix.row(j).transpose();
          const Real b = hs.vector[j];
          const Real gj = b - n.dot(rref);
          if (gj < gbest)
          {
            gbest = gj;
            jbest = j;
          }
        }

        const auto nref = hs.matrix.row(jbest).transpose();
        const auto itc = mesh.getPolytope(cd, cell);
        const auto& cellO = *itc;
        const Geometry::Point qface(cellO, rref);

        const auto& Jinv = qface.getJacobianInverse();
        s_nphys_u = Jinv.transpose() * nref;
        const Real nlen = s_nphys_u.norm();
        if (nlen == 0)
        {
          tau = 0;
          return true;
        }

        const auto nphys = trace_sign * s_nphys_u / nlen;
        decltype(auto) vphys = m_vel(qface);
        const Real vn = vphys.dot(nphys);
        const Real h = std::max<Real>(0, vn) * tau;
        const auto& xface = qface.getPhysicalCoordinates();
        s_xint = xface + (-h - m_eps) * nphys;

        mesh.getPolytopeTransformation(cd, cell).inverse(s_rtmp, s_xint);
        Geometry::Polytope::Project(g).cell(s_rtmp, s_rtmp);
        rref = s_rtmp;

        const Real b = hs.vector[jbest];
        const Real gcur = b - nref.dot(rref);
        const Real ndn = nref.dot(nref);
        if (ndn > 0)
        {
          const Real target = std::max(m_eps, gcur + m_eps);
          const Real alpha = (gcur - target) / ndn;
          rref += alpha * nref;
        }

        tau = 0;
        return true;
      }

    private:
      const Real m_dt;   ///< Time step size
      std::reference_wrapper<const Geometry::Mesh<Context::Local>> m_mesh;  ///< Reference to mesh
      std::reference_wrapper<const Velocity> m_vel;  ///< Reference to velocity field
      const Real m_eps;  ///< Tolerance for boundary treatment
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

        const FirstOrderBoundaryPolicy bp(-dt, mesh, m_velocity);
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
