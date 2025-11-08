/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
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
  template <class Velocity>
  class FirstOrderBoundaryPolicy
  {
    public:
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
      const Real m_dt;
      std::reference_wrapper<const Geometry::Mesh<Context::Local>> m_mesh;
      std::reference_wrapper<const Velocity> m_vel;
      const Real m_eps;
  };

  /**
   * @brief Lagrangian variational advection for scalar fields.
   */
  template <class ... Params>
  class Lagrangian;

  template <class FES, class Data, class Initial, class VectorField, class Step>
  class Lagrangian<
    Variational::TrialFunction<Variational::GridFunction<FES, Data>, FES>,
    Variational::TestFunction<FES>, Initial, VectorField, Step>
  {
    public:
      using FESType =
        FES;

      using DataType =
        Data;

      using InitialType =
        Initial;

      using VectorFieldType =
        VectorField;

      using StepType =
        Step;

      using SolutionType =
        Variational::GridFunction<FES, Data>;

      using TrialFunctionType =
        Variational::TrialFunction<Variational::GridFunction<FES, Data>, FES>;

      using TestFunctionType =
        Variational::TestFunction<FES>;

      template <class U0, class VVel, class S = StepType>
      Lagrangian(TrialFunctionType& u, TestFunctionType& v, U0&& u0, VVel&& vel, S&& st = S{})
        : m_t(0),
          m_u(u), m_v(v),
          m_initial(std::forward<U0>(u0)),
          m_velocity(std::forward<VVel>(vel)),
          m_step(std::forward<S>(st))
      {}

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
      Real m_t;

      std::reference_wrapper<TrialFunctionType> m_u;
      std::reference_wrapper<TestFunctionType> m_v;

      InitialType m_initial;
      VectorFieldType m_velocity;
      StepType m_step;
  };

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
