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

namespace Rodin::Models::Advection
{
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
          m_initial(std::forward<U0>(u0)),   // may be value or ref (e.g., U0 = T or T&)
          m_velocity(std::forward<VVel>(vel)), // T or T&
          m_step(std::forward<S>(st))
      {}

      void step(const Real& dt)
      {
        using namespace Variational;

        auto& u = m_u.get();
        auto& v = m_v.get();

        BoundaryNormal n(u.getFiniteElementSpace().getMesh());

        const RealFunction vn =
          [&](const Geometry::Point& p) -> Real
          {
            static thread_local Math::SpatialVector<Real> s_n;
            s_n = n(p);
            const auto dot = m_velocity(p).dot(s_n);
            if (dot > 0)
              return dot;
            else
              return 0;
          };

        Problem pb(u, v);
        if (m_t > 0)
        {
          pb = Integral(u, v)
             - Integral(u.getSolution(), Flow(dt, v, m_velocity, m_step))
             ;
        }
        else
        {
          pb = Integral(u, v)
             - Integral(m_initial, Flow(dt, v, m_velocity, m_step))
             ;
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
