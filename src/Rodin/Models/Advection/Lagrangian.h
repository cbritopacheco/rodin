#ifndef RODIN_MODELS_ADVECTION_LAGRANGIAN_H
#define RODIN_MODELS_ADVECTION_LAGRANGIAN_H

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/FaceNormal.h"
#include "Rodin/Variational/GridFunction.h"

namespace Rodin::Models::Advection
{
  /**
   * @brief Lagrangian variational advection for scalar fields.
   *
   * The method solves the advection equation:
   * @f[
   * \frac{\partial u}{\partial t} + \beta \cdot \nabla u = 0
   * @f]
   * assuming that the velocity field @f$ \beta @f$ is divergence-free.
   */
  template <class Solution, class VectorField>
  class Lagrangian;

  template <class FES, class Data, class VectorField>
  class Lagrangian<Variational::GridFunction<FES, Data>, VectorField>
  {
    public:
      using FESType = FES;
      using DataType = Data;
      using VectorFieldType = VectorField;
      using SolutionType = Variational::GridFunction<FES, Data>;
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      using SolverType = typename Data::SolverType;

      struct Trace
      {
        Geometry::Point point;
        ScalarType integral;
      };

      template <class RungeKutta>
      class Pullback : public Variational::FunctionBase<Pullback<RungeKutta>>
      {
        public:
          using RungeKuttaType = RungeKutta;

          template <class Velocity, class RK>
          Pullback(const SolutionType& u, Velocity&& velocity,  RK&& rk)
            : m_solution(u),
              m_velocity(std::forward<Velocity>(velocity)),
              m_rk(std::forward<RK>(rk))
          {}

          constexpr
          auto getValue(const Geometry::Point& p) const
          {
            Trace trace = backtrace(p);
            return m_solution.get()(trace.point);
          }

          Trace backtrace(const Real& dt, const Geometry::Point& p) const
          {
            static thread_local Math::SpatialPoint s_rc{{}};
            static thread_local Math::SpatialPoint s_pc{{}};

            struct Exit
            {
              Real t;
              Index face;
            };

            const auto& rk = m_rk;
            const auto& polytope = p.getPolytope();
            const size_t d = polytope.getDimension();
            const Index idx = polytope.getIndex();
            const auto& mesh = polytope.getMesh();
            const auto& conn = mesh.getConnectivity();

            // Reference velocity
            const auto vr = p.getJacobianInverse() * m_velocity(p);

            // Reference point
            s_pc = p.getPhysicalCoordinates();
            s_rc = p.getReferenceCoordinates();

            Real tau = dt;
            while (tau > 0)
            {
              const auto it = mesh.getPolytope(d, idx);
              const Geometry::Polytope::Type g = it->getGeometry();

              const auto& hs = Geometry::Polytope::Traits(g).getHalfSpace();
              Real exitTime = std::numeric_limits<Real>::infinity();
              Index face;
              Index local = 0;
              for (const Index& f : conn.getIncidence(d, d - 1).at(idx))
              {
                const auto& normal = hs.matrix.row(local);
                const Real& b = hs.vector[local];
                const Real tf = (b - normal.dot(s_rc)) / normal.dot(vr);
                assert(tf >= 0);
                if (tf < exitTime)
                {
                  exitTime = tf;
                  face = f;
                }
              }
              assert(std::isfinite(exitTime));

              if (tau < exitTime)
              {
                rk.step(tau, s_rc);
                tau = 0;
                break;
              }
              else
              {
                rk.step(exit, s_rc);
                const auto& inc = conn.getIncidence(d - 1, d).at(face);

                // Snap to the face

                assert(inc.size() == 2);
                for (const Index& cell : inc)
                {
                  if (cell != idx)
                  {
                    mesh.getPolytopeTransformation(d, cell).inverse(s_rc, s_pc);
                    break;
                  }
                }
                tau -= exitTime;
              }
            }
          }

          void abcissa()
          {
          }

        private:
          std::reference_wrapper<const SolutionType> m_solution;
          VectorFieldType m_velocity;
          RungeKutta m_rk;
          Real m_guard;
      };

      template <class Velocity>
      Lagrangian(const SolutionType& u, const Velocity& theta)
        : m_solution(u),
          m_velocity(theta),
          m_pullback(m_solution)
      {}

      void step(const Real& dt)
      {
        const auto& fes = m_solution.get().getFiniteElementSpace();
        Variational::TrialFunction u(fes);
        Variational::TestFunction v(fes);
        Variational::Problem pb(u, v);
        pb = Integral(u, v) - Integral(m_pullback, v);
      }

    private:
      std::reference_wrapper<SolutionType> m_solution;
      VectorFieldType m_velocity;
      Pullback m_pullback;
  };
}

#endif
