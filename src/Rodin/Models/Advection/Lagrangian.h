#ifndef RODIN_MODELS_ADVECTION_LAGRANGIAN_H
#define RODIN_MODELS_ADVECTION_LAGRANGIAN_H

#include <functional>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Models::Advection
{
  /**
   * @brief Lagrangian variational advection for scalar fields.
   */
  template <class Solution, class VectorField, class ... Params>
  class Lagrangian;

  template <class Solution, class VectorField, class RungeKutta>
  class Flow : public Variational::FunctionBase<Flow<Solution, VectorField, RungeKutta>>
  {
    public:
      using SolutionType = Solution;
      using VectorFieldType = VectorField;
      using RungeKuttaType = RungeKutta;


      Flow(
          Real t,
          const SolutionType& u,
          const VectorFieldType& velocity,
          const RungeKutta& rk)
        : m_t(t),
          m_solution(u),
          m_velocity(velocity),
          m_rk(rk)
      {}

      constexpr
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        const auto bt = this->forward(p);
        if (bt)
          return m_solution.get()(this->forward(p));
        else
          return 0;
      }

      std::optional<Geometry::Point> forward(const Geometry::Point& p) const
      {
        static thread_local Index s_cellIdx;
        static thread_local Math::SpatialPoint s_rc0{{}};
        static thread_local Math::SpatialPoint s_pc0{{}};
        static thread_local Math::SpatialPoint s_rc1{{}};
        static thread_local Math::SpatialPoint s_pc1{{}};

        const auto& polytope0 = p.getPolytope();
        const auto& mesh = polytope0.getMesh();
        const size_t d = mesh.getDimension();
        const auto& conn = mesh.getConnectivity();
        const auto& rk = m_rk;

        s_cellIdx = polytope0.getIndex();
        s_pc0 = p.getPhysicalCoordinates();
        s_rc0 = p.getReferenceCoordinates();

        Real tau = m_t;
        while (tau > 0)
        {
          const auto it = mesh.getPolytope(d, s_cellIdx);
          const auto& polytope = *it;
          const Geometry::Point q(polytope, s_rc0, s_pc0);
          const Geometry::Polytope::Type g = mesh.getGeometry(d, s_cellIdx);
          const auto& faces = conn.getIncidence(d, d - 1).at(s_cellIdx);
          const auto& hs = Geometry::Polytope::Traits(g).getHalfSpace();

          decltype(auto) vr = q.getJacobianInverse() * m_velocity(q);

          Real exitTime = std::numeric_limits<Real>::infinity();
          Index face;
          Index local;

          assert(faces.size() > 0);
          for (size_t i = 0; i < faces.size(); i++)
          {
            decltype(auto) normal = hs.matrix.row(i);
            const Real dot = normal.dot(vr);
            if (dot > 0)
            {
              const Real tf = (hs.vector[i] - normal.dot(s_rc0)) / dot;
              assert(tf >= 0);
              if (tf < exitTime)
              {
                local = i;
                face = faces[i];
                exitTime = tf;
              }
            }
          }
          assert(exitTime >= 0);
          assert(std::isfinite(exitTime));
          if (tau < exitTime)
          {
            rk.step(s_rc1, tau, s_rc0,
                [&](const Math::SpatialPoint& rc)
                {
                  const Geometry::Point q(polytope, rc);
                  return q.getJacobianInverse() * m_velocity(q);
                });

            s_rc0 = s_rc1;

            tau = 0;
            break;
          }
          else
          {
            if (mesh.isBoundary(face))
            {
              // Boundary reached
              return {};
            }
            else
            {
              rk.step(s_rc1, exitTime, s_rc0,
                  [&](const Math::SpatialPoint& rc)
                  {
                    const Geometry::Point q(polytope, rc);
                    return q.getJacobianInverse() * m_velocity(q);
                  });
              const auto& cells = conn.getIncidence(d - 1, d).at(face);
              assert(cells.size() == 2);
              const Index next = (cells[0] == s_cellIdx) ? cells[1] : cells[0];
              Geometry::Polytope::Project(g).face(local, s_rc1, s_rc1);
              mesh.getPolytopeTransformation(d, s_cellIdx).transform(s_pc1, s_rc1);
              mesh.getPolytopeTransformation(d, next).inverse(s_rc1, s_pc1);
              Geometry::Polytope::Project(mesh.getGeometry(d, next)).cell(s_rc1, s_rc1);

              s_rc0 = s_rc1;
              s_pc0 = s_pc1;
              s_cellIdx = next;

              tau -= exitTime;
            }
          }
        }

        return Geometry::Point(*mesh.getPolytope(d, s_cellIdx), s_rc0);
      }

    private:
      const Real m_t;
      std::reference_wrapper<const SolutionType> m_solution;
      std::reference_wrapper<const VectorFieldType> m_velocity;
      std::reference_wrapper<RungeKutta> m_rk;
  };

  template <class FES, class Data, class VectorField, class RungeKutta>
  class Lagrangian<Variational::GridFunction<FES, Data>, VectorField, RungeKutta>
  {
    public:
      using FESType = FES;
      using DataType = Data;
      using VectorFieldType = VectorField;
      using SolutionType = Variational::GridFunction<FES, Data>;
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      using SolverType = typename Data::SolverType;

      template <class Velocity, class RK>
      Lagrangian(const SolutionType& u, Velocity&& velocity, RK&& rk)
        : m_solution(u),
          m_velocity(std::forward<Velocity>(velocity)),
          m_rk(std::forward<RK>(rk)),
          m_pullback(m_solution)
      {}

      void step(const Real& dt)
      {
        const auto& fes = m_solution.get().getFiniteElementSpace();
        Variational::TrialFunction u(fes);
        Variational::TestFunction v(fes);
        Variational::Problem pb(u, v);
        pb = Integral(u, v) - Integral(m_solution.get(), Flow(v));
      }

    private:
      std::reference_wrapper<SolutionType> m_solution;
      VectorFieldType m_velocity;
      Flow<SolutionType, VectorFieldType, RungeKutta> m_pullback;
      RungeKutta m_rk;
  };
}

#endif
