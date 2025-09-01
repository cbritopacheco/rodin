#ifndef RODIN_MODELS_ADVECTION_LAGRANGIAN_H
#define RODIN_MODELS_ADVECTION_LAGRANGIAN_H

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Math/Vector.h"
#include <functional>

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
  template <class Solution, class VectorField, class ... Params>
  class Lagrangian;

  template <class Solution, class VectorField, class RungeKutta>
  class Pullback : public Variational::FunctionBase<Pullback<Solution, VectorField, RungeKutta>>
  {
    public:
      using SolutionType = Solution;
      using VectorFieldType = VectorField;
      using RungeKuttaType = RungeKutta;

      Pullback(const SolutionType& u, const VectorFieldType& velocity, const RungeKutta& rk)
        : m_solution(u),
          m_velocity(velocity),
          m_rk(rk)
      {}

      constexpr
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return m_solution.get()(backtrace(p));
      }

      Geometry::Point backtrace(const Real& dt, const Geometry::Point& p) const
      {
        static thread_local Math::SpatialPoint s_rc{{}};
        static thread_local Math::SpatialPoint s_pc{{}};

        const auto& rk = m_rk;
        const auto& polytope = p.getPolytope();
        const size_t d = polytope.getDimension();
        const auto& mesh = polytope.getMesh();
        const auto& conn = mesh.getConnectivity();
        Index cellIdx = polytope.getIndex();
        s_pc = p.getPhysicalCoordinates();
        s_rc = p.getReferenceCoordinates();
        Real tau = dt;
        while (tau > 0)
        {
          const Geometry::Point q(*mesh.getPolytope(d, cellIdx), s_rc, s_pc);
          const auto vr = p.getJacobianInverse() * m_velocity(p);
          const Geometry::Polytope::Type g = mesh.getGeometry(d, cellIdx);
          const auto& faces = conn.getIncidence(d, d - 1).at(cellIdx);
          const auto& hs = Geometry::Polytope::Traits(g).getHalfSpace();

          Real exitTime = std::numeric_limits<Real>::infinity();
          Index face;
          Index local;
          for (size_t i = 0; i < faces.size(); i++)
          {
            const auto& normal = hs.matrix.row(local);
            const Real dot = normal.dot(vr);
            if (dot > 0)
            {
              const Real tf = (hs.vector[local] - normal.dot(s_rc)) / dot;
              assert(tf >= 0);
              if (tf < exitTime)
              {
                local = i;
                face = faces[i];
                exitTime = tf;
              }
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

            const auto& cells = conn.getIncidence(d - 1, d).at(face);
            assert(cells.size() == 2);
            const Index next = (cells[0] == cellIdx) ? cells[1] : cells[0];
            Geometry::Polytope::Project(g).face(local, s_rc, s_rc);
            mesh.getPolytopeTransformation(d, cellIdx).transform(s_pc, s_rc);
            mesh.getPolytopeTransformation(d, next).inverse(s_rc, s_pc);
            Geometry::Polytope::Project(mesh.getGeometry(d, next)).cell(s_rc, s_rc);
            tau -= exitTime;
            cellIdx = next;
          }
        }
        return Geometry::Point(*mesh.getPolytope(d, cellIdx), s_rc);
      }

    private:
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
        pb = Integral(u, v) - Integral(m_pullback, v);
      }

    private:
      std::reference_wrapper<SolutionType> m_solution;
      VectorFieldType m_velocity;
      Pullback<SolutionType, VectorFieldType, RungeKutta> m_pullback;
      RungeKutta m_rk;
  };
}

#endif
