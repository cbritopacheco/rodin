#ifndef RODIN_MODELS_ADVECTION_LAGRANGIAN_H
#define RODIN_MODELS_ADVECTION_LAGRANGIAN_H

#include <functional>

#include "Rodin/Math/RungeKutta/RK4.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/ShapeFunction.h"

namespace Rodin::Models::Advection
{
  /**
   * @brief Lagrangian variational advection for scalar fields.
   */
  template <class Solution, class VectorField, class ... Params>
  class Lagrangian;

  template <class Operand, class VectorField, class ... Params>
  class Flow;

  template <
    class Derived,
    class FES, Variational::ShapeFunctionSpaceType Space,
    class VectorField, class RungeKutta>
  class Flow<Variational::ShapeFunctionBase<Derived, FES, Space>, VectorField, RungeKutta>
  : public Variational::ShapeFunctionBase<
      Flow<Variational::ShapeFunctionBase<Derived, FES, Space>, VectorField, RungeKutta>,
      FES, Space>
  {
    public:
      using Operand = Variational::ShapeFunctionBase<Derived, FES, Space>;
      using VectorFieldType = VectorField;
      using RungeKuttaType = RungeKutta;
      using Parent = Variational::ShapeFunctionBase<
        Flow<Operand, VectorFieldType, RungeKuttaType>, FES, Space>;

      template <class RK>
      Flow(Real t, const Operand& u, const VectorFieldType& velocity, RK&& rk = Math::RungeKutta::RK4())
        : m_t(t),
          m_operand(u.copy()),
          m_velocity(velocity),
          m_rk(std::forward<RK>(rk))
      {}

      Flow(const Flow& other)
        : Parent(other),
          m_t(other.m_t),
          m_operand(other.m_operand->copy()),
          m_velocity(other.m_velocity),
          m_rk(other.m_rk),
          m_p(other.m_p)
      {}

      Flow(Flow&& other)
        : Parent(std::move(other)),
          m_t(std::exchange(other.m_t, 0)),
          m_operand(std::move(other.m_operand)),
          m_velocity(std::move(other.m_velocity)),
          m_rk(std::move(other.m_rk)),
          m_p(std::exchange(other.m_p, std::nullopt))
      {}

      std::optional<Geometry::Point> forward(const Geometry::Point& p) const
      {
        static thread_local Index s_cellIdx;
        static thread_local Math::SpatialPoint s_rc0{{}};
        static thread_local Math::SpatialPoint s_pc0{{}};
        static thread_local Math::SpatialPoint s_rc1{{}};
        static thread_local Math::SpatialPoint s_pc1{{}};

        const auto& polytope0 = p.getPolytope();
        const auto& mesh = polytope0.getMesh();
        const size_t cellDim = mesh.getDimension();
        const auto& conn = mesh.getConnectivity();

        s_pc0 = p.getPhysicalCoordinates();

        if (polytope0.getDimension() == cellDim)
        {
          s_cellIdx = polytope0.getIndex();
          s_rc0 = p.getReferenceCoordinates();
        }
        else if (polytope0.getDimension() == cellDim - 1)
        {
          const Index fidx = polytope0.getIndex();
          const auto& adj = conn.getIncidence(cellDim - 1, cellDim).at(fidx);

          if (mesh.isBoundary(fidx))
          {
            // boundary face → unique adjacent cell
            assert(adj.size() == 1);
            const Index c = adj[0];
            mesh.getPolytopeTransformation(cellDim, c).inverse(s_rc1, s_pc0);
            Geometry::Polytope::Project(mesh.getGeometry(cellDim, c)).cell(s_rc0, s_rc1);
            s_cellIdx = c;
          }
          else
          {
            // interior face → pick side via n_ref·a in adj[0]
            assert(adj.size() == 2);

            const Index c0 = adj[0];

            const auto& faces0 = conn.getIncidence(cellDim, cellDim - 1).at(c0);
            size_t j0 = faces0.size();
            for (size_t k = 0; k < faces0.size(); ++k)
            {
              if (faces0[k] == fidx) { j0 = k; break; }
            }
            assert(j0 < faces0.size());

            mesh.getPolytopeTransformation(cellDim, c0).inverse(s_rc1, s_pc0);
            const auto& cell0 = *mesh.getPolytope(cellDim, c0);
            const Geometry::Point q0(cell0, s_rc1, s_pc0);

            const auto a0 = q0.getJacobianInverse() * m_velocity(q0);

            const auto& hs0 =
              Geometry::Polytope::Traits(mesh.getGeometry(cellDim, c0)).getHalfSpace();
            const auto nref0 = hs0.matrix.row(j0);

            if (nref0.dot(a0) > 0)
            {
              const Index c1 = adj[1];
              mesh.getPolytopeTransformation(cellDim, c1).inverse(s_rc1, s_pc0);
              Geometry::Polytope::Project(mesh.getGeometry(cellDim, c1)).cell(s_rc0, s_rc1);
              s_cellIdx = c1;
            }
            else
            {
              Geometry::Polytope::Project(mesh.getGeometry(cellDim, c0)).cell(s_rc0, s_rc1);
              s_cellIdx = c0;
            }
          }
        }

        Real tau = m_t;
        while (tau > 0)
        {
          const auto it = mesh.getPolytope(cellDim, s_cellIdx);
          const auto& polytope = *it;
          const Geometry::Point q(polytope, s_rc0, s_pc0);

          const auto vr = q.getJacobianInverse() * m_velocity(q);

          const Geometry::Polytope::Type g = mesh.getGeometry(cellDim, s_cellIdx);
          const auto& faces = conn.getIncidence(cellDim, cellDim - 1).at(s_cellIdx);
          const auto& hs = Geometry::Polytope::Traits(g).getHalfSpace();

          Real exitTime = std::numeric_limits<Real>::infinity();
          size_t local = faces.size();

          for (size_t i = 0; i < faces.size(); i++)
          {
            const auto nref = hs.matrix.row(i);
            const Real denom = nref.dot(vr);
            if (denom > 0)
            {
              const Real numer = hs.vector[i] - nref.dot(s_rc0);
              const Real tf = numer / denom; // may be negative; admissibility handled by denom>0
              if (tf < exitTime)
              {
                exitTime = tf;
                local = i;
              }
            }
          }

          if (local == faces.size())
          {
            // no admissible face (e.g., tangent); stop
            return Geometry::Point(*mesh.getPolytope(cellDim, s_cellIdx), s_rc0);
          }

          const Index face = faces[local];

          if (tau < exitTime)
          {
            m_rk.step(
              s_rc1, tau, s_rc0,
              [&](const Math::SpatialPoint& rc)
              {
                const Geometry::Point qp(polytope, rc);
                return qp.getJacobianInverse() * m_velocity(qp);
              });

            s_rc0 = s_rc1;
            tau = 0;
            break;
          }
          else
          {
            if (mesh.isBoundary(face))
            {
              // push-forward: outflow reached → drop
              return {};
            }

            m_rk.step(
              s_rc1, exitTime, s_rc0,
              [&](const Math::SpatialPoint& rc)
              {
                const Geometry::Point qp(polytope, rc);
                return qp.getJacobianInverse() * m_velocity(qp);
              });

            Geometry::Polytope::Project(g).face(local, s_rc1, s_rc1);

            mesh.getPolytopeTransformation(cellDim, s_cellIdx).transform(s_pc1, s_rc1);

            const auto& cells = conn.getIncidence(cellDim - 1, cellDim).at(face);
            assert(cells.size() == 2);
            const Index next = (cells[0] == s_cellIdx) ? cells[1] : cells[0];

            mesh.getPolytopeTransformation(cellDim, next).inverse(s_rc1, s_pc1);
            Geometry::Polytope::Project(mesh.getGeometry(cellDim, next)).cell(s_rc1, s_rc1);

            s_rc0 = s_rc1;
            s_pc0 = s_pc1;
            s_cellIdx = next;

            tau -= exitTime;
          }
        }

        return Geometry::Point(*mesh.getPolytope(cellDim, s_cellIdx), s_rc0);
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return m_operand->getDOFs(polytope);
      }

      constexpr
      decltype(auto) getBasis(size_t local) const
      {
        return m_operand->getBasis(local);
      }

      const Geometry::Point& getPoint() const
      {
        assert(m_p);
        return *m_p;
      }

      Flow& setPoint(const Geometry::Point& p)
      {
        m_p = this->forward(p);
        return *this;
      }

      constexpr
      const auto& getLeaf() const
      {
        return m_operand->getLeaf();
      }

      virtual Flow* copy() const noexcept override
      {
        return new Flow(*this);
      }

    private:
      const Real m_t;
      std::unique_ptr<Operand> m_operand;
      std::reference_wrapper<const VectorFieldType> m_velocity;
      RungeKutta m_rk;
      std::optional<Geometry::Point> m_p;
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
