#ifndef RODIN_MODELS_ADVECTION_LAGRANGIAN_H
#define RODIN_MODELS_ADVECTION_LAGRANGIAN_H

#include <cstdlib>
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

  class DefaultBoundaryPolicy
  {
    public:
      constexpr
      bool operator()(Real&, Index&, Math::SpatialPoint&) const
      {
        return false;
      }
  };

  class DefaultTangentPolicy
  {
    public:
      constexpr
      bool operator()(Real&, Index&, Math::SpatialPoint&) const
      {
        return true;
      }
  };

  template <class Operand, class VectorField, class Step = Math::RungeKutta::RK4, class BoundaryPolicy = DefaultBoundaryPolicy, class TangentPolicy = DefaultTangentPolicy>
  class Flow;

  template <
    class Derived,
    class FES, Variational::ShapeFunctionSpaceType Space,
    class VectorField,
    class Step,
    class BoundaryPolicy,
    class TangentPolicy
    >
  class Flow<Variational::ShapeFunctionBase<Derived, FES, Space>, VectorField, Step, BoundaryPolicy, TangentPolicy>
  : public Variational::ShapeFunctionBase<
      Flow<Variational::ShapeFunctionBase<Derived, FES, Space>, VectorField, Step, BoundaryPolicy, TangentPolicy>, FES, Space>
  {
    public:
      using Operand = Variational::ShapeFunctionBase<Derived, FES, Space>;
      using VectorFieldType = VectorField;
      using StepType = Step;
      using Parent = Variational::ShapeFunctionBase<
        Flow<Operand, VectorFieldType, StepType, BoundaryPolicy, TangentPolicy>, FES, Space>;

      template <class S, class BP, class TP>
      Flow(Real t, const Operand& u, const VectorFieldType& velocity,
            S&& st = S(), BP&& bp = BP(), TP&& tp = TP())
        : m_t(t),
          m_operand(u.copy()),
          m_velocity(velocity),
          m_step(std::forward<S>(st)),
          m_bp(std::forward<BP>(bp)),
          m_tp(std::forward<TP>(tp)),
          m_p(std::nullopt)
      {}

      Flow(const Flow& other)
        : Parent(other),
          m_t(other.m_t),
          m_operand(other.m_operand->copy()),
          m_velocity(other.m_velocity),
          m_step(other.m_step),
          m_bp(other.m_bp),
          m_tp(other.m_tp),
          m_p(other.m_p)
      {}

      Flow(Flow&& other)
        : Parent(std::move(other)),
          m_t(std::exchange(other.m_t, 0)),
          m_operand(std::move(other.m_operand)),
          m_velocity(std::move(other.m_velocity)),
          m_step(std::move(other.m_step)),
          m_bp(std::move(other.m_bp)),
          m_tp(std::move(other.m_tp)),
          m_p(std::exchange(other.m_p, std::nullopt))
      {}

      std::optional<Geometry::Point> forward(const Geometry::Point& p) const
      {
        static thread_local Index s_cellIdx;
        static thread_local Math::SpatialPoint s_rc0{{}};
        static thread_local Math::SpatialPoint s_rc1{{}};
        static thread_local Math::SpatialPoint s_pc{{}};
        const auto& polytope0 = p.getPolytope();
        const auto& mesh = polytope0.getMesh();
        const size_t cellDim = mesh.getDimension();
        const auto& conn = mesh.getConnectivity();

        s_pc = p.getPhysicalCoordinates();

        if (polytope0.getDimension() == cellDim)
        {
          s_cellIdx = polytope0.getIndex();
          s_rc0 = p.getReferenceCoordinates();
        }
        else if (polytope0.getDimension() == cellDim - 1) // Start on a face
        {
          const Index fidx = polytope0.getIndex();
          const auto& adj = conn.getIncidence(cellDim - 1, cellDim).at(fidx);

          if (mesh.isBoundary(fidx)) // Start on a boundary face
          {
            const Index c = adj[0];
            mesh.getPolytopeTransformation(cellDim, c).inverse(s_rc1, s_pc);
            Geometry::Polytope::Project(mesh.getGeometry(cellDim, c)).cell(s_rc0, s_rc1);
            s_cellIdx = c;
          }
          else // Start on an interior face
          {
            assert(adj.size() == 2);
            const Index c0 = adj[0];
            const auto& faces0 = conn.getIncidence(cellDim, cellDim - 1).at(c0);
            size_t j0 = faces0.size();
            for (size_t k = 0; k < faces0.size(); ++k)
            {
              if (faces0[k] == fidx)
              {
                j0 = k;
                break;
              }
            }
            assert(j0 < faces0.size());
            mesh.getPolytopeTransformation(cellDim, c0).inverse(s_rc1, s_pc);
            const auto& cell0 = *mesh.getPolytope(cellDim, c0);
            const Geometry::Point q0(cell0, s_rc1, s_pc);
            const auto a0 = q0.getJacobianInverse() * m_velocity(q0);
            const auto& hs0 =
              Geometry::Polytope::Traits(mesh.getGeometry(cellDim, c0)).getHalfSpace();
            const auto nref0 = hs0.matrix.row(j0);
            if (nref0.dot(a0) > 0) // Flow into cell c1
            {
              const Index c1 = adj[1];
              mesh.getPolytopeTransformation(cellDim, c1).inverse(s_rc1, s_pc);
              Geometry::Polytope::Project(mesh.getGeometry(cellDim, c1)).cell(s_rc0, s_rc1);
              s_cellIdx = c1;
            }
            else // Flow into cell c0
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
          const Geometry::Point q(polytope, s_rc0, s_pc);
          const auto vr = q.getJacobianInverse() * m_velocity(q);
          if (vr.squaredNorm() == 0)
            break;
          const Geometry::Polytope::Type g = mesh.getGeometry(cellDim, s_cellIdx);
          const auto& faces = conn.getIncidence(cellDim, cellDim - 1).at(s_cellIdx);
          const auto& hs = Geometry::Polytope::Traits(g).getHalfSpace();

          Real exitTime = std::numeric_limits<Real>::infinity();
          Index local;
          bool tangent = true;
          for (size_t i = 0; i < faces.size(); i++)
          {
            decltype(auto) nref = hs.matrix.row(i);
            const Real denom = nref.dot(vr);
            if (denom > 0) // Outflow face
            {
              const Real numer = hs.vector[i] - nref.dot(s_rc0);
              const Real tf = numer / denom;
              if (tf < exitTime) // Find minimum exit time
              {
                local = i;
                exitTime = tf;
                tangent = false;
              }
            }
          }

          if (tangent) // Tangential flow
          {
            if(!m_tp(tau, s_cellIdx, s_rc0))
              return {};
          }
          else // Transversal flow
          {
            const Index face = faces[local];

            if (tau < exitTime) // Does not exit cell
            {
              m_step.step(
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
            else // Exits cell
            {
              s_rc0 += exitTime * vr;
              if (mesh.isBoundary(face)) // Exits domain
              {
                if(!m_bp(tau, s_cellIdx, s_rc0))
                  return {};
              }
              else // Crosses face into adjacent cell
              {
                mesh.getPolytopeTransformation(cellDim, s_cellIdx).transform(s_pc, s_rc0);
                const auto& cells = conn.getIncidence(cellDim - 1, cellDim).at(face);
                assert(cells.size() == 2);
                s_cellIdx = (cells[0] == s_cellIdx) ? cells[1] : cells[0];
                mesh.getPolytopeTransformation(cellDim, s_cellIdx).inverse(s_rc0, s_pc);
                Geometry::Polytope::Project(mesh.getGeometry(cellDim, s_cellIdx)).cell(s_rc1, s_rc0);
                s_rc0 = s_rc1;
                tau -= exitTime;
              }
            }
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
      Step m_step;
      BoundaryPolicy m_bp;
      TangentPolicy m_tp;
      std::optional<Geometry::Point> m_p;
  };

  template <class FES, class Data, class VectorField, class Step>
  class Lagrangian<Variational::GridFunction<FES, Data>, VectorField, Step>
  {
    public:
      using FESType = FES;
      using DataType = Data;
      using VectorFieldType = VectorField;
      using SolutionType = Variational::GridFunction<FES, Data>;
      using StepType = Step;
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
      Flow<SolutionType, VectorFieldType, StepType> m_pullback;
      Step m_rk;
  };
}

#endif
