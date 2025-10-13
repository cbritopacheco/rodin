/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FLOW_H
#define RODIN_VARIATIONAL_FLOW_H

#include "ForwardDecls.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Region.h"

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/RungeKutta/RK4.h"
#include "Rodin/Math/RootFinding/NewtonRaphson.h"

#include "Rodin/QF/GenericPolytopeQuadrature.h"
#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Variational
{
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

  template <
    class Derived,
    class FES,
    ShapeFunctionSpaceType Space,
    class VectorField,
    class Step,
    class Root,
    class BoundaryPolicy,
    class TangentPolicy
  > class Flow<ShapeFunctionBase<Derived, FES, Space>, VectorField, Step, Root, BoundaryPolicy, TangentPolicy>
  : public ShapeFunctionBase<
      Flow<ShapeFunctionBase<Derived, FES, Space>, VectorField, Step, Root, BoundaryPolicy, TangentPolicy>, FES, Space>
  {
    public:
      using Operand =
        ShapeFunctionBase<Derived, FES, Space>;

      using VectorFieldType =
        VectorField;

      using FESType =
        FES;

      using StepType =
        Step;

      using RootType =
        Root;

      using BoundaryPolicyType =
        BoundaryPolicy;

      using TangentPolicyType =
        TangentPolicy;

      using ScalarType =
        typename FormLanguage::Traits<FESType>::ScalarType;

      using Parent =
        ShapeFunctionBase<
          Flow<Operand, VectorFieldType, StepType, RootType, BoundaryPolicy, TangentPolicy>, FES, Space>;

      template <
        class VVel,
        class S = StepType, class R = RootType,
        class B = BoundaryPolicy, class T = TangentPolicy>
      Flow(const Real& t,
           const Operand& u,
           VVel&& vel,
           S&& st=S{}, R&& rt=R{}, B&& bp=B{}, T&& tp=T{})
        : Parent(u.getFiniteElementSpace()),
          m_t(t),
          m_operand(u.copy()),
          m_velocity(std::forward<VVel>(vel)),
          m_step(std::forward<S>(st)),
          m_root(std::forward<R>(rt)),
          m_bp(std::forward<B>(bp)),
          m_tp(std::forward<T>(tp)),
          m_p(nullptr)
      {}

      Flow(const Flow& other)
        : Parent(other),
          m_t(other.m_t),
          m_operand(other.m_operand->copy()),
          m_velocity(other.m_velocity),
          m_step(other.m_step),
          m_root(other.m_root),
          m_bp(other.m_bp),
          m_tp(other.m_tp),
          m_p(other.m_p),
          m_trace(other.m_trace)
      {}

      Flow(Flow&& other)
        : Parent(std::move(other)),
          m_t(std::exchange(other.m_t, 0)),
          m_operand(std::move(other.m_operand)),
          m_velocity(std::move(other.m_velocity)),
          m_step(std::move(other.m_step)),
          m_root(std::move(other.m_root)),
          m_bp(std::move(other.m_bp)),
          m_tp(std::move(other.m_tp)),
          m_p(std::exchange(other.m_p, nullptr)),
          m_trace(std::move(other.m_trace))
      {}

      std::optional<Geometry::Point> forward(const Geometry::Point& p) const
      {
        static thread_local Index s_cellIdx;
        static thread_local Math::SpatialPoint s_rc0{{}};
        static thread_local Math::SpatialPoint s_rc1{{}};
        static thread_local Math::SpatialPoint s_rc_tau{{}};
        static thread_local Math::SpatialPoint s_pc{{}};

        const auto& polytope0 = p.getPolytope();
        const auto& mesh = polytope0.getMesh();
        const size_t cellDim = mesh.getDimension();
        const auto& conn = mesh.getConnectivity();

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
            assert(adj.size() == 1);
            const Index c = adj[0];
            mesh.getPolytopeTransformation(cellDim, c).inverse(s_rc1, p.getPhysicalCoordinates());
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
            const auto& pc = p.getPhysicalCoordinates();
            mesh.getPolytopeTransformation(cellDim, c0).inverse(s_rc1, pc);
            const auto& cell0 = *mesh.getPolytope(cellDim, c0);
            const Geometry::Point q0(cell0, s_rc1, pc);
            const auto a0 = q0.getJacobianInverse() * m_velocity(q0);
            const auto& hs0 =
              Geometry::Polytope::Traits(mesh.getGeometry(cellDim, c0)).getHalfSpace();
            const auto nref0 = hs0.matrix.row(j0);
            if (nref0.dot(a0) > 0) // Flow into cell c1
            {
              const Index c1 = adj[1];
              mesh.getPolytopeTransformation(cellDim, c1).inverse(s_rc1, pc);
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
        else
        {
          assert(false);
          return {};
        }

        s_rc1.resizeLike(s_rc0);
        s_rc_tau.resizeLike(s_rc0);

        Real tau = std::abs(m_t);
        const Real sgn = Math::sgn(m_t);
        while (tau > 0)
        {
          const auto it = mesh.getPolytope(cellDim, s_cellIdx);
          const auto& cell = *it;
          const Geometry::Polytope::Type g = mesh.getGeometry(cellDim, s_cellIdx);
          const auto& faces = conn.getIncidence(cellDim, cellDim - 1).at(s_cellIdx);
          const auto& hs = Geometry::Polytope::Traits(g).getHalfSpace();

          // Field in reference space
          const auto vr = [&](const Math::SpatialPoint& rc)
          {
            static thread_local Math::SpatialVector<Real> s_vr;
            const Geometry::Point qp(cell, rc);
            s_vr = sgn * qp.getJacobianInverse() * m_velocity(qp);
            return s_vr;
          };

          std::optional<Real> tmin;
          Index local = faces.size();

          if constexpr (requires { m_step.dense(s_rc0, tau, vr); })
          {
            // Dense-output path (fastest)
            const auto D = m_step.dense(s_rc0, tau, vr);
            s_rc_tau = D.X(tau);

            for (size_t i = 0; i < faces.size(); ++i)
            {
              const auto nref = hs.matrix.row(i);
              const Real bf = hs.vector[i];

              const Real g0 = bf - nref.dot(s_rc0);
              const Real gtau = bf - nref.dot(s_rc_tau);
              if (g0 * gtau < 0) // root in (0, tau)
              {
                const Real t0 = Real(0.5) * tau;
                auto event = [&](Real& t)
                {
                  const auto X = D.X(t);
                  const auto V = D.V(t);
                  return std::pair{ bf - nref.dot(X), -nref.dot(V) };
                };
                if (auto rt = m_root.solve(event, t0, Real(0), tau))
                {
                  const Real t = *rt;
                  if (!tmin.has_value() || (t < *tmin)) { tmin = t; local = i; }
                }
              }
            }
          }
          else
          {
            s_rc_tau = s_rc0;

            m_step.step(s_rc_tau, tau, s_rc0, vr);

            for (size_t i = 0; i < faces.size(); ++i)
            {
              const auto nref = hs.matrix.row(i);
              const Real bf = hs.vector[i];

              const Real f0 = bf - nref.dot(s_rc0);
              const Real ftau = bf - nref.dot(s_rc_tau);

              if (f0 * ftau < 0) // root in (0, tau)
              {
                const Real t0 = Real(0.5) * tau;
                const auto event = [&](Real& t)
                {
                  m_step.step(s_rc1, t, s_rc0, vr);
                  decltype(auto) vrv = vr(s_rc1);
                  return std::pair{ bf - nref.dot(s_rc1), -nref.dot(vrv) };
                };

                if (const auto rt = m_root.solve(event, t0, Real(0), tau))
                {
                  const Real t = *rt;
                  if (!tmin.has_value() || (t < *tmin))
                  {
                    tmin = t;
                    local = i;
                  }
                }
              }
            }
          }

          // No face will be hit -> advance whole remaining time
          if (!tmin.has_value())
          {
            if (!m_tp(tau, s_cellIdx, s_rc0))
              return {};
            m_step.step(s_rc1, tau, s_rc0, vr);
            s_rc0 = s_rc1;
            break;
          }

          assert(local < faces.size());

          // Hit the closest face at t*
          const Real tstar = *tmin;
          const Index face = faces[local];

          // Arrive on face (in ref of current cell)
          m_step.step(s_rc1, tstar, s_rc0, vr);
          s_rc0 = s_rc1;
          tau -= tstar;

          if (mesh.isBoundary(face))
          {
            // Let boundary policy decide (may update tau/s_cellIdx/s_rc0)
            if (!m_bp(tau, s_cellIdx, s_rc0))
              return {};
            // The policy is responsible for keeping tau > 0 if it wants to continue.
            continue;
          }

          // Cell hop
          mesh.getPolytopeTransformation(cellDim, s_cellIdx).transform(s_pc, s_rc0);
          const auto& cells = conn.getIncidence(cellDim - 1, cellDim).at(face);
          assert(cells.size() == 2);
          s_cellIdx = (cells[0] == s_cellIdx) ? cells[1] : cells[0];
          mesh.getPolytopeTransformation(cellDim, s_cellIdx).inverse(s_rc0, s_pc);
          Geometry::Polytope::Project(mesh.getGeometry(cellDim, s_cellIdx)).cell(s_rc1, s_rc0);
          s_rc0 = s_rc1;
        }
        return Geometry::Point(*mesh.getPolytope(cellDim, s_cellIdx), s_rc0);
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        return m_operand->getDOFs(polytope);
      }

      constexpr
      ScalarType getBasis(size_t local) const
      {
        if (m_trace)
          return m_operand->getBasis(local);
        else
          return ScalarType(0);
      }

      const Geometry::Point& getPoint() const
      {
        assert(m_p);
        return *m_p;
      }

      Flow& setPoint(const Geometry::Point& p)
      {
        m_p = &p;
        if (auto tr = this->forward(p))
        {
          m_trace.emplace(std::move(*tr));   // construct in-place
          m_operand->setPoint(*m_trace);
        }
        else
        {
          m_trace.reset();
        }
        return *this;
      }

      const auto& getTrace() const
      {
        return m_trace;
      }

      constexpr
      const auto& getLeaf() const
      {
        return m_operand->getLeaf();
      }

      const auto& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      virtual Flow* copy() const noexcept override
      {
        return new Flow(*this);
      }

    private:
      const Real m_t;
      std::unique_ptr<Operand> m_operand;
      VectorFieldType m_velocity;
      Step m_step;
      Root m_root;
      BoundaryPolicy m_bp;
      TangentPolicy m_tp;
      const Geometry::Point* m_p;
      std::optional<Geometry::Point> m_trace;
  };

  template <class D, class FES, ShapeFunctionSpaceType S, class VVel>
  Flow(const Real&, const ShapeFunctionBase<D, FES, S>&, VVel&&)
    -> Flow<
         ShapeFunctionBase<D, FES, S>,
         VVel,                       // keep T or T&
         Math::RungeKutta::RK4,      // value default
         Math::RootFinding::NewtonRaphson<typename FormLanguage::Traits<FES>::ScalarType>,
         DefaultBoundaryPolicy,
         DefaultTangentPolicy>;

  template <class D, class FES, ShapeFunctionSpaceType S, class VVel, class SStep>
  Flow(const Real&, const ShapeFunctionBase<D, FES, S>&, VVel&&, SStep&&)
    -> Flow<
         ShapeFunctionBase<D,FES,S>,
         VVel,
         SStep,
         Math::RootFinding::NewtonRaphson<typename FormLanguage::Traits<FES>::ScalarType>,
         DefaultBoundaryPolicy,
         DefaultTangentPolicy>;

  template <
    class D,
    class FES,
    ShapeFunctionSpaceType S,
    class VVel, class SStep, class RRoot, class BBP, class TTP>
  Flow(const Real&, const ShapeFunctionBase<D,FES,S>&, VVel&&, SStep&&, RRoot&&, BBP&&, TTP&&)
    -> Flow<ShapeFunctionBase<D, FES, S>, VVel, SStep, RRoot, BBP, TTP>;

  template <
    class Derived,
    class FES,
    ShapeFunctionSpaceType Space,
    class VectorField,
    class Step,
    class Root,
    class BoundaryPolicy,
    class TangentPolicy
  > class QuadratureRule<Flow<ShapeFunctionBase<Derived, FES, Space>, VectorField, Step, Root, BoundaryPolicy, TangentPolicy>>
    : public LinearFormIntegratorBase<
        typename FormLanguage::Traits<
          Flow<ShapeFunctionBase<Derived, FES, Space>, VectorField, Step, Root, BoundaryPolicy, TangentPolicy>>::ScalarType>
  {
    public:
      using FESType =
        FES;

      using IntegrandType =
        Flow<ShapeFunctionBase<Derived, FES, Space>, VectorField, Step, Root, BoundaryPolicy, TangentPolicy>;

      using ScalarType =
        typename FormLanguage::Traits<IntegrandType>::ScalarType;

      using Parent = LinearFormIntegratorBase<ScalarType>;

      QuadratureRule(const IntegrandType& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy()),
          m_polytope(nullptr)
      {}

      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_polytope(other.m_polytope),
          m_qf(other.m_qf),
          m_ps(other.m_ps)
      {}

      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_polytope(std::exchange(other.m_polytope, nullptr)),
          m_qf(std::move(other.m_qf)),
          m_ps(std::move(other.m_ps))
      {}

      ScalarType integrate(size_t local)
      {
        return ScalarType(0);
      }

      QuadratureRule& setPolytope(const Geometry::Polytope& polytope)
      {
        m_polytope = &polytope;
        auto& integrand = *m_integrand;
        auto& scatter = this->getScatter();
        const auto& fes = integrand.getFiniteElementSpace();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& fe = fes.getFiniteElement(d, idx);
        const auto& u = integrand.getOperand();
        m_qf.reset(new QF::GenericPolytopeQuadrature(polytope.getGeometry()));
        m_ps.clear();
        m_ps.reserve(m_qf->getSize());
        for (size_t i = 0; i < m_qf->getSize(); i++)
          m_ps.emplace_back(polytope, m_qf->getPoint(i));
        scatter.clear();
        for (size_t i = 0; i < m_ps.size(); ++i)
        {
          const auto& p = m_ps[i];
          const Real w = m_qf->getWeight(i) * p.getDistortion();
          m_integrand.setPoint(p);
          const auto& trace = m_integrand.getTrace();
          if (!trace)
            continue;
          const auto& cell = trace->getPolytope();
          const size_t cellDim = cell.getDimension();
          const Index cellIdx = cell.getIndex();
          decltype(auto) fe = fes.getFiniteElement(cellDim, cellIdx);
          decltype(auto) up = u.getValue(p);
          for (size_t local = 0; local < fe.getCount(); ++local)
          {
            const Index global = fes.getGlobalIndex({ cellDim, cellIdx }, local);
            decltype(auto) basis = fe.getBasis(local);
            decltype(auto) mapping = fes.getInverseMapping({ cellDim, cellIdx }, basis);
            scatter.push(global, w * up * mapping(*trace));
          }
        }
        return *this;
      }

      constexpr
      const IntegrandType& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      const Geometry::Polytope& getPolytope() const override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      virtual Geometry::Region getRegion() const = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<IntegrandType> m_integrand;

      std::unique_ptr<QF::QuadratureFormulaBase> m_qf;
      std::vector<Geometry::Point> m_ps;

      const Geometry::Polytope* m_polytope;
  };
}

#endif
