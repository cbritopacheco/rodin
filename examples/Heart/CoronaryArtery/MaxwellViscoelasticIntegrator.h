#ifndef RODIN_EXAMPLES_HEART_MAXWELLVISCOELASTICINTEGRATOR_H
#define RODIN_EXAMPLES_HEART_MAXWELLVISCOELASTICINTEGRATOR_H

/**
 * @file MaxwellViscoelasticIntegrator.h
 * @brief History (memory) term of a single Maxwell viscoelastic branch.
 *
 * One Maxwell branch evolves an internal PK1-stress tensor P_v by
 *
 *   dP_v/dt + (1/tau) P_v = (g/tau) P_eq ,
 *
 * with P_eq the equilibrium first Piola-Kirchhoff stress.  A BDF1 step gives
 * the recursion
 *
 *   P_v^{n+1} = ( tau P_v^n + g dt P_eq^{n+1} ) / (tau + dt)
 *             = decay * P_v^n + frac * P_eq^{n+1},
 *   decay = tau/(tau+dt),   frac = g dt/(tau+dt).
 *
 * If the total stress is S = P_eq + P_v, the branch contribution to the solid
 * residual int S:grad(w) splits as
 *
 *   int P_v^{n+1} : grad(w)
 *     = frac * int P_eq^{n+1}:grad(w)        (= frac * InternalForce, CURRENT)
 *     + decay * int P_v^n   :grad(w)         (HISTORY, this integrator).
 *
 * So the driver scales the existing solid internal force / tangent by
 * (1 + frac) and adds this history linear form (a KNOWN force; no tangent,
 * since P_v^n is data).  P_v^n is stored as a P0 (element-constant) tensor
 * field, FLATTENED row-major into a vdim = d*d vector field:
 *
 *   Pv_value[c*d + j] = (P_v)_{c j},   c = spatial row, j = reference column.
 *
 * The same flattening must be used when the driver updates P_v^n -> P_v^{n+1}.
 *
 * The contraction with the test gradient, for a vector test basis Psi_b with
 * reference Jacobian Jref and inverse geometry Jacobian Jinv, is
 *
 *   P_v : grad Psi_b = sum_{c,j} (P_v)_{c j} (grad Psi_b)_{c j},
 *   (grad Psi_b)_{c j} = sum_r Jref(c,r) Jinv(r,j).
 *
 * Implementation mirrors VMSGradDivLinearIntegrator (same base class, same
 * quadrature / Jacobian handling); only the local contraction differs.
 */

#include <cassert>
#include <cstdint>
#include <stdexcept>

#include <Rodin/Math/SpatialVector.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Variational/IntegrationPoint.h>

namespace Rodin::Examples::Heart
{
  /**
   * @brief Linear integrator assembling  int_K  scale * (P_v : grad v).
   *
   * @tparam TestFunction  Vector (displacement) test field type.
   * @tparam StoredTensor  Stored P_v field (vdim = d*d, row-major flattening);
   *                       must expose getValue(IntegrationPoint) -> vector.
   */
  template <class TestFunction, class StoredTensor>
  class MaxwellHistoryLinearIntegrator final
    : public Variational::LinearFormIntegratorBase<
        typename TestFunction::ScalarType>
  {
    public:
      using ScalarType = typename TestFunction::ScalarType;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;

      /**
       * @param[in] v     Vector test field.
       * @param[in] pv    Stored internal stress P_v^n (flattened, vdim = d*d).
       * @param[in] scale Scalar premultiplier (e.g. decay = tau/(tau+dt)).
       */
      MaxwellHistoryLinearIntegrator(
          const TestFunction& v, const StoredTensor& pv, ScalarType scale)
        : Parent(v.getLeaf()),
          m_v(v),
          m_pv(pv),
          m_polytope(nullptr),
          m_scale(scale)
      {}

      MaxwellHistoryLinearIntegrator(const MaxwellHistoryLinearIntegrator&) = default;
      MaxwellHistoryLinearIntegrator(MaxwellHistoryLinearIntegrator&&) = default;

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      MaxwellHistoryLinearIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const auto geometry = polytope.getGeometry();
        const auto idx = polytope.getIndex();

        const auto& fes = m_v.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);

        const size_t nte = m_v.getDOFs(polytope);
        const size_t vdim = fes.getVectorDimension();

        if (vdim != d)
          throw std::runtime_error(
            "MaxwellHistoryLinearIntegrator expects a vector test space whose "
            "dimension matches the cell dimension.");
        if (nte != fe.getCount())
          throw std::runtime_error(
            "MaxwellHistoryLinearIntegrator expects element-local DOFs to match "
            "the finite-element basis count.");

        const size_t pvOrder = getFunctionOrder(m_pv, polytope, size_t(0));
        const size_t qOrder = pvOrder + derivativeOrder(fe.getOrder());

        const auto& qf = QF::PolytopeQuadratureFormula::get(qOrder, geometry);
        const auto& q = polytope.getQuadrature(qf);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const Variational::IntegrationPoint ip(p, &qf, qp);

          const ScalarType wdet =
            static_cast<ScalarType>(qf.getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();
          const auto& rc = p.getReferenceCoordinates();

          // P_v at this point, flattened row-major: pv[c*d + j] = (P_v)_{c j}.
          const auto pv = m_pv.getValue(ip);
          if (pv.size() != d * d)
            throw std::runtime_error(
              "MaxwellHistoryLinearIntegrator expects the stored tensor field to "
              "have vdim = d*d.");

          for (size_t b = 0; b < nte; ++b)
            m_vec(static_cast<Eigen::Index>(b)) +=
              wdet * m_scale * contract(fe.getBasis(b), rc, Jinv, pv, d);
        }

        return *this;
      }

      ScalarType integrate(size_t local) final override
      {
        return m_vec(static_cast<Eigen::Index>(local));
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      MaxwellHistoryLinearIntegrator* copy() const noexcept final override
      {
        return new MaxwellHistoryLinearIntegrator(*this);
      }

    private:
      // P_v : grad Psi = sum_{c,j} pv[c*d+j] * sum_r Jref(c,r) Jinv(r,j).
      template <class Basis, class JInv, class Vector>
      static ScalarType contract(
          const Basis& basis,
          const Math::SpatialVector<Real>& rc,
          const JInv& Jinv,
          const Vector& pv,
          size_t d)
      {
        const auto Jref = basis.getJacobian()(rc);
        ScalarType s = 0;
        for (size_t c = 0; c < d; ++c)
        {
          for (size_t j = 0; j < d; ++j)
          {
            ScalarType gradPhys = 0;
            for (size_t r = 0; r < d; ++r)
              gradPhys +=
                Jref(static_cast<std::uint8_t>(c), static_cast<std::uint8_t>(r))
                * Jinv(static_cast<std::uint8_t>(r), static_cast<std::uint8_t>(j));
            s += pv(static_cast<std::uint8_t>(c * d + j)) * gradPhys;
          }
        }
        return s;
      }

      static size_t derivativeOrder(size_t order) noexcept
      {
        return order == 0 ? 0 : order - 1;
      }

      template <class Function>
      static size_t getFunctionOrder(
          const Function& f, const Geometry::Polytope& polytope, size_t fallback)
      {
        if constexpr (requires { f.getOrder(polytope); })
        {
          const auto order = f.getOrder(polytope);
          if (order.has_value())
            return *order;
        }
        return fallback;
      }

      const TestFunction& m_v;
      const StoredTensor& m_pv;
      const Geometry::Polytope* m_polytope;
      ScalarType m_scale;

      Math::Vector<ScalarType> m_vec;
  };
}

#endif
