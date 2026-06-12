#ifndef RODIN_EXAMPLES_HEART_VMSCONVECTIONINTEGRATOR_H
#define RODIN_EXAMPLES_HEART_VMSCONVECTIONINTEGRATOR_H

/**
 * @file VMSConvectionIntegrator.h
 * @brief Dynamic projected-VMS convective stabilization integrators.
 *
 * This file defines two custom cell integrators used by the coronary
 * Navier--Stokes example.
 *
 * The stabilized weak contribution implemented here is:
 *
 * @f[
 *   \sum_K
 *   \int_K
 *     \tau_K \rho^2
 *     \bigl((\nabla u_h^{n+1})u_h^n\bigr)
 *     \cdot
 *     \bigl((\nabla v_h)u_h^n\bigr)
 *   \, dx
 *   -
 *   \sum_K
 *   \int_K
 *     \rho
 *     \left(
 *       \tau_K \rho \Pi_h[(\nabla u_h^n)u_h^n]
 *       +
 *       u_h^{\prime,n+1}
 *     \right)
 *     \cdot
 *     \bigl((\nabla v_h)u_h^n\bigr)
 *   \, dx .
 * @f]
 *
 * The dynamic subscale @f$u_h^{\prime,n+1}@f$ is not computed inside this
 * header. It is computed in the driver and passed as @c OldSubScale.  The
 * subscale is the (dissipative) response to the unresolved convective
 * residual, @f$u' = -\tau_K R@f$, so the projection-minus-residual term enters
 * with @c Pi MINUS the raw convective acceleration (NOT the other way round):
 *
 * @f[
 *   u_h^{\prime,n+1}
 *   =
 *   \tau_K \rho
 *   \left(
 *     \Pi_h[(\nabla u_h^n)u_h^n]
 *     -
 *     (\nabla u_h^n)u_h^n
 *   \right)
 *   +
 *   \tau_K \frac{\rho}{\Delta t} u_h^{\prime,n}.
 * @f]
 *
 * Equivalently, as written in the driver:
 * @f$ u' = \tau_K \rho\,(\tfrac{1}{\Delta t} u_h^{\prime,n}
 *          - ((\nabla u_h^n)u_h^n - \Pi_h[(\nabla u_h^n)u_h^n])) @f$.
 *
 * The stabilization parameter is:
 *
 * @f[
 *   \tau_K
 *   =
 *   \frac{\texttt{vmsScale}}
 *   {
 *     \sqrt{
 *       \left(2\rho/\Delta t\right)^2
 *       +
 *       \left(c_2 \rho |u_h^n|/h_K\right)^2
 *       +
 *       \left(c_1 \mu_h/h_K^2\right)^2
 *     }
 *   }.
 * @f]
 *
 * This implementation is intentionally specialized:
 *
 * - cell integrals in the mesh spatial dimension,
 * - vector-valued velocity fields whose basis functions provide Jacobians,
 * - explicit frozen convective velocity @f$u_h^n@f$,
 * - projected stabilization parameter @f$\tau_K@f$,
 * - local stabilization length @f$h_K = |K|^{1/d}@f$ for cell dimension
 *   @f$d@f$, computed outside this header as part of @f$\tau_K@f$.
 *
 * The integrators below do not implement a full residual-based VMS
 * Navier--Stokes method. They implement only the convective, projected,
 * term-by-term stabilization used by the coronary example experiments. The
 * advection velocity, projected convective acceleration, projected tau, and
 * dynamic subscale are treated as already known fields. Consequently this
 * contribution is suitable as a lagged Oseen/Picard stabilization term. If it
 * is used inside a Newton solve, it should be understood as a frozen
 * stabilization contribution unless the caller also differentiates these
 * external projected fields.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <Rodin/Math/SpatialVector.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Variational/H1/H1Element.h>
#include <Rodin/Variational/IntegrationPoint.h>

namespace Rodin::Examples::Heart
{
  /**
   * @brief Cell integrator for the bilinear projected-VMS convective term.
   *
   * This class adds streamline-oriented dissipation in the direction of the
   * frozen velocity @f$u_h^n@f$. For @f$\tau_K \ge 0@f$, the square case
   * satisfies:
   *
   * @f[
   *   A_K(w_h,w_h)
   *   =
   *   \int_K
   *     \tau_K \rho^2
   *     \left|(\nabla w_h)u_h^n\right|^2
   *   \, dx
   *   \ge 0 .
   * @f]
   *
   * Thus the term damps unresolved convective derivatives without adding
   * isotropic viscosity. It is not a pressure, divergence, or viscous residual
   * stabilization.
   *
   * This assembles:
   *
   * @f[
   *   A_K(u_h, v_h)
   *   =
   *   \int_K
   *     \tau_K \rho^2
   *     \bigl((\nabla u_h)u_h^n\bigr)
   *     \cdot
   *     \bigl((\nabla v_h)u_h^n\bigr)
   *   \, dx .
   * @f]
   *
   * For vector basis functions @f$\Phi_a@f$ and @f$\Psi_b@f$, the local entry is:
   *
   * @f[
   *   A_{ba}
   *   =
   *   \int_K
   *     \tau_K \rho^2
   *     \left((\nabla \Phi_a)u_h^n\right)
   *     \cdot
   *     \left((\nabla \Psi_b)u_h^n\right)
   *   \, dx .
   * @f]
   *
   * The implementation is finite-element-order independent for vector finite
   * elements whose basis functions expose a Jacobian. Local entries are formed
   * directly from @f$(\nabla \Phi_a)u_h^n@f$ rather than by assuming a scalar
   * P2 block structure.
   *
   * @tparam TrialFunction Trial velocity field type.
   * @tparam TestFunction Test velocity field type.
   * @tparam OldVelocity Frozen convective velocity field type.
   * @tparam ProjectedTau Projected stabilization parameter field type.
   */
  template <class TrialFunction, class TestFunction, class OldVelocity, class ProjectedTau>
  class VMSConvectionBilinearIntegrator final
    : public Variational::LocalBilinearFormIntegratorBase<
        typename TrialFunction::ScalarType>
  {
    public:
      using ScalarType = typename TrialFunction::ScalarType;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructor.
       *
       * @param[in] u Trial velocity @f$u_h@f$.
       * @param[in] v Test velocity @f$v_h@f$.
       * @param[in] uOld Frozen convective velocity @f$u_h^n@f$.
       * @param[in] rho Fluid density @f$\rho@f$.
       * @param[in] tau Projected stabilization parameter @f$\tau_K@f$.
       */
      VMSConvectionBilinearIntegrator(
          const TrialFunction& u,
          const TestFunction& v,
          const OldVelocity& uOld,
          const ProjectedTau& tau,
          ScalarType rho)
        : Parent(u.getLeaf(), v.getLeaf()),
          m_u(u),
          m_v(v),
          m_uOld(uOld),
          m_tau(tau),
          m_polytope(nullptr),
          m_rho(rho)
      {}

      VMSConvectionBilinearIntegrator(
          const VMSConvectionBilinearIntegrator&) = default;

      VMSConvectionBilinearIntegrator(
          VMSConvectionBilinearIntegrator&&) = default;

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Builds the local matrix on one cell.
       */
      VMSConvectionBilinearIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const auto geometry = polytope.getGeometry();
        const auto idx = polytope.getIndex();

        const auto& trialFES = m_u.getFiniteElementSpace();
        const auto& testFES  = m_v.getFiniteElementSpace();

        const auto& trialFE = trialFES.getFiniteElement(d, idx);
        const auto& testFE  = testFES .getFiniteElement(d, idx);

        const size_t ntr = m_u.getDOFs(polytope);
        const size_t nte = m_v.getDOFs(polytope);
        const size_t vdim = trialFES.getVectorDimension();

        if (vdim == 0 || vdim != testFES.getVectorDimension())
          throw std::runtime_error("VMSConvectionBilinearIntegrator expects matching vector-valued spaces.");

        if (ntr != trialFE.getCount() || nte != testFE.getCount())
          throw std::runtime_error("VMSConvectionBilinearIntegrator expects element-local DOFs to match finite-element basis count.");

        /*
         * Conservative quadrature choice for
         *
         *   tau (grad Phi uOld) . (grad Psi uOld).
         *
         * The tau field is often projected from a non-polynomial expression,
         * so exactness is not guaranteed. This still integrates the polynomial
         * part of the lagged convective stabilization without the old P2
         * underintegration.
         */
        const size_t uOldOrder =
          getFunctionOrder(m_uOld, polytope, trialFE.getOrder());
        const size_t tauOrder =
          getFunctionOrder(m_tau, polytope, size_t(0));
        const size_t qOrder =
          tauOrder
          + derivativeOrder(trialFE.getOrder())
          + derivativeOrder(testFE.getOrder())
          + 2 * uOldOrder;

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qOrder, geometry);

        const auto& q = polytope.getQuadrature(qf);

        m_mat.resize(
          static_cast<Eigen::Index>(nte),
          static_cast<Eigen::Index>(ntr));

        m_mat.setZero();

        ScalarType* A = m_mat.data();

        std::vector<Math::SpatialVector<ScalarType>> trDir(ntr);
        std::vector<Math::SpatialVector<ScalarType>> teDir(nte);

        for (auto& t : trDir)
          t.resize(static_cast<std::uint8_t>(vdim));

        for (auto& t : teDir)
          t.resize(static_cast<std::uint8_t>(vdim));

        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const Variational::IntegrationPoint ip(p, &qf, qp);

          const ScalarType wdet =
            static_cast<ScalarType>(qf.getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();
          const auto& rc = p.getReferenceCoordinates();

          const auto uOld = m_uOld.getValue(ip);
          const auto tau = m_tau.getValue(ip);

          if (uOld.size() != d)
            throw std::runtime_error("VMSConvectionBilinearIntegrator expects the frozen velocity dimension to match the cell dimension.");

          /*
           * Directional derivatives along frozen velocity:
           *
           *   trDir[a] = (grad Phi_a) uOld,
           *   teDir[b] = (grad Psi_b) uOld.
           */
          for (size_t a = 0; a < ntr; ++a)
            fillDirectionalDerivative(trDir[a], trialFE.getBasis(a), rc, Jinv, uOld);

          for (size_t b = 0; b < nte; ++b)
            fillDirectionalDerivative(teDir[b], testFE.getBasis(b), rc, Jinv, uOld);

          for (size_t b = 0; b < nte; ++b)
          {
            for (size_t a = 0; a < ntr; ++a)
            {
              /*
               * Local vector-basis entry:
               *
               *   int_K tau rho^2
               *     ((grad Phi_a) uOld)
               *     ·
               *     ((grad Psi_b) uOld).
               */
              const ScalarType kij =
                wdet * tau * m_rho * m_rho * Math::dot(trDir[a], teDir[b]);

              A[b * ntr + a] += kij;
            }
          }
        }

        return *this;
      }

      ScalarType integrate(size_t trial, size_t test) final override
      {
        return m_mat(
          static_cast<Eigen::Index>(test),
          static_cast<Eigen::Index>(trial));
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      VMSConvectionBilinearIntegrator* copy() const noexcept final override
      {
        return new VMSConvectionBilinearIntegrator(*this);
      }

    private:
      /**
       * @brief Computes @f$(\nabla \Phi)uOld@f$ for one vector basis.
       *
       * The vector finite-element basis supplies the reference Jacobian
       * @f$\hat\nabla\Phi@f$. Physical derivatives are obtained by
       * @f$\nabla_x\Phi = \hat\nabla\Phi J^{-1}@f$.
       */
      template <class Basis, class JInv, class Vector>
      static void fillDirectionalDerivative(
          Math::SpatialVector<ScalarType>& out,
          const Basis& basis,
          const Math::SpatialVector<Real>& rc,
          const JInv& Jinv,
          const Vector& uOld)
      {
        const auto Jref = basis.getJacobian()(rc);
        const size_t vdim = out.size();
        const size_t d = uOld.size();

        out.setZero();
        for (size_t c = 0; c < vdim; ++c)
        {
          for (size_t j = 0; j < d; ++j)
          {
            ScalarType gradPhys = 0;
            for (size_t r = 0; r < d; ++r)
            {
              gradPhys +=
                Jref(static_cast<std::uint8_t>(c), static_cast<std::uint8_t>(r))
                * Jinv(static_cast<std::uint8_t>(r), static_cast<std::uint8_t>(j));
            }
            out(static_cast<std::uint8_t>(c)) +=
              gradPhys * uOld(static_cast<std::uint8_t>(j));
          }
        }
      }

      static size_t derivativeOrder(size_t order) noexcept
      {
        return order == 0 ? 0 : order - 1;
      }

      template <class Function>
      static size_t getFunctionOrder(
          const Function& f,
          const Geometry::Polytope& polytope,
          size_t fallback)
      {
        if constexpr (requires { f.getOrder(polytope); })
        {
          const auto order = f.getOrder(polytope);
          if (order.has_value())
            return *order;
        }
        return fallback;
      }

      const TrialFunction& m_u;
      const TestFunction& m_v;
      const OldVelocity& m_uOld;
      const ProjectedTau& m_tau;

      const Geometry::Polytope* m_polytope;
      ScalarType m_rho;

      Eigen::Matrix<
        ScalarType,
        Eigen::Dynamic,
        Eigen::Dynamic,
        Eigen::RowMajor> m_mat;
  };

  /**
   * @brief Cell integrator for the grad-div (continuity) stabilization.
   *
   * Assembles the symmetric positive-semidefinite term
   *
   * @f[
   *   C_K(u_h, v_h)
   *   =
   *   \int_K \tau_C\,(\nabla\cdot u_h)(\nabla\cdot v_h)\, dx ,
   * @f]
   *
   * with a projected, element-local parameter @f$\tau_C@f$ (units of dynamic
   * viscosity, e.g. @f$\tau_C = \gamma\,\rho\,|u_h^n|\,h_K@f$).
   *
   * This is the consistent VMS continuity / "div-u" stabilization for an
   * inf-sup STABLE velocity-pressure pair (Taylor-Hood @f$P_2/P_1@f$): it
   * vanishes identically for divergence-free solutions, so it only improves
   * mass conservation and damps the spurious pressure component, without
   * changing the consistent solution.  For @f$\tau_C \ge 0@f$,
   * @f$C_K(u_h,u_h) = \int_K \tau_C (\nabla\cdot u_h)^2 \ge 0@f$.
   *
   * IMPORTANT: this is NOT a pressure-gradient (PSPG / "grad-p") term.  It is
   * NOT an inf-sup remedy and must not be used to stabilize equal-order
   * (e.g. @f$P_1/P_1@f$) pairs.  For Taylor-Hood it is exactly the right
   * additional stabilization; for equal-order one would instead need a PSPG
   * term @f$\int \tau_M \nabla q\cdot R_M@f$, which is unnecessary here.
   *
   * The local entry for vector basis functions @f$\Phi_a, \Psi_b@f$ is
   *
   * @f[
   *   C_{ba}
   *   =
   *   \int_K \tau_C\,(\nabla\cdot\Phi_a)(\nabla\cdot\Psi_b)\, dx ,
   * @f]
   *
   * with @f$\nabla\cdot\Phi = \sum_c (\hat\nabla\Phi\,J^{-1})_{cc}@f$.
   *
   * @tparam TrialFunction Trial velocity field type.
   * @tparam TestFunction Test velocity field type.
   * @tparam ProjectedTauC Projected grad-div parameter field type.
   */
  template <class TrialFunction, class TestFunction, class ProjectedTauC>
  class VMSGradDivBilinearIntegrator final
    : public Variational::LocalBilinearFormIntegratorBase<
        typename TrialFunction::ScalarType>
  {
    public:
      using ScalarType = typename TrialFunction::ScalarType;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructor.
       *
       * @param[in] u Trial velocity @f$u_h@f$.
       * @param[in] v Test velocity @f$v_h@f$.
       * @param[in] tauC Projected grad-div parameter @f$\tau_C@f$.
       */
      VMSGradDivBilinearIntegrator(
          const TrialFunction& u,
          const TestFunction& v,
          const ProjectedTauC& tauC)
        : Parent(u.getLeaf(), v.getLeaf()),
          m_u(u),
          m_v(v),
          m_tauC(tauC),
          m_polytope(nullptr)
      {}

      VMSGradDivBilinearIntegrator(
          const VMSGradDivBilinearIntegrator&) = default;

      VMSGradDivBilinearIntegrator(
          VMSGradDivBilinearIntegrator&&) = default;

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      VMSGradDivBilinearIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = &polytope;

        const size_t d = polytope.getDimension();
        const auto geometry = polytope.getGeometry();
        const auto idx = polytope.getIndex();

        const auto& trialFES = m_u.getFiniteElementSpace();
        const auto& testFES  = m_v.getFiniteElementSpace();

        const auto& trialFE = trialFES.getFiniteElement(d, idx);
        const auto& testFE  = testFES .getFiniteElement(d, idx);

        const size_t ntr = m_u.getDOFs(polytope);
        const size_t nte = m_v.getDOFs(polytope);
        const size_t vdim = trialFES.getVectorDimension();

        if (vdim == 0 || vdim != testFES.getVectorDimension())
          throw std::runtime_error("VMSGradDivBilinearIntegrator expects matching vector-valued spaces.");

        if (vdim != d)
          throw std::runtime_error("VMSGradDivBilinearIntegrator expects vector dimension to match the cell dimension.");

        if (ntr != trialFE.getCount() || nte != testFE.getCount())
          throw std::runtime_error("VMSGradDivBilinearIntegrator expects element-local DOFs to match finite-element basis count.");

        const size_t tauOrder =
          getFunctionOrder(m_tauC, polytope, size_t(0));
        const size_t qOrder =
          tauOrder
          + derivativeOrder(trialFE.getOrder())
          + derivativeOrder(testFE.getOrder());

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qOrder, geometry);

        const auto& q = polytope.getQuadrature(qf);

        m_mat.resize(
          static_cast<Eigen::Index>(nte),
          static_cast<Eigen::Index>(ntr));

        m_mat.setZero();

        ScalarType* A = m_mat.data();

        std::vector<ScalarType> trDiv(ntr);
        std::vector<ScalarType> teDiv(nte);

        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const Variational::IntegrationPoint ip(p, &qf, qp);

          const ScalarType wdet =
            static_cast<ScalarType>(qf.getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();
          const auto& rc = p.getReferenceCoordinates();

          const ScalarType tauC = m_tauC.getValue(ip);

          for (size_t a = 0; a < ntr; ++a)
            trDiv[a] = divergence(trialFE.getBasis(a), rc, Jinv, vdim, d);

          for (size_t b = 0; b < nte; ++b)
            teDiv[b] = divergence(testFE.getBasis(b), rc, Jinv, vdim, d);

          for (size_t b = 0; b < nte; ++b)
          {
            for (size_t a = 0; a < ntr; ++a)
            {
              A[b * ntr + a] += wdet * tauC * trDiv[a] * teDiv[b];
            }
          }
        }

        return *this;
      }

      ScalarType integrate(size_t trial, size_t test) final override
      {
        return m_mat(
          static_cast<Eigen::Index>(test),
          static_cast<Eigen::Index>(trial));
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      VMSGradDivBilinearIntegrator* copy() const noexcept final override
      {
        return new VMSGradDivBilinearIntegrator(*this);
      }

    private:
      /**
       * @brief Computes @f$\nabla\cdot\Phi = \sum_c (\hat\nabla\Phi\,J^{-1})_{cc}@f$.
       */
      template <class Basis, class JInv>
      static ScalarType divergence(
          const Basis& basis,
          const Math::SpatialVector<Real>& rc,
          const JInv& Jinv,
          size_t vdim,
          size_t d)
      {
        const auto Jref = basis.getJacobian()(rc);
        ScalarType div = 0;
        for (size_t c = 0; c < vdim; ++c)
        {
          for (size_t r = 0; r < d; ++r)
          {
            div +=
              Jref(static_cast<std::uint8_t>(c), static_cast<std::uint8_t>(r))
              * Jinv(static_cast<std::uint8_t>(r), static_cast<std::uint8_t>(c));
          }
        }
        return div;
      }

      static size_t derivativeOrder(size_t order) noexcept
      {
        return order == 0 ? 0 : order - 1;
      }

      template <class Function>
      static size_t getFunctionOrder(
          const Function& f,
          const Geometry::Polytope& polytope,
          size_t fallback)
      {
        if constexpr (requires { f.getOrder(polytope); })
        {
          const auto order = f.getOrder(polytope);
          if (order.has_value())
            return *order;
        }
        return fallback;
      }

      const TrialFunction& m_u;
      const TestFunction& m_v;
      const ProjectedTauC& m_tauC;

      const Geometry::Polytope* m_polytope;

      Eigen::Matrix<
        ScalarType,
        Eigen::Dynamic,
        Eigen::Dynamic,
        Eigen::RowMajor> m_mat;
  };

  /**
   * @brief Linear (projection-subtraction) part of the orthogonal grad-div /
   *        pressure-subscale stabilization.
   *
   * Pairs with VMSGradDivBilinearIntegrator to realize the ORTHOGONAL
   * (Codina OSS) div-u stabilization
   *
   * @f[
   *   \int_K \tau_C\,\bigl(P_h^\perp(\nabla\cdot u_h)\bigr)(\nabla\cdot v_h)
   *   =
   *   \underbrace{\int_K \tau_C(\nabla\cdot u_h)(\nabla\cdot v_h)}_{\text{bilinear}}
   *   -
   *   \underbrace{\int_K \tau_C\,\xi\,(\nabla\cdot v_h)}_{\text{this term}},
   * @f]
   *
   * where @f$\xi = \Pi_h(\nabla\cdot u_h)@f$ is the (lagged) finite-element
   * projection of the velocity divergence onto the pressure space, supplied as
   * @c ProjectedDiv.  Here @f$\tau_C = \rho\,\tau_2@f$ with the continuity
   * parameter @f$\tau_2 = (h/k^2)^2/(c_1\tau_1)@f$.  The pressure subscale is
   * @f$\tilde p = -\tau_C\,P_h^\perp(\nabla\cdot u_h)@f$.
   *
   * In the global residual this term appears with a MINUS sign:
   * @code - VMSGradDivLinearIntegrator(...) @endcode
   *
   * @tparam TestFunction Test velocity field type.
   * @tparam ProjectedDiv Projected divergence field type (scalar).
   * @tparam ProjectedTauC Projected grad-div parameter field type (scalar).
   */
  template <class TestFunction, class ProjectedDiv, class ProjectedTauC>
  class VMSGradDivLinearIntegrator final
    : public Variational::LinearFormIntegratorBase<
        typename TestFunction::ScalarType>
  {
    public:
      using ScalarType = typename TestFunction::ScalarType;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;

      VMSGradDivLinearIntegrator(
          const TestFunction& v,
          const ProjectedDiv& xi,
          const ProjectedTauC& tauC)
        : Parent(v.getLeaf()),
          m_v(v),
          m_xi(xi),
          m_tauC(tauC),
          m_polytope(nullptr)
      {}

      VMSGradDivLinearIntegrator(const VMSGradDivLinearIntegrator&) = default;
      VMSGradDivLinearIntegrator(VMSGradDivLinearIntegrator&&) = default;

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      VMSGradDivLinearIntegrator& setPolytope(
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

        if (vdim == 0)
          throw std::runtime_error("VMSGradDivLinearIntegrator expects a vector-valued test space.");
        if (vdim != d)
          throw std::runtime_error("VMSGradDivLinearIntegrator expects vector dimension to match the cell dimension.");
        if (nte != fe.getCount())
          throw std::runtime_error("VMSGradDivLinearIntegrator expects element-local DOFs to match finite-element basis count.");

        const size_t tauOrder = getFunctionOrder(m_tauC, polytope, size_t(0));
        const size_t xiOrder  = getFunctionOrder(m_xi, polytope, size_t(0));
        const size_t qOrder =
          tauOrder + xiOrder + derivativeOrder(fe.getOrder());

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qOrder, geometry);

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

          const ScalarType tauC = m_tauC.getValue(ip);
          const ScalarType xi   = m_xi.getValue(ip);
          const ScalarType c    = tauC * xi;

          for (size_t b = 0; b < nte; ++b)
          {
            const ScalarType divPsi = divergence(fe.getBasis(b), rc, Jinv, vdim, d);
            m_vec(static_cast<Eigen::Index>(b)) += wdet * c * divPsi;
          }
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

      VMSGradDivLinearIntegrator* copy() const noexcept final override
      {
        return new VMSGradDivLinearIntegrator(*this);
      }

    private:
      template <class Basis, class JInv>
      static ScalarType divergence(
          const Basis& basis,
          const Math::SpatialVector<Real>& rc,
          const JInv& Jinv,
          size_t vdim,
          size_t d)
      {
        const auto Jref = basis.getJacobian()(rc);
        ScalarType div = 0;
        for (size_t c = 0; c < vdim; ++c)
          for (size_t r = 0; r < d; ++r)
            div +=
              Jref(static_cast<std::uint8_t>(c), static_cast<std::uint8_t>(r))
              * Jinv(static_cast<std::uint8_t>(r), static_cast<std::uint8_t>(c));
        return div;
      }

      static size_t derivativeOrder(size_t order) noexcept
      {
        return order == 0 ? 0 : order - 1;
      }

      template <class Function>
      static size_t getFunctionOrder(
          const Function& f,
          const Geometry::Polytope& polytope,
          size_t fallback)
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
      const ProjectedDiv& m_xi;
      const ProjectedTauC& m_tauC;

      const Geometry::Polytope* m_polytope;

      Math::Vector<ScalarType> m_vec;
  };

  /**
   * @brief Cell integrator for the linear projected-VMS convective term.
   *
   * This class provides the explicit/source part paired with
   * VMSConvectionBilinearIntegrator. It represents the finite-element
   * projection and dynamic subscale correction that is subtracted from the
   * bilinear convective stabilization in the global residual.
   *
   * This assembles the right-hand-side-like term:
   *
   * @f[
   *   L_K(v_h)
   *   =
   *   \int_K
   *     \rho
   *     \left(
   *       \tau_K \rho u_h^{proj}
   *       +
   *       u_h^{\prime,n+1}
   *     \right)
   *     \cdot
   *     \bigl((\nabla v_h)u_h^n\bigr)
   *   \, dx .
   * @f]
   *
   * In the global variational expression this term is intended to appear with
   * a minus sign:
   *
   * @code
   * - VMSConvectionLinearIntegrator(...)
   * @endcode
   *
   * Therefore the total VMS contribution is:
   *
   * @f[
   *   A_K(u_h,v_h) - L_K(v_h).
   * @f]
   *
   * With the documented subscale update, this makes the stabilization act on
   * the projected fine-scale part of the convective acceleration. If
   * @f$u_h = u_h^n@f$ and the dynamic history is zero, the projected component
   * cancels the matching part of the bilinear term, leaving only the
   * orthogonal convective correction.
   *
   * @tparam TestFunction Test velocity field type.
   * @tparam OldSubScale Dynamic subscale field type.
   * @tparam OldVelocity Frozen convective velocity field type.
   * @tparam ProjectedVelocity Projected convective acceleration type.
   * @tparam ProjectedTau Projected stabilization parameter field type.
   */
  template <class TestFunction, class OldSubScale, class OldVelocity, class ProjectedVelocity, class ProjectedTau>
  class VMSConvectionLinearIntegrator final
    : public Variational::LinearFormIntegratorBase<
        typename TestFunction::ScalarType>
  {
    public:
      using ScalarType = typename TestFunction::ScalarType;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;

      /**
       * @brief Constructor.
       *
       * @param[in] v Test velocity @f$v_h@f$.
       * @param[in] subOld Dynamic subscale @f$u_h^{\prime,n+1}@f$.
       * @param[in] uOld Frozen convective velocity @f$u_h^n@f$.
       * @param[in] uProj Projected convective acceleration
       *                  @f$\Pi_h[(\nabla u_h^n)u_h^n]@f$.
       * @param[in] rho Fluid density @f$\rho@f$.
       * @param[in] tau Projected stabilization parameter @f$\tau_K@f$.
       */
      VMSConvectionLinearIntegrator(
          const TestFunction& v,
          const OldSubScale& subOld,
          const OldVelocity& uOld,
          const ProjectedVelocity& uProj,
          const ProjectedTau& tau,
          ScalarType rho,
          ScalarType dt)
        : Parent(v.getLeaf()),
          m_v(v),
          m_sub(subOld),
          m_uOld(uOld),
          m_uProj(uProj),
          m_tau(tau),
          m_polytope(nullptr),
          m_rho(rho),
          m_dt(dt)
      {}

      VMSConvectionLinearIntegrator(
          const VMSConvectionLinearIntegrator&) = default;

      VMSConvectionLinearIntegrator(
          VMSConvectionLinearIntegrator&&) = default;

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return *m_polytope;
      }

      /**
       * @brief Builds the local vector on one cell.
       */
      VMSConvectionLinearIntegrator& setPolytope(
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

        if (vdim == 0)
          throw std::runtime_error("VMSConvectionLinearIntegrator expects a vector-valued test space.");

        if (nte != fe.getCount())
          throw std::runtime_error("VMSConvectionLinearIntegrator expects element-local DOFs to match finite-element basis count.");

        /*
         * Conservative quadrature choice for
         *
         *   rho (tau rho uProj + subScale) . ((grad Psi) uOld).
         *
         * The projected fields can be polynomial, while tau is usually a
         * projected representation of a non-polynomial expression.
         */
        const size_t uOldOrder =
          getFunctionOrder(m_uOld, polytope, fe.getOrder());
        const size_t tauOrder =
          getFunctionOrder(m_tau, polytope, size_t(0));
        const size_t projOrder =
          getFunctionOrder(m_uProj, polytope, uOldOrder);
        const size_t subOrder =
          getFunctionOrder(m_sub, polytope, projOrder);
        const size_t coefficientOrder =
          std::max(tauOrder + projOrder, subOrder);
        const size_t qOrder =
          coefficientOrder + derivativeOrder(fe.getOrder()) + uOldOrder;

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qOrder, geometry);

        const auto& q = polytope.getQuadrature(qf);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        std::vector<Math::SpatialVector<ScalarType>> teDir(nte);
        Math::SpatialVector<ScalarType> coefficient(static_cast<std::uint8_t>(vdim));

        for (auto& t : teDir)
          t.resize(static_cast<std::uint8_t>(vdim));

        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);
          const Variational::IntegrationPoint ip(p, &qf, qp);

          const ScalarType wdet =
            static_cast<ScalarType>(qf.getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();
          const auto& rc = p.getReferenceCoordinates();

          const auto uOld  = m_uOld.getValue(ip);
          const auto uProj = m_uProj.getValue(ip);
          const auto subScale = m_sub.getValue(ip);
          const auto tau = m_tau.getValue(ip);

          if (uOld.size() != d)
            throw std::runtime_error("VMSConvectionLinearIntegrator expects the frozen velocity dimension to match the cell dimension.");

          /*
           * Directional derivative of test basis along frozen velocity:
           *
           *   teDir[b] = (grad Psi_b) uOld.
           */
          for (size_t b = 0; b < nte; ++b)
            fillDirectionalDerivative(teDir[b], fe.getBasis(b), rc, Jinv, uOld);

          for (size_t c = 0; c < vdim; ++c)
          {
            coefficient(static_cast<std::uint8_t>(c)) =
              tau * m_rho * (uProj(static_cast<std::uint8_t>(c))
              + 1. / m_dt * subScale(static_cast<std::uint8_t>(c)));
          }

          for (size_t b = 0; b < nte; ++b)
          {
            /*
             * Local vector-basis entry:
             *
             *   int_K rho
             *     (tau rho uProj + sub)
             *     ·
             *     ((grad Psi_b) uOld).
             *
             * The global expression should subtract this integrator.
             */
            m_vec(static_cast<Eigen::Index>(b)) +=
              wdet * m_rho * Math::dot(coefficient, teDir[b]);
          }
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

      VMSConvectionLinearIntegrator* copy() const noexcept final override
      {
        return new VMSConvectionLinearIntegrator(*this);
      }

    private:
      /**
       * @brief Computes @f$(\nabla \Psi)uOld@f$ for one vector basis.
       *
       * The vector finite-element basis supplies the reference Jacobian
       * @f$\hat\nabla\Psi@f$. Physical derivatives are obtained by
       * @f$\nabla_x\Psi = \hat\nabla\Psi J^{-1}@f$.
       */
      template <class Basis, class JInv, class Vector>
      static void fillDirectionalDerivative(
          Math::SpatialVector<ScalarType>& out,
          const Basis& basis,
          const Math::SpatialVector<Real>& rc,
          const JInv& Jinv,
          const Vector& uOld)
      {
        const auto Jref = basis.getJacobian()(rc);
        const size_t vdim = out.size();
        const size_t d = uOld.size();

        out.setZero();
        for (size_t c = 0; c < vdim; ++c)
        {
          for (size_t j = 0; j < d; ++j)
          {
            ScalarType gradPhys = 0;
            for (size_t r = 0; r < d; ++r)
            {
              gradPhys +=
                Jref(static_cast<std::uint8_t>(c), static_cast<std::uint8_t>(r))
                * Jinv(static_cast<std::uint8_t>(r), static_cast<std::uint8_t>(j));
            }
            out(static_cast<std::uint8_t>(c)) +=
              gradPhys * uOld(static_cast<std::uint8_t>(j));
          }
        }
      }

      static size_t derivativeOrder(size_t order) noexcept
      {
        return order == 0 ? 0 : order - 1;
      }

      template <class Function>
      static size_t getFunctionOrder(
          const Function& f,
          const Geometry::Polytope& polytope,
          size_t fallback)
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
      const OldSubScale& m_sub;
      const OldVelocity& m_uOld;
      const ProjectedVelocity& m_uProj;
      const ProjectedTau& m_tau;

      const Geometry::Polytope* m_polytope;
      ScalarType m_rho;
      ScalarType m_dt;

      Math::Vector<ScalarType> m_vec;

  };
}

#endif
