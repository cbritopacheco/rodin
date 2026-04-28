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
 * header. It is computed in the driver and passed as @c OldSubScale:
 *
 * @f[
 *   u_h^{\prime,n+1}
 *   =
 *   \tau_K \rho
 *   \left(
 *     (\nabla u_h^n)u_h^n
 *     -
 *     \Pi_h[(\nabla u_h^n)u_h^n]
 *   \right)
 *   +
 *   \tau_K \frac{\rho}{\Delta t} u_h^{\prime,n}.
 * @f]
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
 * - 3D cells only,
 * - vector-valued quadratic H1 velocity fields,
 * - diagonal component-wise vector assembly,
 * - explicit frozen convective velocity @f$u_h^n@f$,
 * - explicit effective viscosity @f$\mu_h@f$,
 * - local stabilization length @f$h_K = |K|^{1/3}@f$.
 */

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

namespace Rodin::Examples::Heart
{
  /**
   * @brief Cell integrator for the bilinear projected-VMS convective term.
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
   * For scalar basis functions @f$\phi_a@f$ and @f$\psi_b@f$, the local
   * scalar block is:
   *
   * @f[
   *   A_{ba}
   *   =
   *   \int_K
   *     \tau_K \rho^2
   *     \left(\nabla \phi_a \cdot u_h^n\right)
   *     \left(\nabla \psi_b \cdot u_h^n\right)
   *   \, dx .
   * @f]
   *
   * This scalar block is copied onto each diagonal velocity component block.
   * Cross-component coupling is not assembled here because
   * @f$((\nabla u_h)u_h^n)_c@f$ depends only on the @f$c@f$-th velocity
   * component basis functions.
   *
   * @tparam TrialFunction Trial velocity field type.
   * @tparam TestFunction Test velocity field type.
   * @tparam OldVelocity Frozen convective velocity field type.
   * @tparam Viscosity Explicit effective viscosity field type.
   */
  template <class TrialFunction, class TestFunction, class OldVelocity, class Viscosity>
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
       * @param[in] mu Explicit viscosity @f$\mu_h@f$.
       * @param[in] rho Fluid density @f$\rho@f$.
       * @param[in] dt Time step @f$\Delta t@f$.
       * @param[in] c1 Diffusive stabilization constant.
       * @param[in] c2 Convective stabilization constant.
       * @param[in] vmsScale Optional empirical scaling of @f$\tau_K@f$.
       */
      VMSConvectionBilinearIntegrator(
          const TrialFunction& u,
          const TestFunction& v,
          const OldVelocity& uOld,
          const Viscosity& mu,
          ScalarType rho,
          ScalarType dt,
          ScalarType c1 = 4.0,
          ScalarType c2 = 2.0,
          ScalarType vmsScale = 0.05)
        : Parent(u.getLeaf(), v.getLeaf()),
          m_u(u),
          m_v(v),
          m_uOld(uOld),
          m_mu(mu),
          m_c1(c1),
          m_c2(c2),
          m_polytope(nullptr),
          m_rho(rho),
          m_dt(dt),
          m_vmsScale(vmsScale)
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

        if (d != 3)
          throw std::runtime_error("VMSConvectionBilinearIntegrator expects 3D cells.");

        const auto& trialFES = m_u.getFiniteElementSpace();
        const auto& testFES  = m_v.getFiniteElementSpace();

        const auto& trialFE = trialFES.getFiniteElement(d, idx);
        const auto& testFE  = testFES .getFiniteElement(d, idx);

        const size_t ntr = m_u.getDOFs(polytope);
        const size_t nte = m_v.getDOFs(polytope);

        constexpr size_t KTrial = 2;
        constexpr size_t KTest  = 2;

        const Variational::H1Element<KTrial, ScalarType> trialScalarFE(geometry);
        const Variational::H1Element<KTest,  ScalarType> testScalarFE(geometry);

        const size_t ntrS = trialScalarFE.getCount();
        const size_t nteS = testScalarFE.getCount();

        assert(ntrS > 0 && nteS > 0);
        assert(ntr % ntrS == 0);
        assert(nte % nteS == 0);

        const size_t vdim = ntr / ntrS;
        assert(vdim == nte / nteS);

        /*
         * Conservative quadrature choice.
         *
         * The integrand contains products of directional derivatives and
         * non-polynomial coefficients through tau_K if mu/uOld vary in space.
         * This order is intentionally not minimal.
         */
        const size_t qOrder = trialFE.getOrder() + testFE.getOrder();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qOrder, geometry);

        const auto& q = polytope.getQuadrature(qf);

        const auto& trTab = trialScalarFE.getTabulation(qf);
        const auto& teTab = testScalarFE .getTabulation(qf);

        m_mat.resize(
          static_cast<Eigen::Index>(nte),
          static_cast<Eigen::Index>(ntr));

        m_mat.setZero();

        ScalarType* A = m_mat.data();

        std::vector<Math::SpatialVector<ScalarType>> Gtr(ntrS);
        std::vector<Math::SpatialVector<ScalarType>> Gte(nteS);

        Math::SpatialVector<ScalarType> trDir(ntrS);
        Math::SpatialVector<ScalarType> teDir(nteS);

        for (auto& g : Gtr)
          g.resize(static_cast<std::uint8_t>(d));

        for (auto& g : Gte)
          g.resize(static_cast<std::uint8_t>(d));

        /*
         * Cell length scale:
         *
         *   h_K = |K|^{1/3}.
         */
        const ScalarType hK =
          std::pow(
            polytope.getMeasure(),
            ScalarType(1) / static_cast<ScalarType>(d));

        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);

          const ScalarType wdet =
            static_cast<ScalarType>(qf.getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          fillPhysicalGradients3D(qp, Jinv, trTab, Gtr);
          fillPhysicalGradients3D(qp, Jinv, teTab, Gte);

          const auto uOld = m_uOld.getValue(p);
          const ScalarType mu = m_mu.getValue(p);

          const ScalarType speed = std::sqrt(Math::dot(uOld, uOld));

          /*
           * Dynamic VMS time/transport/diffusion stabilization scale.
           */
          const ScalarType invTau =
            std::sqrt(
                Math::pow2(ScalarType(2) * m_rho / m_dt)
              + Math::pow2(m_c2 * m_rho * speed / hK)
              + Math::pow2(m_c1 * mu / (hK * hK)));

          const ScalarType tau = m_vmsScale / invTau;

          /*
           * Directional derivatives along frozen velocity:
           *
           *   trDir[a] = grad phi_a · uOld,
           *   teDir[b] = grad psi_b · uOld.
           */
          for (size_t a = 0; a < ntrS; ++a)
            trDir[a] = Math::dot(Gtr[a], uOld);

          for (size_t b = 0; b < nteS; ++b)
            teDir[b] = Math::dot(Gte[b], uOld);

          for (size_t b = 0; b < nteS; ++b)
          {
            for (size_t a = 0; a < ntrS; ++a)
            {
              /*
               * Local scalar entry:
               *
               *   int_K tau rho^2
               *     (grad phi_a · uOld)
               *     (grad psi_b · uOld).
               */
              const ScalarType kij =
                wdet * tau * m_rho * trDir[a] * m_rho * teDir[b];

              /*
               * Copy the scalar entry on each velocity diagonal block.
               */
              for (size_t c = 0; c < vdim; ++c)
              {
                const size_t row = b * vdim + c;
                const size_t col = a * vdim + c;

                A[row * ntr + col] += kij;
              }
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
       * @brief Maps reference gradients to physical gradients in 3D.
       *
       * If @f$\hat \nabla \phi@f$ is the reference gradient and @f$J@f$ is
       * the geometric Jacobian, then:
       *
       * @f[
       *   \nabla_x \phi = J^{-T} \hat \nabla \phi .
       * @f]
       */
      template <class Tabulation, class JInv>
      static void fillPhysicalGradients3D(
          size_t qp,
          const JInv& Jinv,
          const Tabulation& tab,
          std::vector<Math::SpatialVector<ScalarType>>& G)
      {
        const ScalarType j00 = Jinv(0,0), j10 = Jinv(1,0), j20 = Jinv(2,0);
        const ScalarType j01 = Jinv(0,1), j11 = Jinv(1,1), j21 = Jinv(2,1);
        const ScalarType j02 = Jinv(0,2), j12 = Jinv(1,2), j22 = Jinv(2,2);

        for (size_t a = 0; a < G.size(); ++a)
        {
          const auto g = tab.getGradient(qp, a);

          const ScalarType gx = g[0];
          const ScalarType gy = g[1];
          const ScalarType gz = g[2];

          G[a][0] = j00 * gx + j10 * gy + j20 * gz;
          G[a][1] = j01 * gx + j11 * gy + j21 * gz;
          G[a][2] = j02 * gx + j12 * gy + j22 * gz;
        }
      }

      const TrialFunction& m_u;
      const TestFunction& m_v;
      const OldVelocity& m_uOld;
      const Viscosity& m_mu;

      ScalarType m_c1;
      ScalarType m_c2;

      const Geometry::Polytope* m_polytope;

      Eigen::Matrix<
        ScalarType,
        Eigen::Dynamic,
        Eigen::Dynamic,
        Eigen::RowMajor> m_mat;

      ScalarType m_rho;
      ScalarType m_dt;
      ScalarType m_vmsScale;
  };

  /**
   * @brief Cell integrator for the linear projected-VMS convective term.
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
   * @tparam TestFunction Test velocity field type.
   * @tparam OldSubScale Dynamic subscale field type.
   * @tparam OldVelocity Frozen convective velocity field type.
   * @tparam ProjectedVelocity Projected convective acceleration type.
   * @tparam Viscosity Explicit effective viscosity field type.
   */
  template <
    class TestFunction,
    class OldSubScale,
    class OldVelocity,
    class ProjectedVelocity,
    class Viscosity>
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
       * @param[in] mu Explicit viscosity @f$\mu_h@f$.
       * @param[in] rho Fluid density @f$\rho@f$.
       * @param[in] dt Time step @f$\Delta t@f$.
       * @param[in] c1 Diffusive stabilization constant.
       * @param[in] c2 Convective stabilization constant.
       * @param[in] vmsScale Optional empirical scaling of @f$\tau_K@f$.
       */
      VMSConvectionLinearIntegrator(
          const TestFunction& v,
          const OldSubScale& subOld,
          const OldVelocity& uOld,
          const ProjectedVelocity& uProj,
          const Viscosity& mu,
          ScalarType rho,
          ScalarType dt,
          ScalarType c1 = 4.0,
          ScalarType c2 = 2.0,
          ScalarType vmsScale = 0.05)
        : Parent(v.getLeaf()),
          m_v(v),
          m_sub(subOld),
          m_uOld(uOld),
          m_uProj(uProj),
          m_mu(mu),
          m_c1(c1),
          m_c2(c2),
          m_polytope(nullptr),
          m_rho(rho),
          m_dt(dt),
          m_vmsScale(vmsScale)
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

        if (d != 3)
          throw std::runtime_error("VMSConvectionLinearIntegrator expects 3D cells.");

        const auto& fes = m_v.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);

        const size_t nte = m_v.getDOFs(polytope);

        constexpr size_t KTest = 2;

        const Variational::H1Element<KTest, ScalarType> scalarFE(geometry);

        const size_t nteS = scalarFE.getCount();

        assert(nteS > 0);
        assert(nte % nteS == 0);

        const size_t vdim = nte / nteS;

        /*
         * Conservative quadrature choice. The coefficients are not generally
         * polynomial because tau_K depends on uOld and mu.
         */
        const size_t qOrder = fe.getOrder();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qOrder, geometry);

        const auto& q = polytope.getQuadrature(qf);
        const auto& tab = scalarFE.getTabulation(qf);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        std::vector<Math::SpatialVector<ScalarType>> Gte(nteS);
        Math::SpatialVector<ScalarType> teDir(nteS);

        for (auto& g : Gte)
          g.resize(static_cast<std::uint8_t>(d));

        /*
         * Cell length scale:
         *
         *   h_K = |K|^{1/3}.
         */
        const ScalarType hK =
          std::pow(
            polytope.getMeasure(),
            ScalarType(1) / static_cast<ScalarType>(d));

        for (size_t qp = 0; qp < q.getSize(); ++qp)
        {
          const auto& p = q.getPoint(qp);

          const ScalarType wdet =
            static_cast<ScalarType>(qf.getWeight(qp) * p.getDistortion());

          const auto Jinv = p.getJacobianInverse();

          fillPhysicalGradients3D(qp, Jinv, tab, Gte);

          const auto uOld    = m_uOld.getValue(p);
          const auto uProj   = m_uProj.getValue(p);
          const auto sub     = m_sub.getValue(p);
          const ScalarType mu = m_mu.getValue(p);

          const ScalarType speed = std::sqrt(Math::dot(uOld, uOld));

          /*
           * Same stabilization parameter as in the bilinear integrator.
           */
          const ScalarType invTau =
            std::sqrt(
                Math::pow2(ScalarType(2) * m_rho / m_dt)
              + Math::pow2(m_c2 * m_rho * speed / hK)
              + Math::pow2(m_c1 * mu / (hK * hK)));

          const ScalarType tau = m_vmsScale / invTau;

          /*
           * Directional derivative of test basis along frozen velocity:
           *
           *   teDir[b] = grad psi_b · uOld.
           */
          for (size_t b = 0; b < nteS; ++b)
            teDir[b] = Math::dot(Gte[b], uOld);

          for (size_t b = 0; b < nteS; ++b)
          {
            for (size_t c = 0; c < vdim; ++c)
            {
              const size_t row = b * vdim + c;

              /*
               * Local vector entry:
               *
               *   int_K rho
               *     (tau rho uProj_c + sub_c)
               *     (grad psi_b · uOld).
               *
               * The global expression should subtract this integrator.
               */
              m_vec(static_cast<Eigen::Index>(row)) +=
                wdet
                * m_rho
                * (tau * m_rho * uProj[c] + sub[c])
                * teDir[b];
            }
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
       * @brief Maps reference gradients to physical gradients in 3D.
       *
       * If @f$\hat \nabla \psi@f$ is the reference gradient and @f$J@f$ is
       * the geometric Jacobian, then:
       *
       * @f[
       *   \nabla_x \psi = J^{-T} \hat \nabla \psi .
       * @f]
       */
      template <class Tabulation, class JInv>
      static void fillPhysicalGradients3D(
          size_t qp,
          const JInv& Jinv,
          const Tabulation& tab,
          std::vector<Math::SpatialVector<ScalarType>>& G)
      {
        const ScalarType j00 = Jinv(0,0), j10 = Jinv(1,0), j20 = Jinv(2,0);
        const ScalarType j01 = Jinv(0,1), j11 = Jinv(1,1), j21 = Jinv(2,1);
        const ScalarType j02 = Jinv(0,2), j12 = Jinv(1,2), j22 = Jinv(2,2);

        for (size_t b = 0; b < G.size(); ++b)
        {
          const auto g = tab.getGradient(qp, b);

          const ScalarType gx = g[0];
          const ScalarType gy = g[1];
          const ScalarType gz = g[2];

          G[b][0] = j00 * gx + j10 * gy + j20 * gz;
          G[b][1] = j01 * gx + j11 * gy + j21 * gz;
          G[b][2] = j02 * gx + j12 * gy + j22 * gz;
        }
      }

      const TestFunction& m_v;
      const OldSubScale& m_sub;
      const OldVelocity& m_uOld;
      const ProjectedVelocity& m_uProj;
      const Viscosity& m_mu;

      ScalarType m_c1;
      ScalarType m_c2;

      const Geometry::Polytope* m_polytope;

      Math::Vector<ScalarType> m_vec;

      ScalarType m_rho;
      ScalarType m_dt;
      ScalarType m_vmsScale;
  };
}

#endif
