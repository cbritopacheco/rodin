#ifndef RODIN_EXAMPLES_HEART_VMSCONVECTIONINTEGRATOR_H
#define RODIN_EXAMPLES_HEART_VMSCONVECTIONINTEGRATOR_H

/**
 * @file VMSConvectionIntegrator.h
 * @brief Custom VMS convective stabilization integrators for the coronary example.
 *
 * This file defines local bilinear and linear form integrators used by the
 * coupled 0D--3D coronary artery example. They implement the VMS terms
 * associated with the projected convective residual:
 * @f[
 *   \tau_K
 *   \left(
 *     (\nabla u_h)u_h^{old}
 *     -
 *     u_h^{proj}
 *   \right)
 *   \cdot
 *   ((\nabla v_h)u_h^{old}) .
 * @f]
 *
 * The implementation is intentionally example-specific:
 *
 * - three-dimensional cells only,
 * - quadratic @f$H^1@f$ velocity basis functions,
 * - diagonal component-wise vector assembly,
 * - cell-wise stabilization length
 *   @f$h_K = |K|^{1/\dim K}@f$.
 */

#include <cassert>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Variational/H1/H1Element.h>

namespace Rodin::Examples::Heart
{
  /**
   * @brief Cell integrator for the bilinear VMS convective stabilization term.
   *
   * This class represents the cell-wise contribution
   * @f[
   *   \int_K \tau_K
   *     \bigl((\nabla u_h) u_h^{old}\bigr)
   *     \cdot
   *     \bigl((\nabla v_h) u_h^{old}\bigr)
   *   \, dx .
   * @f]
   *
   * Here @f$u_h@f$ is the trial velocity, @f$v_h@f$ is the test velocity,
   * @f$u_h^{old}@f$ is the frozen convective velocity, and @f$\tau_K@f$ is
   * computed locally as
   * @f[
   *   \tau_K =
   *   \left(
   *     c_1 \frac{\mu_h}{h_K^2}
   *     +
   *     c_2 \frac{\lVert u_h^{old} \rVert}{h_K}
   *   \right)^{-1}.
   * @f]
   *
   * The local matrix entry is assembled as
   * @f[
   *   A_{bi}
   *   =
   *   \int_K
   *     \tau_K
   *     \bigl(\nabla \phi_i \cdot u_h^{old}\bigr)
   *     \bigl(\nabla \psi_b \cdot u_h^{old}\bigr)
   *   \, dx ,
   * @f]
   * copied component-wise on the diagonal velocity blocks.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a local bilinear
   * form integrator:
   * @f[
   * \dfrac
   * {
   *   \vdash
   *   \int_K \tau_K
   *   ((\nabla u_h)u_h^{old})
   *   \cdot
   *   ((\nabla v_h)u_h^{old})
   *   \, dx
   *   :
   *   \texttt{LocalBilinearFormIntegrator}
   * }
   * {
   *   \vdash u_h : H^1_2(\Omega)^3,
   *   \quad
   *   \vdash v_h : H^1_2(\Omega)^3
   * }
   * @f]
   *
   * @tparam TrialFunction Trial velocity shape-function type.
   * @tparam TestFunction Test velocity shape-function type.
   * @tparam OldVelocity Frozen convective velocity field type.
   * @tparam Viscosity Effective viscosity field type.
   */
  template <class TrialFunction, class TestFunction, class OldVelocity, class Viscosity>
  class VMSConvectionBilinearIntegrator final
    : public Variational::LocalBilinearFormIntegratorBase<
        typename TrialFunction::ScalarType>
  {
    public:
      using ScalarType = typename TrialFunction::ScalarType;
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;

      VMSConvectionBilinearIntegrator(
          const TrialFunction& u,
          const TestFunction& v,
          const OldVelocity& uOld,
          const Viscosity& mu,
          ScalarType c1 = 4.0,
          ScalarType c2 = 2.0)
        : Parent(u.getLeaf(), v.getLeaf()),
          m_u(u),
          m_v(v),
          m_uOld(uOld),
          m_mu(mu),
          m_c1(c1),
          m_c2(c2),
          m_polytope(nullptr)
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

        const size_t qOrder = trialFE.getOrder() + testFE.getOrder() - 2;
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
        std::vector<ScalarType> trDir(ntrS);
        std::vector<ScalarType> teDir(nteS);

        for (auto& g : Gtr)
          g.resize(static_cast<std::uint8_t>(d));
        for (auto& g : Gte)
          g.resize(static_cast<std::uint8_t>(d));

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

          const ScalarType tau =
            ScalarType(1)
            / (m_c1 * mu / (hK * hK) + m_c2 * speed / hK);

          for (size_t a = 0; a < ntrS; ++a)
            trDir[a] = Math::dot(Gtr[a], uOld);

          for (size_t b = 0; b < nteS; ++b)
            teDir[b] = Math::dot(Gte[b], uOld);

          for (size_t b = 0; b < nteS; ++b)
          {
            for (size_t a = 0; a < ntrS; ++a)
            {
              const ScalarType kij =
                wdet * tau * trDir[a] * teDir[b];

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
  };

  /**
   * @brief Cell integrator for the linear projected VMS convective term.
   *
   * This class represents the cell-wise contribution
   * @f[
   *   \int_K \tau_K
   *     u_h^{proj}
   *     \cdot
   *     \bigl((\nabla v_h)u_h^{old}\bigr)
   *   \, dx .
   * @f]
   *
   * In the residual used by the coronary example, this term is subtracted:
   * @f[
   *   -
   *   \int_K \tau_K
   *     u_h^{proj}
   *     \cdot
   *     \bigl((\nabla v_h)u_h^{old}\bigr)
   *   \, dx .
   * @f]
   * The sign is therefore expected to be supplied by the variational expression,
   * for example with @code -vmsLinear @endcode.
   *
   * The stabilization parameter is computed locally as
   * @f[
   *   \tau_K =
   *   \left(
   *     c_1 \frac{\mu_h}{h_K^2}
   *     +
   *     c_2 \frac{\lVert u_h^{old} \rVert}{h_K}
   *   \right)^{-1}.
   * @f]
   *
   * The local vector entry is assembled as
   * @f[
   *   F_b
   *   =
   *   \int_K
   *     \tau_K
   *     u_{h,c}^{proj}
   *     \bigl(\nabla \psi_b \cdot u_h^{old}\bigr)
   *   \, dx ,
   * @f]
   * for each velocity component @f$c@f$.
   *
   * Judgement
   * ---------
   *
   * The following judgement specifies that the expression is a local linear
   * form integrator:
   * @f[
   * \dfrac
   * {
   *   \vdash
   *   \int_K \tau_K
   *   u_h^{proj}
   *   \cdot
   *   ((\nabla v_h)u_h^{old})
   *   \, dx
   *   :
   *   \texttt{LinearFormIntegrator}
   * }
   * {
   *   \vdash v_h : H^1_2(\Omega)^3
   * }
   * @f]
   *
   * @tparam TestFunction Test velocity shape-function type.
   * @tparam OldVelocity Frozen convective velocity field type.
   * @tparam ProjectedVelocity Projected convective acceleration field type.
   * @tparam Viscosity Effective viscosity field type.
   */
  template <class TestFunction, class OldVelocity, class ProjectedVelocity, class Viscosity>
  class VMSConvectionLinearIntegrator final
    : public Variational::LinearFormIntegratorBase<
        typename TestFunction::ScalarType>
  {
    public:
      using ScalarType = typename TestFunction::ScalarType;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;

      VMSConvectionLinearIntegrator(
          const TestFunction& v,
          const OldVelocity& uOld,
          const ProjectedVelocity& uProj,
          const Viscosity& mu,
          ScalarType c1 = 4.0,
          ScalarType c2 = 2.0)
        : Parent(v.getLeaf()),
          m_v(v),
          m_uOld(uOld),
          m_uProj(uProj),
          m_mu(mu),
          m_c1(c1),
          m_c2(c2),
          m_polytope(nullptr)
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

        const size_t qOrder = fe.getOrder() - 1;
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(qOrder, geometry);

        const auto& q = polytope.getQuadrature(qf);
        const auto& tab = scalarFE.getTabulation(qf);

        m_vec.resize(static_cast<Eigen::Index>(nte));
        m_vec.setZero();

        std::vector<Math::SpatialVector<ScalarType>> Gte(nteS);
        std::vector<ScalarType> teDir(nteS);

        for (auto& g : Gte)
          g.resize(static_cast<std::uint8_t>(d));

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

          const auto uOld  = m_uOld.getValue(p);
          const auto uProj = m_uProj.getValue(p);
          const ScalarType mu = m_mu.getValue(p);

          const ScalarType speed = std::sqrt(Math::dot(uOld, uOld));

          const ScalarType tau =
            ScalarType(1)
            / (m_c1 * mu / (hK * hK) + m_c2 * speed / hK);

          for (size_t b = 0; b < nteS; ++b)
            teDir[b] = Math::dot(Gte[b], uOld);

          for (size_t b = 0; b < nteS; ++b)
          {
            for (size_t c = 0; c < vdim; ++c)
            {
              const size_t row = b * vdim + c;

              m_vec(static_cast<Eigen::Index>(row)) +=
                wdet * tau * uProj[c] * teDir[b];
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
      const OldVelocity& m_uOld;
      const ProjectedVelocity& m_uProj;
      const Viscosity& m_mu;

      ScalarType m_c1;
      ScalarType m_c2;

      const Geometry::Polytope* m_polytope;

      Math::Vector<ScalarType> m_vec;
  };
}

#endif
