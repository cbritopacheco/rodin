/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIR_H
#define RODIN_ADAPTATION_WNGIR_H

//
// Welsch natural-gradient interface registration (WNGIR).
//
// Classify-then-displace: the classifier provides a topological-
// reference skeleton Γ_ψ,h; the target geometry is the level set φ
// with ∇φ; the displacement u ∈ V_{h,0} (vector Pk, zero trace on ∂Ω)
// moves the skeleton onto φ = 0.
//
// Per frame (cold start u = 0, no warm start, no σ continuation):
//   σ = max(3h, quantile_90(|φ| on Γ_ψ,h at u = 0)),  fixed.
//
// Per nonlinear iteration k:
//   1. Welsch data at skeleton quadrature points X_q:
//        y_q = X_q + u^k(X_q),  r_q = φ(y_q),  g_q = ∇φ(y_q),
//        ω_q = exp(−r_q²/σ²).
//   2. Assemble the negative first variation
//        FΓ(v) = −∫_Γ ω_q r_q g_q·trace(v) dS.
//   3. H¹ bulk metric plus residual-stabilized observation metric
//        M(v,z) = γ_M ∫ v·z + γ_H ℓ_M² ∫ ∇v:∇z.
//        Mobs(v,z) = γ_obs ∫_Γ (|g|² + ε_g + r²/σ²) v·z dS.
//   4. Linearised admissibility at quadrature/validation points (K,q):
//        F = I + ∇u^k,  j = det F,  Q = ‖F‖²_F/(d·det F^{2/d}),
//        a_j(v) = j F^{-T}:∇v,
//        dQdF = (2/d) j^{−2/d}(F − (F:F)/d · F^{-T}),
//        a_Q(v) = dQdF:∇v,
//        s_j(v) = j − j_safe + a_j(v),  s_Q(v) = Q_max − Q − a_Q(v).
//   5. One-SPD active log-barrier metric:
//        B(s) = −log(s/s0) + s/s0 − 1 on (0, s0), else 0,
//        B'(s) = −1/s + 1/s0,  B''(s) = 1/s²,
//      assembled as Gauss–Newton pullback:
//        (M + Mobs + K_adm) v = FΓ.
//   6. Optimal 1-D step rescale β = ⟨d,v⟩_Γ/⟨v,v⟩_Γ ∈ [1, βmax].
//   7. Nonlinear line search on the TRUE geometry:
//      accept the largest α ∈ (0,1] with j > j_ls, Q < Q_max at all
//      validation points; optionally also E_Welsch decrease.
//   8. u^{k+1} = u^k + α β v.
//   9. Stop on: small step, energy stagnation, active RMS/sup below
//      tolerance, iteration budget, or line-search failure.
//
// Order-generic: the displacement space may be vector P1 or vector
// H1 of order k ≥ 2. All basis values/gradients are queried through
// the FES finite-element API and tabulated once per frame at cell and
// facet quadrature points; per-iteration work is pure coefficient
// contractions. The mesh geometry itself is straight (affine cells);
// curved isoparametric geometry would require per-point Jacobians in
// the tabulation pass (single localised change).
//
// Welsch is for partial matching / topology mismatch: the active
// subset (ω ≥ ω_min) has geometric control through the level-set
// residual; the inactive subset is topological-outlier skeleton. The
// residual is φ(X + u(X)) on Γ_ψ,h — ψ never enters.
//

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include <Rodin/Math.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Rodin::Adaptation
{
  struct WNGIRParameters
  {
    Real h = 0;                 ///< reference mesh size (required).
    Real gammaM = 0;            ///< L² weight; ≤0 ⇒ 1/h.
    Real gammaH = 0;            ///< H¹ weight; ≤0 ⇒ 1/h.
    Real ellM = 0;              ///< Sobolev length; ≤0 ⇒ 3h.
    Real gammaObs = 1;          ///< surface observation metric weight.
    bool residualStabilizedObservationMetric = true;
    Real gammaJ = 1;            ///< j-barrier weight.
    Real gammaQ = 1;            ///< Q-barrier weight.
    Real jSafe = 1e-2;          ///< barrier floor on j.
    Real qMax = 10;             ///< barrier + line-search ceiling on Q.
    Real s0J = 0.25;            ///< j-barrier activation width.
    Real s0Q = 2;               ///< Q-barrier activation width.
    /// Hinge quality regularizer (energy, on the RHS). Penalizes only
    /// cells whose relative distortion exceeds a *good-quality*
    /// threshold qStar (chosen ≪ qMax):
    ///   E_q = (gammaQual/2) ∫_Ω max(0, Q_rel − qStar)² dX.
    /// ρ'(Q)=max(0,Q−qStar) is the redistribution force; ρ''(Q)=1 the
    /// Gauss–Newton Hessian (added to the metric). The energy is FLAT
    /// for well-shaped cells (Q ≤ qStar) — zero force, so it does not
    /// fight the fit — and pushes the degrading tail back to qStar.
    /// Its fixed point is the fit-optimal configuration subject to "no
    /// element worse than qStar", reached without touching the cells
    /// already good; displaced strain flows to under-threshold
    /// neighbours (which offer no resistance), so strain redistributes
    /// implicitly. The barrier (metric + line search) remains the hard
    /// validity wall underneath. gammaQual ≤ 0 disables.
    ///
    /// Default preset (on): gammaQual=1 is an O(1) trade-off weight
    /// against the data term; qStar=1.75 is a dimensionless distortion
    /// threshold (Q_rel=1 is ideal, ≪ qMax=10 the hard wall) and is
    /// h-independent — refining the mesh does not rescale it.
    Real gammaQual = 1;
    Real qStar = Real(1.75);
    /// Size hinge (energy, on the RHS). Companion to the shape hinge:
    /// penalizes cells whose Jacobian drops below a *good* level jStar
    /// (chosen ≫ jSafe):
    ///   E_s = (gammaSize/2) ∫_Ω max(0, jStar − j)² dX.
    /// ρ'(j) = −max(0, jStar − j) (restoring force pushing j up toward
    /// jStar); ρ''(j) = 1 for j < jStar (GN Hessian → metric). Flat for
    /// j ≥ jStar (no fit-fighting). Unlike the shape hinge on Q_rel —
    /// which is scale-invariant and undefined for j ≤ 0 — the size
    /// hinge depends only on j and the cofactor action a_j, so it is
    /// evaluated for ALL cells, INCLUDING inverted ones (j ≤ 0), giving
    /// a strong restoring force that pulls inverted elements back to
    /// validity. The metric barrier + line search remain the hard wall.
    /// gammaSize ≤ 0 disables.
    ///
    /// Default preset (on): gammaSize=1 (O(1) trade-off weight);
    /// jStar=0.3 is a dimensionless volume-ratio threshold (j=1 at
    /// identity, ≫ jSafe=1e-2) and is h-independent.
    Real gammaSize = 1;
    Real jStar = Real(0.3);
    Real omegaMin = 0.1;        ///< active-set threshold on ω.
    Real alphaMin = 1e-4;       ///< line-search floor.
    bool energyLineSearch = true;
    Real jMinRatio = 1e-8;      ///< hard inadmissibility floor.
    Real jLineSearchRatio = 1e-2;
    Real activeRMSTol = 0;      ///< ≤0 ⇒ 4h².
    Real activeSupTol = 0;      ///< ≤0 ⇒ 10h².
    Real energyStagTol = 1e-4;
    Real stepTol = 0;           ///< ≤0 ⇒ 1e-4·h.
    Real pointLocationTolerance = 0; ///< ≤0 ⇒ 1e-10 for moved FE evaluation.
    Real cgRelativeTolerance = 1e-6; ///< relative residual tolerance for CG.
    std::size_t cgMaxIterations = 0; ///< 0 ⇒ min(2000, max(100, 2*ndofs)).
    std::size_t andersonMemory = 3;  ///< 0 disables safeguarded Anderson.
    std::size_t andersonStart = 2;   ///< first iteration where AA may be used.
    Real andersonDamping = 1;        ///< first AA damping trial.
    Real andersonMinDamping = 0.125; ///< smallest AA damping trial.
    std::size_t maxIterations = 200;
    std::size_t quadratureOrder = 0; ///< 0 ⇒ 2·(FE order).
    bool hasInterfaceAttribute = false;
    Geometry::Attribute interfaceAttribute = 0;
    bool trace = false;
    /// If true, also add the nonlinear barrier first variation to the
    /// RHS. Default false: admissibility enters as metric stiffness and
    /// true-geometry line-search validation, not as a competing energy.
    bool includeAdmissibilityGradient = false;
    /// Optimal 1-D rescale of the lifted step along itself:
    ///   β = ⟨d, v⟩_Γ / ⟨v, v⟩_Γ  (surface inner products),
    /// clamped to [1, betaMax]; line search starts at β·v instead of
    /// v. The H¹ lift systematically under-scales the skeleton trace
    /// (gain ≈ surface-weight / M-diagonal ≈ 1/20 at default γ), so
    /// without β the iteration is linearly convergent with ρ ≈ 0.95.
    /// β recovers Newton-matched magnitude while preserving the lift's
    /// smooth admissibility-aware shape. Since β only scales the same
    /// descent direction, the nonlinear line search remains the final
    /// admissibility and energy-decrease guard.
    Real betaMax = 50;
  };

  struct WNGIRReport
  {
    std::size_t iterations = 0;
    Real sigma = 0;
    Real lastAlpha = 0;
    Real acceptedStep = 0;
    Real minJ = 1;
    Real maxQRel = 1;
    Real activeRMS = 0;
    Real activeSup = 0;
    Real activeFraction = 0;
    Real energy = 0;
    const char* exitReason = "iter-budget";
    // Wall-clock breakdown (seconds, accumulated over iterations).
    Real tAssembly = 0;   ///< WNGIR variational problem assembly.
    Real tFactor = 0;     ///< CG setup/preconditioner.
    Real tSolve = 0;      ///< CG iterations.
    Real tLineSearch = 0; ///< true-geometry admissibility + energy LS.
    std::size_t linearIterations = 0;
    Real linearError = 0;
    std::size_t andersonTried = 0;
    std::size_t andersonAccepted = 0;
    Real lastAndersonTheta = 0;
  };

  namespace Detail
  {
    inline Real referenceCellMargin(
        Geometry::Polytope::Type geometry,
        const Math::SpatialPoint& rc,
        std::size_t& mostViolatedFace)
    {
      const Geometry::Polytope::Traits traits(geometry);
      const auto& hs = traits.getHalfSpace();
      Real margin = std::numeric_limits<Real>::infinity();
      mostViolatedFace = 0;
      for (std::size_t j = 0; j < static_cast<std::size_t>(hs.vector.size()); ++j)
      {
        const Real phi = hs.vector[j] - rc.dot(hs.matrix.row(j).transpose());
        if (phi < margin)
        {
          margin = phi;
          mostViolatedFace = j;
        }
      }
      return margin;
    }

    inline Geometry::Point makeTranslatedPoint(
        const Geometry::Point& source,
        const Math::SpatialVector<Real>& yPhysical,
        Real tolerance)
    {
      const auto& sourcePolytope = source.getPolytope();
      const auto& mesh = sourcePolytope.getMesh();
      const std::size_t cd = mesh.getDimension();
      const auto& conn = mesh.getConnectivity();
      const Real tol = tolerance > Real(0)
        ? tolerance
        : Real(64) * std::numeric_limits<Real>::epsilon();

      Index cell = Index(-1);
      Math::SpatialPoint rc;

      auto cellMargin =
        [&](Index candidate, Math::SpatialPoint& out) -> Real
        {
          mesh.getPolytopeTransformation(cd, candidate).inverse(
              out, yPhysical);
          std::size_t face = 0;
          return referenceCellMargin(mesh.getGeometry(cd, candidate), out, face);
        };

      if (sourcePolytope.getDimension() == cd)
      {
        cell = sourcePolytope.getIndex();
        cellMargin(cell, rc);
      }
      else if (sourcePolytope.getDimension() + 1 == cd)
      {
        const Index face = sourcePolytope.getIndex();
        const auto& adjacent = conn.getIncidence(cd - 1, cd).at(face);
        if (adjacent.empty())
          throw std::runtime_error(
              "WNGIR translated-point evaluation failed: face has no adjacent cell.");

        Real bestMargin = -std::numeric_limits<Real>::infinity();
        Math::SpatialPoint bestRc;
        for (const Index candidate : adjacent)
        {
          Math::SpatialPoint candidateRc;
          const Real margin = cellMargin(candidate, candidateRc);
          if (margin > bestMargin)
          {
            bestMargin = margin;
            bestRc = candidateRc;
            cell = candidate;
          }
        }
        rc = bestRc;
      }
      else
      {
        throw std::runtime_error(
            "WNGIR translated-point evaluation requires a cell or face source point.");
      }

      for (std::size_t hop = 0; hop < 64; ++hop)
      {
        mesh.getPolytopeTransformation(cd, cell).inverse(rc, yPhysical);

        std::size_t mostViolatedFace = 0;
        const Real margin = referenceCellMargin(
            mesh.getGeometry(cd, cell), rc, mostViolatedFace);
        if (margin >= -tol)
        {
          const auto it = mesh.getPolytope(cd, cell);
          return Geometry::Point(*it, rc, yPhysical);
        }

        const auto& faces = conn.getIncidence(cd, cd - 1).at(cell);
        if (mostViolatedFace >= faces.size())
          break;

        const Index face = faces[mostViolatedFace];
        if (mesh.isBoundary(face))
          break;

        const auto& adjacent = conn.getIncidence(cd - 1, cd).at(face);
        if (adjacent.size() != 2)
          break;

        const Index next = (adjacent[0] == cell) ? adjacent[1] : adjacent[0];
        if (next == cell)
          break;
        cell = next;
      }

      const auto it = mesh.getPolytope(cd, cell);
      mesh.getPolytopeTransformation(cd, cell).inverse(rc, yPhysical);
      return Geometry::Point(*it, rc, yPhysical);
    }

    template <class Function, class Vector>
    decltype(auto) evaluateTranslatedPoint(
        const Function& f,
        const Geometry::Point& source,
        const Vector& displacement,
        Real tolerance)
    {
      Math::SpatialVector<Real> y(source.getPolytope().getMesh().getDimension());
      const auto& x = source.getPhysicalCoordinates();
      for (std::size_t r = 0; r < static_cast<std::size_t>(y.size()); ++r)
        y(static_cast<Eigen::Index>(r)) =
          x(static_cast<Eigen::Index>(r))
          + displacement(static_cast<Eigen::Index>(r));
      const auto p = makeTranslatedPoint(source, y, tolerance);
      if constexpr (requires { f.getValue(p); })
        return f.getValue(p);
      else
        return f(p);
    }

    inline Math::SpatialMatrix<Real> makeZeroMatrix(std::size_t dim)
    {
      Math::SpatialMatrix<Real> out(dim, dim);
      out.setZero();
      return out;
    }

    template <class FE, class ReferencePoint, class JacobianInverse>
    Math::SpatialMatrix<Real> physicalJacobian(
        const FE& fe,
        std::size_t local,
        const ReferencePoint& rc,
        const JacobianInverse& Jinv,
        std::size_t dim)
    {
      const auto jref = fe.getBasis(local).getJacobian()(rc);
      auto jp = makeZeroMatrix(dim);
      for (std::size_t r = 0; r < dim; ++r)
        for (std::size_t c = 0; c < dim; ++c)
          for (std::size_t a = 0; a < dim; ++a)
            jp(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c))
              += jref(r, a) * Jinv(a, c);
      return jp;
    }

    template <class Displacement>
    Math::SpatialMatrix<Real> deformationGradient(
        const Displacement& u,
        const Geometry::Polytope& polytope,
        const Variational::IntegrationPoint& ip,
        std::size_t dim)
    {
      const auto& fes = u.getFiniteElementSpace();
      const auto& fe = fes.getFiniteElement(
          polytope.getDimension(), polytope.getIndex());
      const auto& dofs = fes.getDOFs(
          polytope.getDimension(), polytope.getIndex());
      const auto Jinv = ip.getPoint().getJacobianInverse();
      const auto& rc = ip.getPoint().getReferenceCoordinates();

      auto F = Math::SpatialMatrix<Real>::Identity(dim, dim);
      for (std::size_t l = 0; l < fe.getCount(); ++l)
        F += u.getData()(dofs[l])
          * physicalJacobian(fe, l, rc, Jinv, dim);
      return F;
    }

    template <class PhiDerived, class GradDerived,
              class TrialFunction, class TestFunction, class Displacement>
    class WNGIRSurfaceObservationMetric final
      : public Variational::LocalBilinearFormIntegratorBase<
          typename TrialFunction::ScalarType>
    {
      public:
        using ScalarType = typename TrialFunction::ScalarType;
        using Parent =
          Variational::LocalBilinearFormIntegratorBase<ScalarType>;
        using PhiType = Variational::RealFunctionBase<PhiDerived>;
        using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

        WNGIRSurfaceObservationMetric(
            const PhiType& phi,
            const GradType& grad,
            const TrialFunction& du,
            const TestFunction& v,
            const Displacement& current,
            const WNGIRParameters& parameters,
            Real sigma2)
          : Parent(du.getLeaf(), v.getLeaf()),
            m_phi(phi.copy()),
            m_grad(grad.copy()),
            m_du(du),
            m_v(v),
            m_current(current),
            m_parameters(parameters),
            m_sigma2(sigma2)
        {}

        WNGIRSurfaceObservationMetric(
            const WNGIRSurfaceObservationMetric& other)
          : Parent(other),
            m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
            m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
            m_du(other.m_du),
            m_v(other.m_v),
            m_current(other.m_current),
            m_parameters(other.m_parameters),
            m_sigma2(other.m_sigma2),
            m_polytope(other.m_polytope)
        {}

        const Geometry::Polytope& getPolytope() const final override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        WNGIRSurfaceObservationMetric& setPolytope(
            const Geometry::Polytope& polytope) final override
        {
          m_polytope = &polytope;

          const std::size_t faceDim = polytope.getDimension();
          const Index faceIdx = polytope.getIndex();
          const auto& trialFES = m_du.getFiniteElementSpace();
          const auto& testFES = m_v.getFiniteElementSpace();
          const auto& trialFE = trialFES.getFiniteElement(faceDim, faceIdx);
          const auto& testFE = testFES.getFiniteElement(faceDim, faceIdx);
          const auto& params = m_parameters.get();
          const std::size_t qOrder = params.quadratureOrder > 0
            ? params.quadratureOrder
            : std::max<std::size_t>(
                2, 2 * std::max(trialFE.getOrder(), testFE.getOrder()));
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(qOrder, polytope.getGeometry());
          const auto& quad = polytope.getQuadrature(qf);

          const std::size_t ntr = trialFE.getCount();
          const std::size_t nte = testFE.getCount();
          m_matrix.resize(
              static_cast<Eigen::Index>(nte),
              static_cast<Eigen::Index>(ntr));
          m_matrix.setZero();

          constexpr Real epsG = Real(1e-12);
          for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
          {
            const auto& pt = quad.getPoint(qp);
            const Variational::IntegrationPoint ip(pt, &qf, qp);
            const Real w = qf.getWeight(qp) * pt.getDistortion();
            const auto uq = m_current.get().getValue(ip);
            Math::SpatialVector<Real> displacement(polytope.getMesh().getDimension());
            displacement.setZero();
            for (std::size_t r = 0; r < static_cast<std::size_t>(displacement.size()); ++r)
              displacement(static_cast<Eigen::Index>(r)) = uq(r);

            const auto moved =
              makeTranslatedPoint(pt, pt.getPhysicalCoordinates() + displacement,
                                  params.pointLocationTolerance);
            const Real r = m_phi->getValue(moved);
            const auto g = m_grad->getValue(moved);
            const Real obsWeight =
                g.dot(g) + epsG
              + (params.residualStabilizedObservationMetric
                  ? (r * r) / m_sigma2 : Real(0));

            const auto& rc = pt.getReferenceCoordinates();
            for (std::size_t te = 0; te < nte; ++te)
            {
              const auto testValue = testFE.getBasis(te)(rc);
              for (std::size_t tr = 0; tr < ntr; ++tr)
              {
                const auto trialValue = trialFE.getBasis(tr)(rc);
                m_matrix(
                    static_cast<Eigen::Index>(te),
                    static_cast<Eigen::Index>(tr))
                  += w * params.gammaObs * obsWeight
                   * trialValue.dot(testValue);
              }
            }
          }
          return *this;
        }

        ScalarType integrate(std::size_t tr, std::size_t te) final override
        {
          return m_matrix(
              static_cast<Eigen::Index>(te),
              static_cast<Eigen::Index>(tr));
        }

        Geometry::Region getRegion() const final override
        {
          return Geometry::Region::Faces;
        }

        WNGIRSurfaceObservationMetric* copy() const noexcept final override
        {
          return new WNGIRSurfaceObservationMetric(*this);
        }

      private:
        std::unique_ptr<PhiType> m_phi;
        std::unique_ptr<GradType> m_grad;
        TrialFunction m_du;
        TestFunction m_v;
        std::reference_wrapper<const Displacement> m_current;
        std::reference_wrapper<const WNGIRParameters> m_parameters;
        Real m_sigma2;
        const Geometry::Polytope* m_polytope = nullptr;
        Math::Matrix<Real> m_matrix;
    };

    template <class PhiDerived, class GradDerived,
              class TestFunction, class Displacement>
    class WNGIRSurfaceForce final
      : public Variational::LinearFormIntegratorBase<
          typename TestFunction::ScalarType>
    {
      public:
        using ScalarType = typename TestFunction::ScalarType;
        using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
        using PhiType = Variational::RealFunctionBase<PhiDerived>;
        using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

        WNGIRSurfaceForce(
            const PhiType& phi,
            const GradType& grad,
            const TestFunction& v,
            const Displacement& current,
            const WNGIRParameters& parameters,
            Real sigma2)
          : Parent(v.getLeaf()),
            m_phi(phi.copy()),
            m_grad(grad.copy()),
            m_v(v),
            m_current(current),
            m_parameters(parameters),
            m_sigma2(sigma2)
        {}

        WNGIRSurfaceForce(const WNGIRSurfaceForce& other)
          : Parent(other),
            m_phi(other.m_phi ? other.m_phi->copy() : nullptr),
            m_grad(other.m_grad ? other.m_grad->copy() : nullptr),
            m_v(other.m_v),
            m_current(other.m_current),
            m_parameters(other.m_parameters),
            m_sigma2(other.m_sigma2),
            m_polytope(other.m_polytope)
        {}

        const Geometry::Polytope& getPolytope() const final override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        WNGIRSurfaceForce& setPolytope(
            const Geometry::Polytope& polytope) final override
        {
          m_polytope = &polytope;

          const std::size_t faceDim = polytope.getDimension();
          const Index faceIdx = polytope.getIndex();
          const auto& testFES = m_v.getFiniteElementSpace();
          const auto& testFE = testFES.getFiniteElement(faceDim, faceIdx);
          const auto& params = m_parameters.get();
          const std::size_t qOrder = params.quadratureOrder > 0
            ? params.quadratureOrder
            : std::max<std::size_t>(2, 2 * testFE.getOrder());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(qOrder, polytope.getGeometry());
          const auto& quad = polytope.getQuadrature(qf);

          const std::size_t nte = testFE.getCount();
          m_vector.resize(static_cast<Eigen::Index>(nte));
          m_vector.setZero();

          for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
          {
            const auto& pt = quad.getPoint(qp);
            const Variational::IntegrationPoint ip(pt, &qf, qp);
            const Real w = qf.getWeight(qp) * pt.getDistortion();
            const auto uq = m_current.get().getValue(ip);
            Math::SpatialVector<Real> displacement(polytope.getMesh().getDimension());
            displacement.setZero();
            for (std::size_t r = 0; r < static_cast<std::size_t>(displacement.size()); ++r)
              displacement(static_cast<Eigen::Index>(r)) = uq(r);

            const auto moved =
              makeTranslatedPoint(pt, pt.getPhysicalCoordinates() + displacement,
                                  params.pointLocationTolerance);
            const Real r = m_phi->getValue(moved);
            const auto g = m_grad->getValue(moved);
            const Real omega = std::exp(-r * r / m_sigma2);
            const auto dGamma = (-omega * r) * g;

            const auto& rc = pt.getReferenceCoordinates();
            for (std::size_t te = 0; te < nte; ++te)
            {
              const auto testValue = testFE.getBasis(te)(rc);
              m_vector(static_cast<Eigen::Index>(te))
                += w * dGamma.dot(testValue);
            }
          }
          return *this;
        }

        ScalarType integrate(std::size_t local) final override
        {
          return m_vector(static_cast<Eigen::Index>(local));
        }

        Geometry::Region getRegion() const final override
        {
          return Geometry::Region::Faces;
        }

        WNGIRSurfaceForce* copy() const noexcept final override
        {
          return new WNGIRSurfaceForce(*this);
        }

      private:
        std::unique_ptr<PhiType> m_phi;
        std::unique_ptr<GradType> m_grad;
        TestFunction m_v;
        std::reference_wrapper<const Displacement> m_current;
        std::reference_wrapper<const WNGIRParameters> m_parameters;
        Real m_sigma2;
        const Geometry::Polytope* m_polytope = nullptr;
        Math::Vector<Real> m_vector;
    };

    template <class TrialFunction, class TestFunction, class Displacement>
    class WNGIRAdmissibilityMetric final
      : public Variational::LocalBilinearFormIntegratorBase<
          typename TrialFunction::ScalarType>
    {
      public:
        using ScalarType = typename TrialFunction::ScalarType;
        using Parent =
          Variational::LocalBilinearFormIntegratorBase<ScalarType>;

        WNGIRAdmissibilityMetric(
            const TrialFunction& du,
            const TestFunction& v,
            const Displacement& current,
            const WNGIRParameters& parameters)
          : Parent(du.getLeaf(), v.getLeaf()),
            m_du(du),
            m_v(v),
            m_current(current),
            m_parameters(parameters)
        {}

        WNGIRAdmissibilityMetric(
            const WNGIRAdmissibilityMetric&) = default;

        const Geometry::Polytope& getPolytope() const final override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        WNGIRAdmissibilityMetric& setPolytope(
            const Geometry::Polytope& polytope) final override
        {
          m_polytope = &polytope;

          const std::size_t dim = polytope.getDimension();
          const Real d = static_cast<Real>(dim);
          const auto geometry = polytope.getGeometry();
          const auto idx = polytope.getIndex();

          const auto& trialFES = m_du.getFiniteElementSpace();
          const auto& testFES = m_v.getFiniteElementSpace();
          const auto& trialFE = trialFES.getFiniteElement(dim, idx);
          const auto& testFE = testFES.getFiniteElement(dim, idx);
          const auto& params = m_parameters.get();
          const std::size_t qOrder = params.quadratureOrder > 0
            ? params.quadratureOrder
            : std::max<std::size_t>(2, 2 * trialFE.getOrder());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(qOrder, geometry);
          const auto& quad = polytope.getQuadrature(qf);

          const std::size_t ntr = trialFE.getCount();
          const std::size_t nte = testFE.getCount();
          m_matrix.resize(
              static_cast<Eigen::Index>(nte),
              static_cast<Eigen::Index>(ntr));
          m_matrix.setZero();

          for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
          {
            const auto& pt = quad.getPoint(qp);
            const Variational::IntegrationPoint ip(pt, &qf, qp);
            const Real w = qf.getWeight(qp) * pt.getDistortion();
            const auto F =
              deformationGradient(m_current.get(), polytope, ip, dim);
            const Real jK = F.determinant();
            // Near-singular: cofactor via inverse is unreliable; skip.
            if (std::abs(jK) < Real(1e-14))
              continue;
            // Shape terms (Q_rel and its barriers/hinge) need j > 0; the
            // size hinge needs only j and a_j, so it runs for ALL j —
            // including inverted cells, which it pulls back to validity.
            const bool shapeOK = jK > Real(0);
            Real frob2 = 0;
            Real qK = 0;
            Real sJ0 = 0, sQ0 = 0;
            bool jActive = false, qActive = false, qualActive = false;
            if (shapeOK)
            {
              frob2 = F.squaredNorm();
              qK = frob2 / (d * std::pow(jK, Real(2) / d));
              sJ0 = jK - params.jSafe;
              sQ0 = params.qMax - qK;
              jActive = params.gammaJ > Real(0)
                && sJ0 > Real(0) && sJ0 < params.s0J;
              qActive = params.gammaQ > Real(0)
                && sQ0 > Real(0) && sQ0 < params.s0Q;
              qualActive = params.gammaQual > Real(0) && qK > params.qStar;
            }
            const bool sizeActive =
              params.gammaSize > Real(0) && jK < params.jStar;
            if (!jActive && !qActive && !qualActive && !sizeActive)
              continue;

            const auto Jinv = pt.getJacobianInverse();
            const auto& rc = pt.getReferenceCoordinates();
            const auto FinvT = F.inverse().transpose();
            const auto dQdF = shapeOK
              ? Math::SpatialMatrix<Real>(
                  (Real(2) / d) * std::pow(jK, -Real(2) / d)
                  * (F - (frob2 / d) * FinvT))
              : makeZeroMatrix(dim);

            std::vector<Real> aJTrial(ntr), aQTrial(ntr);
            std::vector<Real> aJTest(nte), aQTest(nte);
            auto fillActions =
              [&](const auto& fe, std::vector<Real>& aJ, std::vector<Real>& aQ)
            {
              for (std::size_t l = 0; l < fe.getCount(); ++l)
              {
                const auto jp = physicalJacobian(fe, l, rc, Jinv, dim);
                Real accJ = 0;
                Real accQ = 0;
                for (std::size_t r = 0; r < dim; ++r)
                  for (std::size_t c = 0; c < dim; ++c)
                  {
                    const auto rr = static_cast<Eigen::Index>(r);
                    const auto cc = static_cast<Eigen::Index>(c);
                    accJ += jK * FinvT(rr, cc) * jp(rr, cc);
                    accQ += dQdF(rr, cc) * jp(rr, cc);
                  }
                aJ[l] = accJ;
                aQ[l] = accQ;
              }
            };
            fillActions(trialFE, aJTrial, aQTrial);
            fillActions(testFE, aJTest, aQTest);

            if (jActive)
            {
              const Real bpp = Real(1) / (sJ0 * sJ0);
              for (std::size_t te = 0; te < nte; ++te)
                for (std::size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(
                      static_cast<Eigen::Index>(te),
                      static_cast<Eigen::Index>(tr))
                    += w * params.gammaJ * bpp
                    * aJTrial[tr] * aJTest[te];
            }
            if (qActive)
            {
              const Real bpp = Real(1) / (sQ0 * sQ0);
              for (std::size_t te = 0; te < nte; ++te)
                for (std::size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(
                      static_cast<Eigen::Index>(te),
                      static_cast<Eigen::Index>(tr))
                    += w * params.gammaQ * bpp
                    * aQTrial[tr] * aQTest[te];
            }
            if (qualActive)
            {
              // Gauss–Newton Hessian of the quality hinge: ρ''(Q)=1.
              for (std::size_t te = 0; te < nte; ++te)
                for (std::size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(
                      static_cast<Eigen::Index>(te),
                      static_cast<Eigen::Index>(tr))
                    += w * params.gammaQual
                    * aQTrial[tr] * aQTest[te];
            }
            if (sizeActive)
            {
              // GN Hessian of the size hinge: ρ''(j)=1 for j < jStar.
              for (std::size_t te = 0; te < nte; ++te)
                for (std::size_t tr = 0; tr < ntr; ++tr)
                  m_matrix(
                      static_cast<Eigen::Index>(te),
                      static_cast<Eigen::Index>(tr))
                    += w * params.gammaSize
                    * aJTrial[tr] * aJTest[te];
            }
          }
          return *this;
        }

        ScalarType integrate(std::size_t tr, std::size_t te) final override
        {
          return m_matrix(
              static_cast<Eigen::Index>(te),
              static_cast<Eigen::Index>(tr));
        }

        Geometry::Region getRegion() const final override
        {
          return Geometry::Region::Cells;
        }

        WNGIRAdmissibilityMetric* copy() const noexcept final override
        {
          return new WNGIRAdmissibilityMetric(*this);
        }

      private:
        TrialFunction m_du;
        TestFunction m_v;
        std::reference_wrapper<const Displacement> m_current;
        std::reference_wrapper<const WNGIRParameters> m_parameters;
        const Geometry::Polytope* m_polytope = nullptr;
        Math::Matrix<Real> m_matrix;
    };

    template <class TestFunction, class Displacement>
    class WNGIRAdmissibilityGradient final
      : public Variational::LinearFormIntegratorBase<
          typename TestFunction::ScalarType>
    {
      public:
        using ScalarType = typename TestFunction::ScalarType;
        using Parent = Variational::LinearFormIntegratorBase<ScalarType>;

        WNGIRAdmissibilityGradient(
            const TestFunction& v,
            const Displacement& current,
            const WNGIRParameters& parameters)
          : Parent(v.getLeaf()),
            m_v(v),
            m_current(current),
            m_parameters(parameters)
        {}

        WNGIRAdmissibilityGradient(
            const WNGIRAdmissibilityGradient&) = default;

        const Geometry::Polytope& getPolytope() const final override
        {
          assert(m_polytope);
          return *m_polytope;
        }

        WNGIRAdmissibilityGradient& setPolytope(
            const Geometry::Polytope& polytope) final override
        {
          m_polytope = &polytope;

          const std::size_t dim = polytope.getDimension();
          const Real d = static_cast<Real>(dim);
          const auto geometry = polytope.getGeometry();
          const auto idx = polytope.getIndex();

          const auto& testFES = m_v.getFiniteElementSpace();
          const auto& testFE = testFES.getFiniteElement(dim, idx);
          const std::size_t qOrder = m_parameters.get().quadratureOrder > 0
            ? m_parameters.get().quadratureOrder
            : std::max<std::size_t>(2, 2 * testFE.getOrder());
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(qOrder, geometry);
          const auto& quad = polytope.getQuadrature(qf);

          const std::size_t nte = testFE.getCount();
          m_vector.resize(static_cast<Eigen::Index>(nte));
          m_vector.setZero();

          for (std::size_t qp = 0; qp < quad.getSize(); ++qp)
          {
            const auto& pt = quad.getPoint(qp);
            const Variational::IntegrationPoint ip(pt, &qf, qp);
            const Real w = qf.getWeight(qp) * pt.getDistortion();
            const auto F =
              deformationGradient(m_current.get(), polytope, ip, dim);
            const Real jK = F.determinant();
            if (std::abs(jK) < Real(1e-14))
              continue;
            const auto& params = m_parameters.get();
            // Shape force needs j > 0; size force runs for ALL j.
            const bool shapeOK = jK > Real(0);
            Real frob2 = 0, qK = 0, sJ0 = 0, sQ0 = 0;
            bool jActive = false, qActive = false, qualActive = false;
            if (shapeOK)
            {
              frob2 = F.squaredNorm();
              qK = frob2 / (d * std::pow(jK, Real(2) / d));
              sJ0 = jK - params.jSafe;
              sQ0 = params.qMax - qK;
              jActive = params.includeAdmissibilityGradient
                && params.gammaJ > Real(0)
                && sJ0 > Real(0) && sJ0 < params.s0J;
              qActive = params.includeAdmissibilityGradient
                && params.gammaQ > Real(0)
                && sQ0 > Real(0) && sQ0 < params.s0Q;
              qualActive = params.gammaQual > Real(0) && qK > params.qStar;
            }
            const bool sizeActive =
              params.gammaSize > Real(0) && jK < params.jStar;
            if (!jActive && !qActive && !qualActive && !sizeActive)
              continue;

            const auto Jinv = pt.getJacobianInverse();
            const auto& rc = pt.getReferenceCoordinates();
            const auto FinvT = F.inverse().transpose();
            const auto dQdF = shapeOK
              ? Math::SpatialMatrix<Real>(
                  (Real(2) / d) * std::pow(jK, -Real(2) / d)
                  * (F - (frob2 / d) * FinvT))
              : makeZeroMatrix(dim);

            for (std::size_t te = 0; te < nte; ++te)
            {
              const auto jp = physicalJacobian(testFE, te, rc, Jinv, dim);
              Real aJ = 0;
              Real aQ = 0;
              for (std::size_t r = 0; r < dim; ++r)
                for (std::size_t c = 0; c < dim; ++c)
                {
                  const auto rr = static_cast<Eigen::Index>(r);
                  const auto cc = static_cast<Eigen::Index>(c);
                  aJ += jK * FinvT(rr, cc) * jp(rr, cc);
                  aQ += dQdF(rr, cc) * jp(rr, cc);
                }

              Real val = 0;
              if (jActive)
              {
                const Real bp = -Real(1) / sJ0 + Real(1) / params.s0J;
                val += -params.gammaJ * bp * aJ;
              }
              if (qActive)
              {
                const Real bp = -Real(1) / sQ0 + Real(1) / params.s0Q;
                val += params.gammaQ * bp * aQ;
              }
              if (qualActive)
              {
                // −DE_q density: ρ'(Q)=max(0,Q−qStar), dQ-direction aQ.
                val += -params.gammaQual * (qK - params.qStar) * aQ;
              }
              if (sizeActive)
              {
                // −DE_s density: +γ_s·max(0,jStar−j)·a_j pushes j up.
                val += params.gammaSize * (params.jStar - jK) * aJ;
              }
              m_vector(static_cast<Eigen::Index>(te)) += w * val;
            }
          }
          return *this;
        }

        ScalarType integrate(std::size_t local) final override
        {
          return m_vector(static_cast<Eigen::Index>(local));
        }

        Geometry::Region getRegion() const final override
        {
          return Geometry::Region::Cells;
        }

        WNGIRAdmissibilityGradient* copy() const noexcept final override
        {
          return new WNGIRAdmissibilityGradient(*this);
        }

      private:
        TestFunction m_v;
        std::reference_wrapper<const Displacement> m_current;
        std::reference_wrapper<const WNGIRParameters> m_parameters;
        const Geometry::Polytope* m_polytope = nullptr;
        Math::Vector<Real> m_vector;
    };
  }

  /// Reusable-scratch overload. The trial/test functions and the two
  /// `Variational::Problem` objects (step + bulk) are owned by the
  /// caller (`class WNGIR`) and threaded in by reference so they are
  /// constructed once and reused across iterations and frames. Only
  /// their bodies are reassigned/reassembled per iteration.
  template <class Mesh, class Displacement, class PhiDerived, class GradDerived,
            class TrialT, class TestT, class ProblemT>
  WNGIRReport solveWNGIR(
      const Mesh& mesh,
      Displacement& u,
      const std::vector<Rodin::Index>& interfaceFacets,
      const Variational::RealFunctionBase<PhiDerived>& phi,
      const Variational::VectorFunctionBase<Real, GradDerived>& grad,
      const WNGIRParameters& p,
      TrialT& duStep,
      TestT& vStep,
      TrialT& duBulk,
      TestT& vBulk,
      ProblemT& stepProblem,
      ProblemT& bulkProblem)
  {
    using Vec = Math::SpatialVector<Real>;
    using Mat = Math::SpatialMatrix<Real>;
    using Rodin::Index;
    auto zeroVec = [](std::size_t dim) {
      Vec out(dim);
      out.setZero();
      return out;
    };
    auto zeroMat = [](std::size_t dim) {
      Mat out(dim, dim);
      out.setZero();
      return out;
    };

    WNGIRReport rep;
    const Real h = p.h;
    const Real gammaM = p.gammaM > Real(0) ? p.gammaM : Real(1) / h;
    const Real gammaH = p.gammaH > Real(0) ? p.gammaH : Real(1) / h;
    const Real ellM = p.ellM > Real(0) ? p.ellM : Real(3) * h;
    const Real activeRMSTol =
      p.activeRMSTol > Real(0) ? p.activeRMSTol : Real(4) * h * h;
    const Real activeSupTol =
      p.activeSupTol > Real(0) ? p.activeSupTol : Real(10) * h * h;
    const Real stepTol =
      p.stepTol > Real(0) ? p.stepTol : Real(1e-4) * h;
    constexpr Real epsG = Real(1e-12);

    const auto& fes = u.getFiniteElementSpace();
    const std::size_t meshDim = mesh.getDimension();
    const Real d = static_cast<Real>(meshDim);

    // =================================================================
    // Per-frame tabulation (geometry is fixed; only u changes).
    // =================================================================

    // ---- Facet tables: trace basis at facet quadrature points ----
    struct FacetQP
    {
      Real w = 0;                 ///< quadrature weight × |J_facet|.
      Vec X;                      ///< physical coordinates.
      Math::SpatialPoint rc;      ///< reference coordinates on the face.
      std::vector<Vec> val;       ///< trace basis values (per local).
    };
    struct FacetTable
    {
      Index facet = 0;
      std::vector<Index> dofs;
      std::vector<FacetQP> qps;
    };
    std::vector<FacetTable> facetTables;
    facetTables.reserve(interfaceFacets.size());
    for (const Index facetIdx : interfaceFacets)
    {
      const auto face = mesh.getFace(facetIdx);
      const auto& fe = fes.getFiniteElement(meshDim - 1, facetIdx);
      const std::size_t nLocal = fe.getCount();
      const std::size_t feOrder = fe.getOrder();
      const std::size_t qOrder = p.quadratureOrder > 0
        ? p.quadratureOrder
        : std::max<std::size_t>(2, 2 * feOrder);
      const auto& qf = QF::PolytopeQuadratureFormula::get(
          qOrder, face->getGeometry());
      const auto& quad = face->getQuadrature(qf);

      FacetTable ft;
      ft.facet = facetIdx;
      ft.dofs.resize(nLocal);
      for (std::size_t l = 0; l < nLocal; ++l)
        ft.dofs[l] = fes.getGlobalIndex({meshDim - 1, facetIdx}, l);
      ft.qps.resize(quad.getSize());
      for (std::size_t q = 0; q < quad.getSize(); ++q)
      {
        const auto& pt = quad.getPoint(q);
        auto& fq = ft.qps[q];
        fq.w = qf.getWeight(q) * pt.getDistortion();
        const auto& pc = pt.getPhysicalCoordinates();
        fq.X = zeroVec(meshDim);
        for (std::size_t r = 0; r < meshDim; ++r)
          fq.X(static_cast<Eigen::Index>(r)) = pc(r);
        const auto& rc = pt.getReferenceCoordinates();
        fq.rc = rc;
        fq.val.resize(nLocal);
        for (std::size_t l = 0; l < nLocal; ++l)
        {
          const auto bv = fe.getBasis(l)(rc);
          fq.val[l] = zeroVec(meshDim);
          for (std::size_t r = 0; r < meshDim; ++r)
            fq.val[l](static_cast<Eigen::Index>(r)) = bv(r);
        }
      }
      facetTables.push_back(std::move(ft));
    }

    auto sourcePoint =
      [&](const FacetTable& ft, const FacetQP& fq)
    {
      const auto face = mesh.getFace(ft.facet);
      return Geometry::Point(*face, fq.rc, fq.X);
    };

    // ---- Fixed Welsch scale: σ = max(3h, q90(|φ| on Γ at u=0)) ----
    std::vector<Real> r0abs;
    for (const auto& ft : facetTables)
      for (const auto& fq : ft.qps)
      {
        Vec zero = zeroVec(meshDim);
        const auto src = sourcePoint(ft, fq);
        r0abs.push_back(std::abs(
              Detail::evaluateTranslatedPoint(
                phi, src, zero, p.pointLocationTolerance)));
      }
    Real sigma = Real(3) * h;
    if (!r0abs.empty())
    {
      const std::size_t k90 = static_cast<std::size_t>(
          Real(0.9) * static_cast<Real>(r0abs.size() - 1));
      std::nth_element(r0abs.begin(), r0abs.begin() + k90, r0abs.end());
      sigma = std::max(sigma, r0abs[k90]);
    }
    const Real sigma2 = sigma * sigma;
    rep.sigma = sigma;
    if (!p.hasInterfaceAttribute)
    {
      rep.exitReason = "missing-interface-attribute";
      return rep;
    }

    // ---- Cell tables: validation geometry for line search/diagnostics.
    struct CellQP
    {
      Real w = 0;                  ///< quadrature weight × |J_cell|.
      std::vector<Mat> jac;        ///< physical basis Jacobians.
    };
    struct CellTable
    {
      std::vector<Index> dofs;
      std::vector<CellQP> qps;
    };
    std::vector<CellTable> cellTables;
    cellTables.reserve(mesh.getCellCount());
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const Index cellIdx = cellIt->getIndex();
      const auto& fe = fes.getFiniteElement(meshDim, cellIdx);
      const std::size_t nLocal = fe.getCount();
      const std::size_t feOrder = fe.getOrder();
      const std::size_t qOrder = p.quadratureOrder > 0
        ? p.quadratureOrder
        : std::max<std::size_t>(2, 2 * feOrder);
      const auto& qf = QF::PolytopeQuadratureFormula::get(
          qOrder, cellIt->getGeometry());
      const auto& quad = cellIt->getQuadrature(qf);

      CellTable ct;
      ct.dofs.resize(nLocal);
      for (std::size_t l = 0; l < nLocal; ++l)
        ct.dofs[l] = fes.getGlobalIndex({meshDim, cellIdx}, l);

      ct.qps.resize(quad.getSize());
      for (std::size_t q = 0; q < quad.getSize(); ++q)
      {
        const auto& pt = quad.getPoint(q);
        auto& cq = ct.qps[q];
        cq.w = qf.getWeight(q) * pt.getDistortion();
        const auto& rc = pt.getReferenceCoordinates();
        // Physical Jacobians: phys = ref · Jinv, with the (affine)
        // cell Jacobian inverse from the point.
        const auto Jinv = pt.getJacobianInverse();
        cq.jac.resize(nLocal);
        for (std::size_t l = 0; l < nLocal; ++l)
        {
          const auto jref = fe.getBasis(l).getJacobian()(rc);
          Mat jp = zeroMat(meshDim);
          for (std::size_t r = 0; r < meshDim; ++r)
            for (std::size_t c = 0; c < meshDim; ++c)
              for (std::size_t a = 0; a < meshDim; ++a)
                jp(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c))
                  += jref(r, a) * Jinv(a, c);
          cq.jac[l] = jp;
        }
      }
      cellTables.push_back(std::move(ct));
    }

    // =================================================================
    // Per-iteration helpers (pure contractions of the tables).
    // =================================================================

    // u (or any coefficient vector) evaluated at a facet QP.
    auto fieldAtFacetQP =
      [&](const Math::Vector<Real>& coef,
          const FacetTable& ft, const FacetQP& fq) -> Vec
    {
      Vec out = zeroVec(meshDim);
      for (std::size_t l = 0; l < ft.dofs.size(); ++l)
        out += coef(ft.dofs[l]) * fq.val[l];
      return out;
    };

    // F = I + ∇u at a cell QP.
    auto deformationAtCellQP =
      [&](const Math::Vector<Real>& coef,
          const CellTable& ct, const CellQP& cq) -> Mat
    {
      Mat F = Mat::Identity(meshDim, meshDim);
      for (std::size_t l = 0; l < ct.dofs.size(); ++l)
      {
        const Real c = coef(ct.dofs[l]);
        F += c * cq.jac[l];
      }
      return F;
    };

    // Closed-form admissibility over all validation points.
    struct FastAdm
    {
      Real minJ = std::numeric_limits<Real>::infinity();
      Real maxQ = 0;
      std::size_t inadmissibleCount = 0;
    };
    auto fastAdmissibility = [&](const Math::Vector<Real>& uData)
    {
      FastAdm a;
      for (const auto& ct : cellTables)
        for (const auto& cq : ct.qps)
        {
          const Mat F = deformationAtCellQP(uData, ct, cq);
          const Real j = F.determinant();
          if (j < a.minJ) a.minJ = j;
          if (j <= p.jMinRatio) ++a.inadmissibleCount;
          if (j > Real(0))
          {
            const Real q =
              F.squaredNorm() / (d * std::pow(j, Real(2) / d));
            if (q > a.maxQ) a.maxQ = q;
          }
        }
      return a;
    };

    // Welsch surface energy + active-set diagnostics on Γ_ψ,h.
    struct SurfaceState
    {
      Real energy = 0;
      Real activeLen = 0, totalLen = 0;
      Real activeRMS = 0, activeSup = 0;
    };
    auto surfaceState = [&](const Math::Vector<Real>& uData)
    {
      SurfaceState s;
      Real sq = 0;
      for (const auto& ft : facetTables)
        for (const auto& fq : ft.qps)
        {
          const Vec y = fq.X + fieldAtFacetQP(uData, ft, fq);
          const auto src = sourcePoint(ft, fq);
          const Vec displacement = y - fq.X;
          const Real r =
            Detail::evaluateTranslatedPoint(
                phi, src, displacement, p.pointLocationTolerance);
          const Real omega = std::exp(-r * r / sigma2);
          s.totalLen += fq.w;
          s.energy += fq.w * Real(0.5) * sigma2 * (Real(1) - omega);
          if (omega >= p.omegaMin)
          {
            s.activeLen += fq.w;
            sq += fq.w * r * r;
            s.activeSup = std::max(s.activeSup, std::abs(r));
          }
        }
      s.activeRMS = s.activeLen > Real(0)
        ? std::sqrt(sq / s.activeLen) : Real(0);
      return s;
    };

    Real ePrev = surfaceState(u.getData()).energy;
    using Clock = std::chrono::steady_clock;
    auto secondsSince = [](Clock::time_point t0) -> Real
    {
      return std::chrono::duration<Real>(Clock::now() - t0).count();
    };

    struct AndersonState
    {
      std::size_t memory = 0;
      std::vector<Math::Vector<Real>> u;
      std::vector<Math::Vector<Real>> g;
      std::vector<Math::Vector<Real>> f;

      void push(
          const Math::Vector<Real>& uk,
          const Math::Vector<Real>& gk)
      {
        u.push_back(uk);
        g.push_back(gk);
        f.push_back(gk - uk);
        while (u.size() > memory + 1)
        {
          u.erase(u.begin());
          g.erase(g.begin());
          f.erase(f.begin());
        }
      }

      bool canAccelerate() const
      {
        return memory > 0 && f.size() >= 2;
      }

      Math::Vector<Real> candidate() const
      {
        const std::size_t p = std::min<std::size_t>(memory, f.size() - 1);
        const std::size_t start = f.size() - p - 1;
        const Eigen::Index n = f.back().size();
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> dF(n, p);
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> dG(n, p);
        for (std::size_t j = 0; j < p; ++j)
        {
          dF.col(static_cast<Eigen::Index>(j)) =
            f[start + j + 1] - f[start + j];
          dG.col(static_cast<Eigen::Index>(j)) =
            g[start + j + 1] - g[start + j];
        }
        const Eigen::Matrix<Real, Eigen::Dynamic, 1> gamma =
          dF.colPivHouseholderQr().solve(f.back());
        return g.back() - dG * gamma;
      }
    };

    AndersonState anderson;
    anderson.memory = p.andersonMemory;

    // =================================================================
    // Constant bulk metric, pre-assembled once per frame.
    // =================================================================
    // M_bulk = γ_M ∫ v·z + γ_H ℓ² ∫ ∇v:∇z is independent of u, so it is
    // assembled and Dirichlet-eliminated once and reused across all
    // iterations. Only the surface observation metric M_obs, the
    // admissibility metric K_adm, and the RHS depend on u and are
    // reassembled per iteration. Both this form and the per-iteration
    // form carry the same homogeneous Dirichlet elimination, so each
    // constrained boundary row is identity in both; the operator sum is
    // ≤ 2·I there, harmless since the boundary RHS is 0 ⇒ u_bd = 0.
    auto zeroForcing = Variational::VectorFunction(
        meshDim,
        [&](const Geometry::Point&)
        {
          return zeroVec(meshDim);
        });
    // duBulk/vBulk/bulkProblem are caller-owned scratch (reused).
    bulkProblem =
        Variational::Integral(gammaM * duBulk, vBulk)
      + Variational::Integral(
          gammaH * ellM * ellM * Variational::Jacobian(duBulk),
          Variational::Jacobian(vBulk))
      + Variational::Integral(zeroForcing, vBulk)
      + Variational::DirichletBC(duBulk, zeroForcing);
    bulkProblem.assemble();
    const Math::SparseMatrix<Real> bulkOperator =
        bulkProblem.getLinearSystem().getOperator();

    // =================================================================
    // Nonlinear iteration.
    // =================================================================
    for (; rep.iterations < p.maxIterations; ++rep.iterations)
    {
      // ---- One variational assembly: bulk + surface + admissibility ----
      auto tic = Clock::now();
      auto zeroBoundary = Variational::VectorFunction(
          meshDim,
          [&](const Geometry::Point&)
          {
            return zeroVec(meshDim);
          });

      Detail::WNGIRAdmissibilityMetric admMetric(duStep, vStep, u, p);
      Detail::WNGIRAdmissibilityGradient admGradient(vStep, u, p);
      Detail::WNGIRSurfaceObservationMetric obsMetric(
          phi, grad, duStep, vStep, u, p, sigma2);
      obsMetric.over(p.interfaceAttribute);
      Detail::WNGIRSurfaceForce surfaceForce(
          phi, grad, vStep, u, p, sigma2);
      surfaceForce.over(p.interfaceAttribute);

      // Only the u-dependent terms are reassembled here: the surface
      // observation metric, the admissibility + hinge-quality metric,
      // and the Welsch first-variation RHS. The constant bulk metric is
      // added to the operator below from the pre-assembled
      // `bulkOperator`. `stepProblem` is caller-owned scratch; the body
      // is reassigned. The gradient integrator carries the (optional)
      // barrier first variation AND the hinge-quality force, so it is
      // on the RHS whenever either is active.
      const bool useGradientForce =
        p.includeAdmissibilityGradient
        || p.gammaQual > Real(0)
        || p.gammaSize > Real(0);
      if (useGradientForce)
      {
        stepProblem =
            obsMetric
          + admMetric
          - surfaceForce
          - admGradient
          + Variational::DirichletBC(duStep, zeroBoundary);
      }
      else
      {
        stepProblem =
            obsMetric
          + admMetric
          - surfaceForce
          + Variational::DirichletBC(duStep, zeroBoundary);
      }
      stepProblem.assemble();
      rep.tAssembly += secondsSince(tic);

      // ---- Solve assembled WNGIR step problem ----
      // Operator = (M_obs + K_adm + boundary identity) + M_bulk.
      tic = Clock::now();
      const auto& stepSystem = stepProblem.getLinearSystem();
      Math::SparseMatrix<Real> stepOperator = stepSystem.getOperator();
      stepOperator += bulkOperator;
      using CGSolver = Eigen::ConjugateGradient<
        Eigen::SparseMatrix<Real>,
        Eigen::Lower | Eigen::Upper,
        Eigen::DiagonalPreconditioner<Real>>;
      CGSolver cg;
      const Eigen::Index ndofs = stepOperator.rows();
      const std::size_t cgMaxIterations = p.cgMaxIterations > 0
        ? p.cgMaxIterations
        : std::min<std::size_t>(
            2000,
            std::max<std::size_t>(
              100,
              2 * static_cast<std::size_t>(std::max<Eigen::Index>(ndofs, 1))));
      cg.setMaxIterations(static_cast<int>(cgMaxIterations));
      cg.setTolerance(p.cgRelativeTolerance > Real(0)
          ? p.cgRelativeTolerance : Real(1e-6));
      cg.compute(stepOperator);
      if (cg.info() != Eigen::Success)
      {
        rep.exitReason = "solve-cg-setup-failed";
        break;
      }
      rep.tFactor += secondsSince(tic);
      tic = Clock::now();
      Math::Vector<Real> vK = cg.solve(stepSystem.getVector());
      rep.linearIterations += static_cast<std::size_t>(
          std::max<int>(cg.iterations(), 0));
      rep.linearError = cg.error();
      if (cg.info() != Eigen::Success)
      {
        rep.exitReason = "solve-cg-failed";
        break;
      }
      if (!vK.allFinite())
      {
        rep.exitReason = "solve-nonfinite";
        break;
      }
      rep.tSolve += secondsSince(tic);

      // ---- 6. Optimal step rescale β = ⟨d,v⟩_Γ/⟨v,v⟩_Γ ----
      Real beta = Real(1);
      Real rawDemonsRMS = Real(0);
      Real liftedTraceRMS = Real(0);
      Real scaledTraceRMS = Real(0);
      if (p.betaMax > Real(1) || p.trace)
      {
        Real bNum = 0, bDen = 0, dDen = 0, gammaMeasure = 0;
        for (const auto& ft : facetTables)
          for (const auto& fq : ft.qps)
          {
            const Vec y = fq.X + fieldAtFacetQP(u.getData(), ft, fq);
            const auto src = sourcePoint(ft, fq);
            const Vec displacement = y - fq.X;
            const Real r =
              Detail::evaluateTranslatedPoint(
                  phi, src, displacement, p.pointLocationTolerance);
            const Vec g =
              Detail::evaluateTranslatedPoint(
                  grad, src, displacement, p.pointLocationTolerance);
            const Real obsWeight =
                g.dot(g) + epsG
              + (p.residualStabilizedObservationMetric
                  ? (r * r) / sigma2 : Real(0));
            const Real omega = std::exp(-r * r / sigma2);
            const Vec dVec = (-omega * r / obsWeight) * g;
            const Vec v = fieldAtFacetQP(vK, ft, fq);
            bNum += fq.w * dVec.dot(v);
            bDen += fq.w * v.dot(v);
            dDen += fq.w * dVec.dot(dVec);
            gammaMeasure += fq.w;
          }
        if (gammaMeasure > Real(0))
        {
          rawDemonsRMS = std::sqrt(std::max<Real>(Real(0), dDen) / gammaMeasure);
          liftedTraceRMS = std::sqrt(std::max<Real>(Real(0), bDen) / gammaMeasure);
        }
        if (p.betaMax > Real(1) && bDen > Real(0) && std::isfinite(bNum))
        {
          beta = std::clamp(bNum / bDen, Real(1), p.betaMax);
          vK *= beta;
        }
        scaledTraceRMS = beta * liftedTraceRMS;
      }

      const Real maxStep = vK.cwiseAbs().maxCoeff();
      if (maxStep <= stepTol)
      {
        rep.exitReason = "step-below-stepTol";
        break;
      }

      // ---- 7. Nonlinear line search on TRUE geometry ----
      tic = Clock::now();
      const Math::Vector<Real> previousU = u.getData();
      Real alpha = Real(1);
      bool accepted = false;
      std::size_t backtracks = 0;
      FastAdm adm{};
      Real eTrial = std::numeric_limits<Real>::infinity();
      while (alpha >= p.alphaMin)
      {
        const Math::Vector<Real> uTrial = previousU + alpha * vK;
        adm = fastAdmissibility(uTrial);
        const bool jOK =
          adm.inadmissibleCount == 0
          && adm.minJ > p.jLineSearchRatio;
        const bool qOK = adm.maxQ < p.qMax;
        bool eOK = true;
        if (jOK && qOK && p.energyLineSearch)
        {
          eTrial = surfaceState(uTrial).energy;
          eOK = std::isfinite(eTrial) && eTrial <= ePrev;
        }
        if (jOK && qOK && eOK)
        {
          u.getData() = uTrial;
          accepted = true;
          break;
        }
        alpha *= Real(0.5);
        ++backtracks;
      }
      rep.tLineSearch += secondsSince(tic);
      if (!accepted)
      {
        u.getData() = previousU;
        rep.exitReason = "line-search-failure";
        break;
      }

      Math::Vector<Real> acceptedU = u.getData();
      FastAdm acceptedAdm = adm;
      SurfaceState acceptedSurf = surfaceState(acceptedU);
      Real acceptedEnergy =
        p.energyLineSearch && std::isfinite(eTrial)
          ? eTrial : acceptedSurf.energy;

      // ---- 8. Safeguarded Anderson acceleration of the accepted map ----
      bool aaTried = false;
      bool aaAccepted = false;
      Real aaTheta = Real(0);
      anderson.push(previousU, acceptedU);
      if (anderson.canAccelerate()
          && rep.iterations + 1 >= p.andersonStart)
      {
        aaTried = true;
        ++rep.andersonTried;
        const Math::Vector<Real> aaFull = anderson.candidate();
        const Math::Vector<Real> aaDelta = aaFull - acceptedU;
        Real theta = std::clamp(
            p.andersonDamping, Real(0), Real(1));
        const Real thetaMin = std::clamp(
            p.andersonMinDamping, Real(0), theta);
        while (theta >= thetaMin && theta > Real(0))
        {
          const Math::Vector<Real> aaTrial =
            acceptedU + theta * aaDelta;
          const FastAdm aaAdm = fastAdmissibility(aaTrial);
          const bool jOK =
            aaAdm.inadmissibleCount == 0
            && aaAdm.minJ > p.jLineSearchRatio;
          const bool qOK = aaAdm.maxQ < p.qMax;
          if (jOK && qOK)
          {
            const SurfaceState aaSurf = surfaceState(aaTrial);
            const bool eOK =
              !p.energyLineSearch
              || (std::isfinite(aaSurf.energy)
                  && aaSurf.energy <= acceptedEnergy);
            if (eOK)
            {
              acceptedU = aaTrial;
              acceptedAdm = aaAdm;
              acceptedSurf = aaSurf;
              acceptedEnergy = aaSurf.energy;
              aaAccepted = true;
              aaTheta = theta;
              u.getData() = acceptedU;
              ++rep.andersonAccepted;
              break;
            }
          }
          theta *= Real(0.5);
        }
      }

      rep.lastAlpha = alpha;
      rep.lastAndersonTheta = aaTheta;
      rep.acceptedStep = (acceptedU - previousU).cwiseAbs().maxCoeff();
      rep.minJ = acceptedAdm.minJ;
      rep.maxQRel = acceptedAdm.maxQ;

      // ---- Diagnostics + stopping ----
      const auto surf = acceptedSurf;
      const Real eNow = acceptedEnergy;
      rep.activeRMS = surf.activeRMS;
      rep.activeSup = surf.activeSup;
      rep.activeFraction = surf.totalLen > Real(0)
        ? surf.activeLen / surf.totalLen : Real(0);
      rep.energy = eNow;
      if (p.trace)
        std::cout << "      wngir it=" << std::setw(3) << rep.iterations
                  << "  E=" << std::scientific << std::setprecision(3)
                  << eNow
                  << "  actRMS=" << surf.activeRMS
                  << "  actRMS/h=" << (h > Real(0) ? surf.activeRMS / h : Real(0))
                  << "  actSup=" << surf.activeSup
                  << "  actFrac=" << rep.activeFraction
                  << "  dΓ/h=" << (h > Real(0) ? rawDemonsRMS / h : Real(0))
                  << "  vΓ/h=" << (h > Real(0) ? liftedTraceRMS / h : Real(0))
                  << "  β=" << beta
                  << "  βvΓ/h=" << (h > Real(0) ? scaledTraceRMS / h : Real(0))
                  << "  step/h=" << (h > Real(0) ? rep.acceptedStep / h : Real(0))
                  << "  aa=" << (aaAccepted ? 1 : 0)
                  << "  aaTry=" << (aaTried ? 1 : 0)
                  << "  aaθ=" << aaTheta
                  << "  cgIt=" << rep.linearIterations
                  << "  cgErr=" << rep.linearError
                  << "  α=" << alpha
                  << "  bt=" << backtracks
                  << "  min_j=" << rep.minJ
                  << "  max_Q=" << rep.maxQRel
                  << '\n';

      if (surf.activeRMS <= activeRMSTol)
      {
        rep.exitReason = "active-rms-converged";
        ++rep.iterations;
        break;
      }
      if (surf.activeSup <= activeSupTol)
      {
        rep.exitReason = "active-sup-converged";
        ++rep.iterations;
        break;
      }
      const Real eRel =
        std::abs(ePrev - eNow) / std::max(ePrev, Real(1e-30));
      if (eRel < p.energyStagTol)
      {
        rep.exitReason = "energy-stagnation";
        ++rep.iterations;
        break;
      }
      ePrev = eNow;
    }
    return rep;
  }

  template <class Displacement>
  class WNGIR
  {
      using FESType = std::decay_t<
        decltype(std::declval<Displacement&>().getFiniteElementSpace())>;
      using TrialType = std::decay_t<
        decltype(Variational::TrialFunction(std::declval<const FESType&>()))>;
      using TestType = std::decay_t<
        decltype(Variational::TestFunction(std::declval<const FESType&>()))>;
      using ProblemType = std::decay_t<decltype(Variational::Problem(
          std::declval<TrialType&>(), std::declval<TestType&>()))>;

    public:
      explicit WNGIR(Displacement& u)
        : m_u(&u),
          m_duStep(u.getFiniteElementSpace()),
          m_vStep(u.getFiniteElementSpace()),
          m_duBulk(u.getFiniteElementSpace()),
          m_vBulk(u.getFiniteElementSpace()),
          m_stepProblem(m_duStep, m_vStep),
          m_bulkProblem(m_duBulk, m_vBulk)
      {}

      void setParameters(const WNGIRParameters& parameters)
      {
        m_parameters = parameters;
      }

      const WNGIRParameters& getParameters() const
      {
        return m_parameters;
      }

      const WNGIRReport& getReport() const
      {
        return m_report;
      }

      template <class Mesh, class PhiDerived, class GradDerived>
      WNGIRReport solve(
          const Mesh& mesh,
          const std::vector<Rodin::Index>& interfaceFacets,
          const Variational::RealFunctionBase<PhiDerived>& phi,
          const Variational::VectorFunctionBase<Real, GradDerived>& grad)
      {
        m_report = solveWNGIR(
            mesh, *m_u, interfaceFacets, phi, grad, m_parameters,
            m_duStep, m_vStep, m_duBulk, m_vBulk,
            m_stepProblem, m_bulkProblem);
        return m_report;
      }

    private:
      Displacement* m_u;
      TrialType m_duStep;
      TestType m_vStep;
      TrialType m_duBulk;
      TestType m_vBulk;
      ProblemType m_stepProblem;
      ProblemType m_bulkProblem;
      WNGIRParameters m_parameters;
      WNGIRReport m_report;
  };

  template <class Displacement>
  WNGIR(Displacement&) -> WNGIR<Displacement>;
}

#endif
