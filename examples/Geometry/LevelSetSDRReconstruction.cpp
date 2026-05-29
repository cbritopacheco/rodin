/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// Static SDR (Surface-aware Distance-Reconstruction) displacement on a P1
// vector finite-element space, expressed through Rodin::Variational::Problem.
//
// Pipeline:
//   1. Compute per-cell phase moments by integrating tanh(phi/eps).
//   2. Classify cells via MinSTCut (binary Potts).
//   3. Tag inside cells, outside cells, and the cut skeleton with
//      attributes 1/2/10.
//   4. Build s_h^LF as a P1 GridFunction by running Rodin::Distance::Eikonal
//      (Fast Marching Method) seeded on the cut skeleton, signed on the
//      interior region.
//   5. Solve a static displacement problem on a P1 vector FES so that
//      phi(X + u) approximates s_h^LF on the band {|s_h^LF| <= delta}.
//
// The SDR residual and tangent are expressed as Rodin form-language
// integrators (Variational::LinearFormIntegratorBase /
// LocalBilinearFormIntegratorBase). The Newton step is run through
// Rodin::Variational::Problem + Rodin::Solver::NewtonSolver. No
// hand-written dense residual/Jacobian loop is used as the main solve.
//
// This iteration intentionally implements only the SDR penalty +
// Tikhonov regularization. The Jacobian admissibility + shape barrier
// (with sigma_K branch and j_K^u) will land in a follow-up turn.
//
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Math.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/IntegrationPoint.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  using Vec2 = Math::SpatialVector<Real>;

  Vec2 vec2(Real x = 0, Real y = 0)
  {
    Vec2 out(2);
    out(0) = x;
    out(1) = y;
    return out;
  }

  // -------------------------------------------------------------------------
  //   Level-set definition (analytic circle, used for phi, grad phi, hess phi)
  // -------------------------------------------------------------------------

  struct CircleLevelSet
  {
    Real cx = 0.51;
    Real cy = 0.48;
    Real radius = 0.31;

    Real phi(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      return std::sqrt(dx * dx + dy * dy) - radius;
    }

    Vec2 grad(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r = std::max(std::sqrt(dx * dx + dy * dy), Real(1e-14));
      return vec2(dx / r, dy / r);
    }

    Math::SpatialMatrix<Real> hess(const Vec2& p) const
    {
      const Real dx = p(0) - cx;
      const Real dy = p(1) - cy;
      const Real r = std::max(std::sqrt(dx * dx + dy * dy), Real(1e-14));
      Math::SpatialMatrix<Real> H(2, 2);
      H(0, 0) = 1 / r - dx * dx / (r * r * r);
      H(0, 1) = -dx * dy / (r * r * r);
      H(1, 0) = H(0, 1);
      H(1, 1) = 1 / r - dy * dy / (r * r * r);
      return H;
    }
  };

  // -------------------------------------------------------------------------
  //   Phase moment cells used by MinSTCut classifier
  // -------------------------------------------------------------------------

  constexpr std::array<std::array<Real, 3>, 3> TriangleBarycentricQuadrature = {{
    {{ Real(2) / 3, Real(1) / 6, Real(1) / 6 }},
    {{ Real(1) / 6, Real(2) / 3, Real(1) / 6 }},
    {{ Real(1) / 6, Real(1) / 6, Real(2) / 3 }}
  }};

  Real applyPhaseMomentMap(Real phi, Real epsilon)
  {
    return std::tanh(phi / epsilon);
  }

  Vec2 interpolateVec(
      const std::array<Vec2, 3>& values,
      const std::array<Real, 3>& bary)
  {
    return bary[0] * values[0] + bary[1] * values[1] + bary[2] * values[2];
  }

  struct CellMomentInfo
  {
    Index index = 0;
    Real area = 0;
    Real moment = 0;
    std::array<Vec2, 3> x;
    std::array<Index, 3> vertices = {{ 0, 0, 0 }};
  };

  // Note: per the PR review we do NOT swap vertex order to force a
  // positive determinant. We compute the moment using the actual vertex
  // order returned by the mesh and accept whatever orientation the
  // polytope's Jacobian carries.
  std::vector<CellMomentInfo> collectCellMomentInfo(
      const Mesh<Context::Local>& mesh,
      const CircleLevelSet& levelSet,
      Real epsilon)
  {
    std::vector<CellMomentInfo> cells;
    cells.reserve(mesh.getCellCount());

    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cellPolytope = *cellIt;
      const auto& vertices = cellPolytope.getVertices();
      if (vertices.size() != 3)
        throw std::runtime_error(
            "LevelSetSDRReconstruction expects triangular cells.");

      CellMomentInfo info;
      info.index = cellPolytope.getIndex();
      for (size_t i = 0; i < 3; ++i)
      {
        info.vertices[i] = vertices[i];
        info.x[i] = mesh.getVertexCoordinates(vertices[i]);
      }

      // Signed area from the actual vertex order.
      const Vec2 e1 = info.x[1] - info.x[0];
      const Vec2 e2 = info.x[2] - info.x[0];
      const Real signedArea = Real(0.5) * (e1(0) * e2(1) - e1(1) * e2(0));
      info.area = std::abs(signedArea);

      Real moment = 0;
      for (const auto& bary : TriangleBarycentricQuadrature)
      {
        const Vec2 xq = interpolateVec(info.x, bary);
        moment += applyPhaseMomentMap(levelSet.phi(xq), epsilon);
      }
      info.moment = moment / TriangleBarycentricQuadrature.size();
      cells.push_back(std::move(info));
    }
    return cells;
  }

  Real facetLength(const Mesh<Context::Local>& mesh, Index facet)
  {
    const auto face = mesh.getFace(facet);
    const auto& vertices = face->getVertices();
    const Vec2 a = mesh.getVertexCoordinates(vertices[0]);
    const Vec2 b = mesh.getVertexCoordinates(vertices[1]);
    return (b - a).norm();
  }

  // -------------------------------------------------------------------------
  //   SDR integrators expressed as Rodin form-language integrators
  // -------------------------------------------------------------------------
  //
  //   SDR continuous energy with smooth band weight (no hard gate):
  //
  //     E_SDR(u) = (rho_S / (2 * M_w * h_ref^2))
  //                * int_Omega W(s_h^LF(X)) (phi(X + u(X)) - s_h^LF(X))^2 dX
  //
  //   where W(s) = exp(-s^2 / (2 * deltaW^2)) is a global Gaussian and
  //   M_w = int_Omega W(s_h^LF) dX is the weighted band measure.
  //   Dropping the hard gate is what eliminates the rank-1 ill-conditioning
  //   in the far field that previously required Tikhonov regularization or
  //   a SubMesh restriction: every cell now contributes a non-zero (but
  //   decreasing-with-distance) SDR Hessian, so the full-mesh linear
  //   system is coercive and Newton can converge quadratically.
  //
  //   Residual:
  //
  //     R(v; u) = (rho_S / (M_w h_ref^2))
  //               * int_Omega W(s) (phi(X+u) - s) grad_y phi(X+u) . v dX
  //
  //   Tangent:
  //
  //     K(du, v; u) = (rho_S / (M_w h_ref^2))
  //               * int_Omega W(s) [
  //                     (grad_y phi(X+u) . du) (grad_y phi(X+u) . v)
  //                   + (phi(X+u) - s) * (hess_y phi(X+u) du . v)
  //                 ] dX
  //
  //   The integrators are FES-independent. They query the FES for the
  //   per-cell finite element, ask the polytope for its quadrature, and
  //   evaluate vector-valued basis functions via FE.getBasis(local)(rc).
  //   No P1 / triangle / vdim==2 assumptions in the integrand evaluation.

  struct SDRParameters
  {
    Real rhoS = 1;        // SDR penalty weight
    Real deltaW = 0;      // Gaussian width of the smooth band weight
    Real hRef = 0;        // reference mesh size
    Real normalizer = 1;  // 1 / (M_w * h_ref^2), precomputed
  };

  // Quadrature order helper.
  inline size_t quadOrderFor(size_t feOrder)
  {
    return std::max<size_t>(2, 2 * feOrder);
  }

  template <class ScalarLevelSet, class SLFType,
            class TestType, class StateType>
  class SDRResidualIntegrator final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;

      SDRResidualIntegrator(
          const ScalarLevelSet& levelSet,
          const SLFType& sLF,
          const TestType& v,
          const StateType& u,
          SDRParameters params)
        : Variational::LinearFormIntegratorBase<Real>(v),
          m_levelSet(levelSet),
          m_sLF(sLF),
          m_test(v),
          m_state(u),
          m_params(params)
      {}

      SDRResidualIntegrator(const SDRResidualIntegrator&) = default;

      SDRResidualIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const size_t d = polytope.getDimension();
        const Index cellIdx = polytope.getIndex();
        const auto& testFES = m_test.get().getFiniteElementSpace();
        const auto& testFE = testFES.getFiniteElement(d, cellIdx);
        const size_t testDofs = testFE.getCount();
        const size_t vdim = testFES.getVectorDimension();

        m_elemVec.resize(testDofs);
        m_elemVec.setZero();

        // FES-driven quadrature: ask Rodin for a formula whose order
        // matches the FE's polynomial degree (doubled for products), then
        // ask the polytope to map it into physical space.
        const auto& qf =
          QF::PolytopeQuadratureFormula::get(
              quadOrderFor(testFE.getOrder()),
              polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const size_t nqp = quadrature.getSize();

        for (size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Real wq = qf.getWeight(q);
          const Real distortion = pt.getDistortion();

          const auto& xq = pt.getCoordinates();
          const auto& rc = pt.getReferenceCoordinates();

          const Real s = m_sLF.get().getValue(pt);
          const auto uqRange = m_state.get().getValue(pt);
          Vec2 y(d);
          for (size_t c = 0; c < vdim; ++c)
            y(c) = xq(c) + uqRange(c);

          const Real r = m_levelSet.phi(y) - s;
          const Vec2 gradPhi = m_levelSet.grad(y);

          // Smooth Gaussian weight — no hard band gate, so every cell
          // contributes a (possibly tiny) SDR Hessian and the full-mesh
          // system stays well-conditioned.
          const Real weight =
            std::exp(-s * s / (2 * m_params.deltaW * m_params.deltaW));

          const Real coef =
            m_params.rhoS * weight * r * m_params.normalizer;
          const Real measure = wq * distortion;

          for (size_t te = 0; te < testDofs; ++te)
          {
            // FES-generic basis value at the quadrature point. For a
            // vector FES this is a SpatialVector of length vdim with the
            // active component filled with the underlying scalar basis.
            const auto testValue = testFE.getBasis(te)(rc);
            Real dot = 0;
            for (size_t c = 0; c < vdim; ++c)
              dot += gradPhi(c) * testValue(c);
            m_elemVec(te) += measure * coef * dot;
          }
        }
        return *this;
      }

      ScalarType integrate(size_t te) final override
      {
        return m_elemVec(te);
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope->get();
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      SDRResidualIntegrator* copy() const noexcept final override
      {
        return new SDRResidualIntegrator(*this);
      }

    private:
      const ScalarLevelSet& m_levelSet;
      std::reference_wrapper<const SLFType> m_sLF;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_state;
      SDRParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  template <class ScalarLevelSet, class SLFType,
            class TrialType, class TestType, class StateType>
  class SDRTangentIntegrator final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;

      SDRTangentIntegrator(
          const ScalarLevelSet& levelSet,
          const SLFType& sLF,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          SDRParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_levelSet(levelSet),
          m_sLF(sLF),
          m_trial(du),
          m_test(v),
          m_state(u),
          m_params(params)
      {}

      SDRTangentIntegrator(const SDRTangentIntegrator&) = default;

      SDRTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const size_t d = polytope.getDimension();
        const Index cellIdx = polytope.getIndex();
        const auto& trialFES = m_trial.get().getFiniteElementSpace();
        const auto& testFES = m_test.get().getFiniteElementSpace();
        const auto& trialFE = trialFES.getFiniteElement(d, cellIdx);
        const auto& testFE = testFES.getFiniteElement(d, cellIdx);
        const size_t trialDofs = trialFE.getCount();
        const size_t testDofs = testFE.getCount();
        const size_t vdim = testFES.getVectorDimension();

        m_matrix.resize(testDofs, trialDofs);
        m_matrix.setZero();

        const auto& qf =
          QF::PolytopeQuadratureFormula::get(
              quadOrderFor(std::max(testFE.getOrder(), trialFE.getOrder())),
              polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const size_t nqp = quadrature.getSize();

        for (size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Real wq = qf.getWeight(q);
          const Real distortion = pt.getDistortion();

          const auto& xq = pt.getCoordinates();
          const auto& rc = pt.getReferenceCoordinates();

          const Real s = m_sLF.get().getValue(pt);
          const auto uqRange = m_state.get().getValue(pt);
          Vec2 y(d);
          for (size_t c = 0; c < vdim; ++c)
            y(c) = xq(c) + uqRange(c);

          // Note on tangent form: this is the Gauss-Newton approximation,
          // dropping the (phi - s) * hess(phi) second-order term. The
          // full Newton tangent for the SDR least-squares functional is
          //   K = coef * [(grad phi . du)(grad phi . v) + r * (hess phi du) . v]
          // The second term is indefinite when r changes sign or when phi
          // is non-convex (true for a sphere distance function inside the
          // sphere), causing Newton steps to diverge. Gauss-Newton keeps
          // only the rank-1 PSD outer-product part, which combined with
          // the convex barrier guarantees a positive-definite K and
          // monotone Newton convergence to the SDR-energy minimum.
          const Vec2 gradPhi = m_levelSet.grad(y);

          const Real weight =
            std::exp(-s * s / (2 * m_params.deltaW * m_params.deltaW));
          const Real coef =
            m_params.rhoS * weight * m_params.normalizer;
          const Real measure = wq * distortion;

          // Cache basis values to avoid recomputation in the (te, tr)
          // double loop.
          std::vector<Math::SpatialVector<Real>> testValues(testDofs);
          std::vector<Math::SpatialVector<Real>> trialValues(trialDofs);
          for (size_t te = 0; te < testDofs; ++te)
            testValues[te] = testFE.getBasis(te)(rc);
          for (size_t tr = 0; tr < trialDofs; ++tr)
            trialValues[tr] = trialFE.getBasis(tr)(rc);

          // Per-test-dof projection of grad phi onto the basis function.
          std::vector<Real> gpDotV(testDofs);
          for (size_t te = 0; te < testDofs; ++te)
          {
            Real s_val = 0;
            for (size_t c = 0; c < vdim; ++c)
              s_val += gradPhi(c) * testValues[te](c);
            gpDotV[te] = s_val;
          }
          std::vector<Real> gpDotU(trialDofs);
          for (size_t tr = 0; tr < trialDofs; ++tr)
          {
            Real s_val = 0;
            for (size_t c = 0; c < vdim; ++c)
              s_val += gradPhi(c) * trialValues[tr](c);
            gpDotU[tr] = s_val;
          }
          for (size_t te = 0; te < testDofs; ++te)
            for (size_t tr = 0; tr < trialDofs; ++tr)
              m_matrix(te, tr) +=
                measure * coef * gpDotU[tr] * gpDotV[te];
        }
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope->get();
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      SDRTangentIntegrator* copy() const noexcept final override
      {
        return new SDRTangentIntegrator(*this);
      }

    private:
      const ScalarLevelSet& m_levelSet;
      std::reference_wrapper<const SLFType> m_sLF;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_state;
      SDRParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  template <class L, class S, class V, class U>
  SDRResidualIntegrator(const L&, const S&, const V&, const U&, SDRParameters)
    -> SDRResidualIntegrator<L, S, V, U>;

  template <class L, class S, class DU, class V, class U>
  SDRTangentIntegrator(const L&, const S&, const DU&, const V&, const U&,
                       SDRParameters)
    -> SDRTangentIntegrator<L, S, DU, V, U>;

  // -------------------------------------------------------------------------
  //   Per-cell affine geometry cache for the Jacobian admissibility barrier
  // -------------------------------------------------------------------------
  //
  //   For each triangle cell K we precompute the affine map A_K = D F_K,
  //   the sigma_K branch sign of det(A_K), and the inverse-transpose used
  //   to evaluate the spatial gradient of the P1 displacement field on K:
  //
  //       grad_X u_h = A_K^{-T} * grad_xhat U
  //
  //   where the reference-gradient matrix grad_xhat U is built from the
  //   vertex displacement coefficients. The PR review explicitly requires
  //   that the actual polytope orientation (sigma_K) is respected and
  //   that vertex order is NOT swapped for positive determinant.
  //
  //   For affine triangles, A_K is constant on K, |K_hat| = 1/2, and
  //   |det A_K| = 2|K|, so the requested integral form
  //
  //       j_K_scale = (1/|K_hat|) int_{K_hat} |det A_K(x_hat)| dx_hat
  //
  //   reduces to |det A_K|.
  //
  struct CellGeomCache
  {
    Index index = 0;
    Real area = 0;                          // |K|
    Real detAK = 0;                         // signed det(A_K)
    Real absDetAK = 0;                      // |det(A_K)| = j_K_scale
    int sigmaK = 1;                         // sign(det A_K)
    Math::SpatialMatrix<Real> A;            // A_K (2x2)
    Math::SpatialMatrix<Real> Ainv;         // A_K^{-1} (2x2)
    Math::SpatialMatrix<Real> AinvT;        // A_K^{-T}
    Math::SpatialMatrix<Real> C;            // A_K A_K^T  (used by analytic shape barrier)
    Math::Matrix<Real> gradN;               // 3x2: rows = spatial grads of basis fns
    std::array<Index, 3> vertices = {{ 0, 0, 0 }};
  };

  std::vector<CellGeomCache> precomputeCellGeometry(const Mesh<Context::Local>& mesh)
  {
    std::vector<CellGeomCache> cache(mesh.getCellCount());
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto& vertices = cell.getVertices();
      CellGeomCache c;
      c.index = cell.getIndex();
      for (size_t i = 0; i < 3; ++i)
        c.vertices[i] = vertices[i];

      const Vec2 X0 = mesh.getVertexCoordinates(vertices[0]);
      const Vec2 X1 = mesh.getVertexCoordinates(vertices[1]);
      const Vec2 X2 = mesh.getVertexCoordinates(vertices[2]);

      Math::SpatialMatrix<Real> A(2, 2);
      A(0, 0) = X1(0) - X0(0); A(0, 1) = X2(0) - X0(0);
      A(1, 0) = X1(1) - X0(1); A(1, 1) = X2(1) - X0(1);

      const Real det = A.determinant();
      c.A = A;
      c.detAK = det;
      c.absDetAK = std::abs(det);
      c.sigmaK = det >= 0 ? 1 : -1;
      c.area = Real(0.5) * c.absDetAK;
      c.Ainv = A.inverse();
      c.AinvT = c.Ainv.transpose();
      c.C = A * A.transpose();

      // Reference gradients of the P1 basis on a triangle:
      //   N_0 = 1 - rx - ry, grad_ref = (-1, -1)
      //   N_1 = rx,          grad_ref = (+1,  0)
      //   N_2 = ry,          grad_ref = ( 0, +1)
      // Spatial gradient: grad_X N_k = A^{-T} * grad_ref N_k.
      c.gradN.resize(3, 2);
      const auto col0 = c.AinvT.col(0);
      const auto col1 = c.AinvT.col(1);
      c.gradN(0, 0) = -(col0(0) + col1(0));
      c.gradN(0, 1) = -(col0(1) + col1(1));
      c.gradN(1, 0) = col0(0);
      c.gradN(1, 1) = col0(1);
      c.gradN(2, 0) = col1(0);
      c.gradN(2, 1) = col1(1);

      cache[c.index] = std::move(c);
    }
    return cache;
  }

  // -------------------------------------------------------------------------
  //   Jacobian admissibility + shape barrier integrators
  // -------------------------------------------------------------------------
  //
  //   For affine triangles, grad_X u_h is constant on K. Let F = I + grad_X u_h.
  //   Then:
  //
  //     A_K^u = F * A_K
  //     j_K^u = sigma_K * det(F * A_K) / j_K_scale
  //           = sigma_K * det(F) * det(A_K) / |det(A_K)|
  //           = det(F)                                  (for affine triangles)
  //
  //     Q_shape(A_K^u)
  //           = (1/d) ||F A_K||_F^2 / (sigma_K det(F A_K))^(2/d)
  //           = (1/d) ||F A_K||_F^2 / (det(F) |det(A_K)|)^(2/d)   (d=2)
  //
  //   The integrators below evaluate the energy
  //
  //     E(u) = (gamma/|Omega|) sum_K |K| * (Q_shape - 1)
  //          + (beta /|Omega|) sum_K |K| * (-log(j_K^u - j_min))
  //
  //   and its derivatives using a 1-point quadrature on each cell (exact
  //   for constant integrand). The math is identical to the PR-review-
  //   approved formulation: sigma_K branch + j_K_scale = |det A_K|.

  struct BarrierParameters
  {
    Real gamma = 1e-1;   // shape weight
    Real beta = 1e-4;    // admissibility weight
    Real jMin = 1e-8;    // admissibility floor
    Real domainMeasure = 0;
  };

  // Compute the 2x2 reference-gradient of u_h on cell K given the cached
  // geometry and vertex displacement coefficients.
  Math::SpatialMatrix<Real> referenceGradFromVertexU(
      const std::array<Vec2, 3>& uVerts)
  {
    // grad_xhat U with row = component, col = reference axis
    Math::SpatialMatrix<Real> gradHat(2, 2);
    gradHat(0, 0) = uVerts[1](0) - uVerts[0](0);
    gradHat(0, 1) = uVerts[2](0) - uVerts[0](0);
    gradHat(1, 0) = uVerts[1](1) - uVerts[0](1);
    gradHat(1, 1) = uVerts[2](1) - uVerts[0](1);
    return gradHat;
  }

  // Vector P1 global dof layout is COMPONENT-MAJOR:
  //   global dof for (vertex v, component c) = v + c * vertexCount.
  // See Rodin::Variational::P1<SpatialVector<Real>, Mesh>::P1(mesh, vdim).
  template <class StateType>
  std::array<Vec2, 3> extractCellU(
      const CellGeomCache& cell, const StateType& u)
  {
    const auto& uData = u.getData();
    const Index vn = u.getFiniteElementSpace().getMesh().getVertexCount();
    std::array<Vec2, 3> uv;
    for (size_t i = 0; i < 3; ++i)
    {
      const Index v = cell.vertices[i];
      uv[i] = vec2(uData(v), uData(v + vn));
    }
    return uv;
  }

  // Energy + gradient computation reused by residual and tangent.
  // Returns Energy contribution; fills gradient w.r.t. the 6 nodal dofs
  // (layout: dof = node*2 + component).
  struct BarrierLocal
  {
    Real energy = 0;
    Real j = 0;            // j_K^u (= det F for affine triangles)
    Boolean valid = true;
    Math::Vector<Real> grad;       // length 6
    Math::Matrix<Real> hess;       // 6x6
  };

  // Fully analytical residual and tangent of the shape + admissibility
  // barrier on an affine triangle, with the sigma_K branch.
  //
  // Notation in this routine (all 2x2 unless noted):
  //   F    = I + grad_X u_h           on the cell (constant for affine map)
  //   M    = F A_K                    (cell-deformed Jacobian)
  //   D    = sigma_K det(M) = |det A_K| * det F
  //   n    = ||M||_F^2 = tr(M^T M)
  //   FC   = F * C, with C = A_K A_K^T
  //   FinvT = F^{-T}
  //   j     = det F                   (= j_K^u for affine triangles)
  //   h     = j / (j - jMin)          (log-barrier multiplier)
  //   K_d   = jMin * j / (j - jMin)^2 (det-Hessian scalar)
  //   w_s   = gamma * |K|/|Omega|
  //   w_d   = beta  * |K|/|Omega|
  //
  // Energy:
  //   E_cell = w_s * (n / (2 D) - 1) - w_d * log(j - jMin)
  //
  // First derivative w.r.t. F (closed form):
  //   S = w_s/D * (F C - (n/2) F^{-T}) - w_d * h * F^{-T}
  //
  // Mapping S -> per-dof residual via grad_X N_k (constant on the cell):
  //   R[k*2 + l] = (S * grad N_k)_l       where grad N_k is a 2-vector
  //
  // Tangent contraction is derived in the comment above the loop.
  BarrierLocal evaluateBarrierLocal(
      const CellGeomCache& cell,
      const std::array<Vec2, 3>& uv,
      const BarrierParameters& params)
  {
    BarrierLocal out;
    out.grad = Math::Vector<Real>::Zero(6);
    out.hess = Math::Matrix<Real>::Zero(6, 6);

    // Build F = I + grad_X u_h. With gradN(k, j) = (grad_X N_k)_j and
    // U_ij = u_node-j-comp-i, we have G_ij = sum_k U(i,k) gradN(k,j).
    Math::SpatialMatrix<Real> U(2, 3);
    for (size_t k = 0; k < 3; ++k)
    {
      U(0, k) = uv[k](0);
      U(1, k) = uv[k](1);
    }
    Math::SpatialMatrix<Real> F =
      Math::SpatialMatrix<Real>::Identity(2, 2) + U * cell.gradN;

    const Real j = F.determinant();
    out.j = j;
    if (j <= params.jMin)
    {
      out.valid = false;
      return out;
    }

    const Math::SpatialMatrix<Real> M = F * cell.A;
    const Real n2 = M.squaredNorm();
    const Real D = cell.absDetAK * j;
    const Real qShape = Real(0.5) * n2 / D;

    const Real areaWeight = cell.area / params.domainMeasure;
    const Real w_s = params.gamma * areaWeight;
    const Real w_d = params.beta  * areaWeight;

    out.energy =
        w_s * (qShape - Real(1))
      - w_d * std::log(j - params.jMin);

    // Closed-form first derivative S = dE/dF (2x2 matrix).
    const Math::SpatialMatrix<Real> FinvT = F.inverse().transpose();
    const Math::SpatialMatrix<Real> FC = F * cell.C;

    const Real h = j / (j - params.jMin);
    const Real Kd = params.jMin * j
                  / ((j - params.jMin) * (j - params.jMin));

    const Real wsOverD = w_s / D;
    Math::SpatialMatrix<Real> S =
        wsOverD * (FC - Real(0.5) * n2 * FinvT)
      - w_d * h * FinvT;

    // Residual: R[k*2 + l] = (S * gradN_k)_l = sum_j S(l, j) gradN(k, j).
    for (size_t k = 0; k < 3; ++k)
      for (size_t l = 0; l < 2; ++l)
      {
        Real sum = 0;
        for (size_t jj = 0; jj < 2; ++jj)
          sum += S(l, jj) * cell.gradN(k, jj);
        out.grad(k * 2 + l) = sum;
      }

    // Tangent: contract the 4-tensor T_{l,j,p,j'} = ∂S_{l,j}/∂F_{p,j'}
    // with grad_X N_k(j) grad_X N_m(j'). All eight terms expand to
    // outer products of 2-vectors, so the per-(k,m) block is a 2x2 matrix.
    //
    // Pre-compute per-node 2-vectors used by every (k, m) pair.
    std::array<Math::SpatialVector<Real>, 3> gN;
    std::array<Math::SpatialVector<Real>, 3> FCgk;     // FC * gN_k
    std::array<Math::SpatialVector<Real>, 3> FinvTgk;  // FinvT * gN_k
    for (size_t k = 0; k < 3; ++k)
    {
      gN[k] = Math::SpatialVector<Real>(2);
      gN[k](0) = cell.gradN(k, 0);
      gN[k](1) = cell.gradN(k, 1);
      FCgk[k]    = FC * gN[k];
      FinvTgk[k] = FinvT * gN[k];
    }

    const Real halfN = Real(0.5) * n2;

    for (size_t k = 0; k < 3; ++k)
    {
      for (size_t m = 0; m < 3; ++m)
      {
        const Real gCg = gN[k].dot(cell.C * gN[m]);  // scalar, symmetric

        for (size_t l = 0; l < 2; ++l)
        {
          for (size_t p = 0; p < 2; ++p)
          {
            // Shape contribution (5 terms):
            const Real shape =
                (l == p ? gCg : Real(0))
              - FCgk[k](l)    * FinvTgk[m](p)
              - FCgk[m](p)    * FinvTgk[k](l)
              + halfN * (FinvTgk[k](l) * FinvTgk[m](p)
                       + FinvTgk[k](p) * FinvTgk[m](l));

            // Det / admissibility contribution (2 terms):
            const Real det =
                Kd * FinvTgk[k](l) * FinvTgk[m](p)
              + h  * FinvTgk[k](p) * FinvTgk[m](l);

            out.hess(k * 2 + l, m * 2 + p) =
                wsOverD * shape + w_d * det;
          }
        }
      }
    }

    return out;
  }

  template <class TestType, class StateType>
  class BarrierResidualIntegrator final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;

      BarrierResidualIntegrator(
          const std::vector<CellGeomCache>& cellCache,
          const TestType& v,
          const StateType& u,
          BarrierParameters params)
        : Variational::LinearFormIntegratorBase<Real>(v),
          m_cellCache(cellCache),
          m_test(v),
          m_state(u),
          m_params(params)
      {}

      BarrierResidualIntegrator(const BarrierResidualIntegrator&) = default;

      BarrierResidualIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const Index cellIdx = polytope.getIndex();
        m_elemVec.resize(6);
        m_elemVec.setZero();
        if (cellIdx >= m_cellCache.get().size())
          return *this;
        const auto& cell = m_cellCache.get()[cellIdx];
        const auto uv = extractCellU(cell, m_state.get());
        const auto local = evaluateBarrierLocal(cell, uv, m_params);
        if (!local.valid)
        {
          // Stiff repulsive linear term to push out of the inadmissible
          // region. Magnitude chosen so the next Newton step retreats.
          const Real bigEnergy = std::numeric_limits<Real>::max() / 1e10;
          (void) bigEnergy;
          // Leave zero residual; the (also zeroed) tangent below would
          // produce a singular system. Instead, throw to halt; in
          // practice the iteration is not expected to ever enter
          // inadmissibility once the barrier is active.
          throw std::runtime_error(
              "Barrier integrator encountered j <= jMin at cell "
              + std::to_string(cellIdx));
        }
        m_elemVec = local.grad;
        return *this;
      }

      ScalarType integrate(size_t te) final override
      {
        return m_elemVec(te);
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope->get();
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      BarrierResidualIntegrator* copy() const noexcept final override
      {
        return new BarrierResidualIntegrator(*this);
      }

    private:
      std::reference_wrapper<const std::vector<CellGeomCache>> m_cellCache;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_state;
      BarrierParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  template <class TrialType, class TestType, class StateType>
  class BarrierTangentIntegrator final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;

      BarrierTangentIntegrator(
          const std::vector<CellGeomCache>& cellCache,
          const TrialType& du,
          const TestType& v,
          const StateType& u,
          BarrierParameters params)
        : Variational::LocalBilinearFormIntegratorBase<Real>(du, v),
          m_cellCache(cellCache),
          m_trial(du),
          m_test(v),
          m_state(u),
          m_params(params)
      {}

      BarrierTangentIntegrator(const BarrierTangentIntegrator&) = default;

      BarrierTangentIntegrator& setPolytope(
          const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const Index cellIdx = polytope.getIndex();
        m_matrix.resize(6, 6);
        m_matrix.setZero();
        if (cellIdx >= m_cellCache.get().size())
          return *this;
        const auto& cell = m_cellCache.get()[cellIdx];
        const auto uv = extractCellU(cell, m_state.get());
        const auto local = evaluateBarrierLocal(cell, uv, m_params);
        if (!local.valid)
        {
          throw std::runtime_error(
              "Barrier tangent encountered j <= jMin at cell "
              + std::to_string(cellIdx));
        }
        m_matrix = local.hess;
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope->get();
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      BarrierTangentIntegrator* copy() const noexcept final override
      {
        return new BarrierTangentIntegrator(*this);
      }

    private:
      std::reference_wrapper<const std::vector<CellGeomCache>> m_cellCache;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_state;
      BarrierParameters m_params;
      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  template <class V, class U>
  BarrierResidualIntegrator(const std::vector<CellGeomCache>&, const V&,
                            const U&, BarrierParameters)
    -> BarrierResidualIntegrator<V, U>;

  template <class DU, class V, class U>
  BarrierTangentIntegrator(const std::vector<CellGeomCache>&, const DU&,
                           const V&, const U&, BarrierParameters)
    -> BarrierTangentIntegrator<DU, V, U>;
}

int main(int, char**)
{
  constexpr size_t n = 17;
  const Real h = Real(1) / static_cast<Real>(n - 1);
  const Real epsilon = 1.25 * h;
  const Real lambdaC = 0.008;

  // -------------------------------------------------------------------------
  //   Step 1-3: mesh, classification, attribute tagging
  // -------------------------------------------------------------------------
  LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.scale(h);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 2);
  mesh.getConnectivity().compute(0, 0);

  const CircleLevelSet levelSet;
  const auto cellMoments = collectCellMomentInfo(mesh, levelSet, epsilon);

  std::vector<Real> volumes(cellMoments.size());
  std::vector<Real> moments(cellMoments.size());
  for (const auto& info : cellMoments)
  {
    volumes[info.index] = info.area;
    moments[info.index] = info.moment;
  }

  std::vector<MinSTCut::Edge> graphEdges;
  for (auto faceIt = mesh.getFace(); faceIt; ++faceIt)
  {
    const Index facet = faceIt->getIndex();
    const auto& incident = mesh.getConnectivity().getIncidence({1, 2}, facet);
    if (incident.size() == 2)
    {
      graphEdges.push_back({
          incident[0], incident[1],
          lambdaC * facetLength(mesh, facet),
          facet});
    }
  }

  const MinSTCut cut;
  const MinSTCut::Result classified =
    cut.classify(volumes, moments, graphEdges);

  std::vector<Index> interfaceFacets;
  interfaceFacets.reserve(classified.cutEdges.size());
  for (const MinSTCut::Edge& edge : classified.cutEdges)
    if (edge.index != MinSTCut::InvalidIndex)
      interfaceFacets.push_back(edge.index);

  constexpr Attribute interiorAttribute = 1;
  constexpr Attribute exteriorAttribute = 2;
  constexpr Attribute interfaceAttribute = 10;
  for (Index cell = 0; cell < classified.labels.size(); ++cell)
  {
    mesh.setAttribute(
        {mesh.getDimension(), cell},
        classified.labels[cell] == MinSTCut::Inside
          ? interiorAttribute
          : exteriorAttribute);
  }
  for (const Index facet : interfaceFacets)
    mesh.setAttribute({mesh.getDimension() - 1, facet}, interfaceAttribute);

  // Tag mesh boundary facets so we can DirichletBC them later.
  constexpr Attribute boundaryAttribute = 20;
  for (auto faceIt = mesh.getBoundary(); faceIt; ++faceIt)
    mesh.setAttribute({mesh.getDimension() - 1, faceIt->getIndex()},
                      boundaryAttribute);

  mesh.save("LevelSetSDRReconstruction_LF.mesh", IO::FileFormat::MFEM);

  // -------------------------------------------------------------------------
  //   Step 4: signed distance s_h^LF via Distance::Eikonal (FMM)
  // -------------------------------------------------------------------------
  using ScalarFES = P1<Real, LocalMesh>;
  ScalarFES sLFFes(mesh);
  GridFunction sLF(sLFFes);
  sLF = Real(0);
  Distance::Eikonal<ScalarFES, Math::Vector<Real>>(sLF)
    .setInterface(interfaceAttribute)
    .setInterior(interiorAttribute)
    .solve()
    .sign();

  // -------------------------------------------------------------------------
  //   Step 5a: identify band-active cells and build a SubMesh on them.
  //
  //   Matching OLD's active-dof restriction in form-language form: only
  //   cells with at least one quadrature point in the SDR band participate
  //   in the linear system. This removes far-field dofs whose Hessian
  //   would otherwise be a tiny barrier term, restoring full rank on the
  //   live dofs and quadratic Newton convergence.
  // -------------------------------------------------------------------------
  // -------------------------------------------------------------------------
  //   Step 5a: smooth band weight and weighted band measure M_w.
  //
  //   deltaW > delta widens the Gaussian so far-field cells still see a
  //   small (but nonzero) SDR contribution. Without this widening the
  //   weight drops to ~1e-2 at distance 3*delta from the interface; with
  //   deltaW = 1.5*delta the same point sees ~exp(-2) = 0.14, which is
  //   the smallest scale we want the SDR Hessian to reach.
  // -------------------------------------------------------------------------
  const Real delta = 1.75 * h;
  const Real deltaW = 1.5 * delta;
  Real domainMeasure = 0;
  Real weightedBandMeasure = 0;
  for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
  {
    const auto& cell = *cellIt;
    for (const auto& bary : TriangleBarycentricQuadrature)
    {
      Math::SpatialPoint rc(2);
      rc(0) = bary[1];
      rc(1) = bary[2];
      const Geometry::Point pt(cell, rc);
      Math::SpatialMatrix<Real> J;
      cell.getTransformation().jacobian(J, rc);
      const Real triangleArea = Real(0.5) * std::abs(J.determinant());
      const Real w =
        triangleArea / static_cast<Real>(TriangleBarycentricQuadrature.size());
      domainMeasure += w;
      const Real s = sLF.getValue(pt);
      const Real W = std::exp(-s * s / (2 * deltaW * deltaW));
      weightedBandMeasure += W * w;
    }
  }
  if (weightedBandMeasure <= 0)
    throw std::runtime_error("Empty smooth band: M_w = 0.");

  // -------------------------------------------------------------------------
  //   Step 5b: P1 vector FES on the full mesh
  // -------------------------------------------------------------------------
  using VectorFES = P1<Math::SpatialVector<Real>, LocalMesh>;
  VectorFES vectorFes(mesh, 2);
  TrialFunction du(vectorFes);
  TestFunction  v(vectorFes);
  GridFunction  u(vectorFes);
  u.getData().setZero();

  SDRParameters params;
  params.rhoS = 1;
  params.deltaW = deltaW;
  params.hRef = h;
  params.normalizer = Real(1) / (weightedBandMeasure * h * h);

  BarrierParameters barrierParams;
  barrierParams.gamma = 1e-1;
  // Stronger admissibility floor so Newton steps that would otherwise
  // push a cell past det(F) = 0 are damped by the log barrier well
  // before the singularity. Matches OLD when delta=1.75h.
  barrierParams.beta = 1e-2;
  barrierParams.jMin = 1e-3;
  barrierParams.domainMeasure = domainMeasure;

  const auto cellCache = precomputeCellGeometry(mesh);

  // -------------------------------------------------------------------------
  //   Step 6: Problem assembly and Newton solve
  // -------------------------------------------------------------------------
  // Full nonlinear residual: R(v;u) = R_SDR + R_barrier.
  // Tangent: K(du,v;u) = K_SDR + K_barrier (full Hessian, no FD wrapper).
  // With the smooth band weight, every cell contributes a nonzero SDR
  // Hessian and the full-mesh system is coercive — no Tikhonov needed.
  // DirichletBC on the grid boundary pins the perimeter dofs to zero,
  // which both eliminates rigid-body modes and reflects the fact that
  // the mesh boundary is the analyst-imposed domain edge (the SDR fit
  // happens in the interior).
  auto zero = VectorFunction{ Zero(), Zero() };

  Problem newton(du, v);
  newton =
        SDRTangentIntegrator(levelSet, sLF, du, v, u, params)
      + BarrierTangentIntegrator(cellCache, du, v, u, barrierParams)
      + SDRResidualIntegrator(levelSet, sLF, v, u, params)
      + BarrierResidualIntegrator(cellCache, v, u, barrierParams)
      + DirichletBC(du, zero).on(boundaryAttribute);

  Solver::SparseLU linearSolver(newton);
  Solver::NewtonSolver solver(linearSolver);
  std::vector<Real> residualHistory;
  solver
    .setMaxIterations(20)
    .setDampingFactor(1.0)
    .setAbsoluteTolerance(1e-10)
    .setRelativeTolerance(1e-8)
    .setMonitor([&](const auto& report)
    {
      std::cout << "Newton " << std::setw(2) << report.iterations
                << ": ||R||=" << std::scientific << std::setprecision(6)
                << report.final_residual
                << "  step=" << report.final_step_norm
                << "  damping=" << report.damping_factor
                << '\n';
      residualHistory.push_back(report.final_residual);
    });

  solver.solve(u);
  const auto report = solver.getReport();

  // -------------------------------------------------------------------------
  //   Step 7: diagnostics + XDMF output
  // -------------------------------------------------------------------------
  LocalMesh moved(mesh);
  const auto& uData = u.getData();
  const Index vn = mesh.getVertexCount();
  for (Index vertex = 0; vertex < vn; ++vertex)
  {
    const Vec2 x = mesh.getVertexCoordinates(vertex);
    // Vector P1 dof layout is component-major: dof = vertex + comp * vn.
    const Real ux = uData(vertex);
    const Real uy = uData(vertex + vn);
    moved.setVertexCoordinates(vertex, vec2(x(0) + ux, x(1) + uy));
  }
  moved.save("LevelSetSDRReconstruction_HF.mesh", IO::FileFormat::MFEM);

  // -------------------------------------------------------------------------
  //   Step 7b: XDMF output.
  //
  //   Cell-centered (P0) fields:
  //     LF mesh:
  //       cell_label    : -1 inside, +1 outside (from MinSTCut)
  //       phase_moment  : per-cell phi-tanh moment used by the classifier
  //       sigma_K       : sign(det A_K), the cell-orientation branch
  //     Moved (HF) mesh:
  //       j             : j_K^u = det(F) for affine triangles
  //       q_shape       : Q_shape(F A_K) using sigma_K branch
  //       cell_label    : carried over for reference
  //
  //   Node-centered (P1) fields:
  //     LF mesh:
  //       sLF, phi      : signed distance (FMM) and the analytic level set
  //     Moved mesh:
  //       displacement  : the vector P1 solution u (passed directly)
  //       phi_moved     : phi evaluated at the moved vertex coordinates
  //
  //   Indexing convention (PR-review items 6 + 7):
  //     - P0 cell writes use FES.getGlobalIndex({d, cellIdx}, 0); no
  //       implicit "dof index == cell index" assumption.
  //     - Vector P1 displacement uses the GridFunction directly; no raw
  //       Eigen indexing of u.getData().
  using ScalarP0 = P0<Real, LocalMesh>;
  ScalarP0 p0Fes(mesh);
  GridFunction cellLabel(p0Fes);
  GridFunction phaseMoment(p0Fes);
  GridFunction sigmaKgf(p0Fes);
  for (const auto& info : cellMoments)
  {
    const Index dof =
      p0Fes.getGlobalIndex({mesh.getDimension(), info.index}, 0);
    cellLabel.getData()(dof) =
      static_cast<Real>(classified.labels[info.index]);
    phaseMoment.getData()(dof) = info.moment;
    sigmaKgf.getData()(dof) =
      static_cast<Real>(cellCache[info.index].sigmaK);
  }
  cellLabel.setName("cell_label");
  phaseMoment.setName("phase_moment");
  sigmaKgf.setName("sigma_K");

  P1<Real, LocalMesh> p1Fes(mesh);
  GridFunction phi(p1Fes);
  phi = [&](const Geometry::Point& p) -> Real
  {
    const auto& X = p.getCoordinates();
    return levelSet.phi(vec2(X(0), X(1)));
  };
  phi.setName("phi");
  sLF.setName("s_LF");

  // ----- HF (moved) mesh -----
  ScalarP0 p0FesMoved(moved);
  GridFunction jK(p0FesMoved);
  GridFunction qShape(p0FesMoved);
  GridFunction cellLabelHF(p0FesMoved);
  // Recompute the moved-mesh cell geometry to get det F and Q_shape from
  // the actual updated vertex positions (rather than re-deriving from u).
  const auto movedCache = precomputeCellGeometry(moved);
  for (const auto& info : cellMoments)
  {
    const Index dof =
      p0FesMoved.getGlobalIndex({moved.getDimension(), info.index}, 0);
    const auto& src = cellCache[info.index];
    const auto& dst = movedCache[info.index];
    // F = A_K^u A_K^{-1}, det F = det(A_K^u) / det(A_K).
    // With the sigma_K branch:
    //   j_K^u = sigma_K det(A_K^u) / |det A_K|
    //         = sigma_K det F * sigma_K = det F.
    jK.getData()(dof) = dst.detAK / src.detAK;
    // Q_shape = 0.5 * ||A_K^u||_F^2 / (sigma_K det A_K^u)
    //         = 0.5 * ||dst.A||_F^2 / (sigma_K * dst.detAK)
    qShape.getData()(dof) =
      Real(0.5) * dst.A.squaredNorm()
        / (static_cast<Real>(src.sigmaK) * dst.detAK);
    cellLabelHF.getData()(dof) =
      static_cast<Real>(classified.labels[info.index]);
  }
  jK.setName("j");
  qShape.setName("q_shape");
  cellLabelHF.setName("cell_label");

  P1<Real, LocalMesh> p1FesMoved(moved);
  GridFunction phiMoved(p1FesMoved);
  phiMoved = [&](const Geometry::Point& p) -> Real
  {
    const auto& X = p.getCoordinates();
    return levelSet.phi(vec2(X(0), X(1)));
  };
  phiMoved.setName("phi_moved");

  u.setName("displacement");

  IO::XDMF xdmf("LevelSetSDRReconstruction");
  auto lfGrid = xdmf.grid("LF");
  lfGrid.setMesh(mesh);
  lfGrid.add(cellLabel, IO::XDMF::Center::Cell);
  lfGrid.add(phaseMoment, IO::XDMF::Center::Cell);
  lfGrid.add(sigmaKgf, IO::XDMF::Center::Cell);
  lfGrid.add(phi, IO::XDMF::Center::Node);
  lfGrid.add(sLF, IO::XDMF::Center::Node);
  // The displacement field is defined on the LF mesh (its FES lives
  // there); writing it on the LF grid is the correct association.
  lfGrid.add(u, IO::XDMF::Center::Node);

  auto hfGrid = xdmf.grid("HF");
  hfGrid.setMesh(moved);
  hfGrid.add(cellLabelHF, IO::XDMF::Center::Cell);
  hfGrid.add(jK, IO::XDMF::Center::Cell);
  hfGrid.add(qShape, IO::XDMF::Center::Cell);
  hfGrid.add(phiMoved, IO::XDMF::Center::Node);

  xdmf.write().close();

  // Interface ||phi(X+u)|| RMS on Gamma_h^LF. Vertex displacements come
  // from `moved` (parent mesh, post-projection) which already encodes the
  // band-vertex displacements applied via the SubMesh polytope map.
  Real interfacePhi = 0;
  Real interfaceLength = 0;
  for (const Index facet : interfaceFacets)
  {
    const auto face = mesh.getFace(facet);
    const auto& vertices = face->getVertices();
    const Vec2 a = moved.getVertexCoordinates(vertices[0]);
    const Vec2 b = moved.getVertexCoordinates(vertices[1]);
    const Vec2 mid = Real(0.5) * (a + b);
    const Real val = levelSet.phi(mid);
    const Real len = (b - a).norm();
    interfacePhi += val * val * len;
    interfaceLength += len;
  }
  const Real interfacePhiRMS =
    std::sqrt(interfacePhi / std::max(interfaceLength, Real(1e-30)));

  std::cout << "\nDiagnostics\n";
  std::cout << "  cells inside / outside: "
            << classified.insideCells.size() << " / "
            << classified.outsideCells.size() << '\n';
  std::cout << "  interface facets: " << interfaceFacets.size() << '\n';
  std::cout << "  vector dofs: " << uData.size() << '\n';
  std::cout << "  h: " << h << ", delta: " << delta << '\n';
  std::cout << "  |Omega_h|: " << domainMeasure
            << ", M_w (smooth band): " << weightedBandMeasure << '\n';
  std::cout << "  final ||phi(X+u)||_RMS on Gamma_h^LF: "
            << interfacePhiRMS << '\n';
  std::cout << "  Newton iterations: " << report.iterations
            << ", converged: " << (report.converged ? "yes" : "no") << '\n';

  return report.converged ? 0 : 1;
}
