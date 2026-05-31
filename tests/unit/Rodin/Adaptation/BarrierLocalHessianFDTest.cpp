/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// Finite-difference verification of the analytical local Hessian
// returned by Rodin::Adaptation::evaluateBarrierLocal.
//
// For a correct analytical Hessian H and a central FD approximation of
// the gradient,
//
//   FD(eps; p) = ( grad(uv + eps p) - grad(uv - eps p) ) / (2 eps),
//
// the error
//
//   err(eps) = || FD(eps; p) - H p || / max(1, || H p ||)
//
// behaves like
//
//   err(eps) ~ C1 * eps^2 + C2 * eps_machine / eps
//
// for sufficiently smooth energy. Hence err(eps) DECREASES quadratically
// from eps = 1e-3 down to ~1e-5 (the truncation regime), reaches a
// minimum at ~1e-5, and INCREASES again like 1/eps as roundoff
// dominates. This V-shape is the canonical signature of a correct
// analytical derivative.
//
// The tests below exercise that signature on a small deterministic
// configuration: a single, uniform-grid triangle, with a small
// admissible displacement and a deterministic unit perturbation.
//
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include <Rodin/Adaptation.h>
#include <Rodin/Geometry.h>
#include <Rodin/Math/SpatialVector.h>
#include <Rodin/Math/Vector.h>

using namespace Rodin;
using namespace Rodin::Adaptation;
using namespace Rodin::Geometry;

namespace
{
  // A tiny mesh of triangles. UniformGrid produces a single right
  // triangle per cell with consistent orientation.
  Mesh<Context::Local> makeSmallTriangleMesh(std::size_t n = 4)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(
        Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(0, 0);
    return mesh;
  }

  Mesh<Context::Local> makeUnitSquareTriangleMesh(std::size_t n)
  {
    auto mesh = makeSmallTriangleMesh(n);
    mesh.scale(Real(1) / static_cast<Real>(n - 1));
    return mesh;
  }

  // Default scalar weights for the test (gamma and beta are passed
  // alongside `BarrierParameters` to `evaluateBarrierLocal` now that
  // they're form-language-friendly RealFunctionBase inputs at the
  // integrator boundary).
  constexpr Real kGamma = 1e-1;
  constexpr Real kBeta  = 1e-2;

  // For the FD test we INTENTIONALLY pick jSafe large enough that the
  // floor barrier is active at the test displacement (j ~ 1 < jSafe).
  // This exercises the active-branch closed-form Hessian; the inactive
  // branch is trivial (identically zero B-contribution, so only the
  // shape block is exercised — that block is unchanged from the previous
  // implementation).
  BarrierParameters defaultParams(Real domainMeasure)
  {
    BarrierParameters p;
    p.jMin = Real(0.1);
    p.jSafe = Real(1.5);
    p.domainMeasure = domainMeasure;
    return p;
  }

  // Deterministic small displacement of the three cell vertices.
  std::array<Math::SpatialVector<Real>, 3> deterministicUv(
      std::size_t seed, Real scale)
  {
    std::array<Math::SpatialVector<Real>, 3> uv;
    for (std::size_t k = 0; k < 3; ++k)
    {
      uv[k] = Math::SpatialVector<Real>(2);
      uv[k](0) = scale * std::sin(Real(1) + k + seed);
      uv[k](1) = scale * std::cos(Real(2) + k - seed);
    }
    return uv;
  }

  // Deterministic length-1 perturbation in the 6-dof local space.
  Math::Vector<Real> deterministicDir(std::size_t seed)
  {
    Math::Vector<Real> p(6);
    for (std::size_t i = 0; i < 6; ++i)
      p(i) = std::sin(Real(0.7) * (i + 1) + Real(seed));
    p /= p.norm();
    return p;
  }

  // Central-FD of grad(uv +- eps p), returned as a 6-vector.
  Math::Vector<Real> centralFD(
      const CellGeomCache& cell,
      const std::array<Math::SpatialVector<Real>, 3>& uv,
      const BarrierParameters& params,
      const Math::Vector<Real>& p,
      Real eps)
  {
    std::array<Math::SpatialVector<Real>, 3> uvPlus = uv;
    std::array<Math::SpatialVector<Real>, 3> uvMinus = uv;
    for (std::size_t k = 0; k < 3; ++k)
      for (std::size_t l = 0; l < 2; ++l)
      {
        uvPlus[k](l) += eps * p(k * 2 + l);
        uvMinus[k](l) -= eps * p(k * 2 + l);
      }
    const auto localPlus =
      evaluateBarrierLocal(cell, uvPlus, kGamma, kBeta, params);
    const auto localMinus =
      evaluateBarrierLocal(cell, uvMinus, kGamma, kBeta, params);
    return (localPlus.grad - localMinus.grad) / (2 * eps);
  }

  // Relative error || FD - H p || / max(1, || H p ||).
  Real relativeError(
      const Math::Vector<Real>& fd, const Math::Vector<Real>& Hp)
  {
    const Real denom = std::max<Real>(1, Hp.norm());
    return (fd - Hp).norm() / denom;
  }
}

namespace Rodin::Tests::Unit
{
  // Sanity-check that the cache itself was built and is admissible at
  // the test displacement.
  TEST(Rodin_Adaptation_Barrier, EvaluateAtSmallDisplacementIsValid)
  {
    auto mesh = makeSmallTriangleMesh();
    auto [cellCache, cellToLocal] = precomputeCellGeometry(mesh);
    ASSERT_FALSE(cellCache.empty());
    const auto& cell = cellCache[0];
    auto uv = deterministicUv(0, /*scale=*/Real(1e-3));
    const auto local = evaluateBarrierLocal(cell, uv, kGamma, kBeta, defaultParams(1));
    EXPECT_TRUE(local.valid);
    EXPECT_GT(local.j, Real(0));
  }

  TEST(Rodin_Adaptation_Barrier, IdentityJacobianRatioIsMeshIndependent)
  {
    const std::array<std::size_t, 4> ns{ 20, 30, 40, 60 };
    const Real jMinRatio = Real(1e-8);
    const Real jSafeRatio = Real(1e-3);
    const Real lineSearchSafetyMargin = Real(10);
    const Real jLineSearchRatio =
      std::max(jMinRatio, lineSearchSafetyMargin * jSafeRatio);

    std::cout << "\nJacobian-ratio refinement check\n";
    for (const std::size_t n : ns)
    {
      auto mesh = makeUnitSquareTriangleMesh(n);
      auto [cellCache, cellToLocal] = precomputeCellGeometry(mesh);
      ASSERT_FALSE(cellCache.empty());

      Real minJIdentity = std::numeric_limits<Real>::infinity();
      Real maxIdentityError = 0;
      for (const auto& cell : cellCache)
      {
        std::array<Math::SpatialVector<Real>, 3> uv;
        for (auto& x : uv)
        {
          x = Math::SpatialVector<Real>(2);
          x(0) = 0;
          x(1) = 0;
        }

        BarrierParameters params;
        params.jMin = jMinRatio;
        params.jSafe = jSafeRatio;
        params.domainMeasure = Real(1);
        const auto local =
          evaluateBarrierLocal(cell, uv, kGamma, kBeta, params);
        ASSERT_TRUE(local.valid);
        minJIdentity = std::min(minJIdentity, local.j);
        maxIdentityError =
          std::max(maxIdentityError, std::abs(local.j - Real(1)));
      }

      std::cout << "  n=" << n
                << "  min_j_identity=" << minJIdentity
                << "  jLineSearch=" << jLineSearchRatio
                << '\n';
      EXPECT_LT(maxIdentityError, Real(1e-12)) << "n=" << n;
      EXPECT_EQ(jLineSearchRatio, Real(1e-2));
    }
  }

  // The minimum FD error across a logarithmic eps sweep must reach the
  // precision a correct analytical Hessian should achieve.
  TEST(Rodin_Adaptation_Barrier, CentralFDMinimumErrorIsAtMachinePrecisionFloor)
  {
    auto mesh = makeSmallTriangleMesh();
    auto [cellCache, cellToLocal] = precomputeCellGeometry(mesh);
    const auto params = defaultParams(/*domainMeasure=*/1);

    // Several cells, several perturbation directions.
    const std::array<std::size_t, 3> seeds{ 0, 1, 2 };
    const std::array<std::size_t, 3> cellIdx{
      0, cellCache.size() / 2, cellCache.size() - 1
    };

    for (const std::size_t s : seeds)
      for (const std::size_t ci : cellIdx)
      {
        const auto& cell = cellCache[ci];
        const auto uv = deterministicUv(s, Real(1e-3));
        const auto local0 = evaluateBarrierLocal(cell, uv, kGamma, kBeta, params);
        ASSERT_TRUE(local0.valid);
        const auto p = deterministicDir(s + ci);
        const Math::Vector<Real> Hp = local0.hess * p;

        // Sweep eps from coarse to fine; minimum should be ~ 1e-9 or
        // better. Anything above 1e-5 indicates a Hessian bug.
        Real minErr = std::numeric_limits<Real>::infinity();
        for (int e = 3; e <= 9; ++e)
        {
          const Real eps = std::pow(Real(10), -Real(e));
          const auto fd = centralFD(cell, uv, params, p, eps);
          minErr = std::min(minErr, relativeError(fd, Hp));
        }
        EXPECT_LT(minErr, Real(1e-7))
            << "seed=" << s << " cell=" << ci
            << " minErr=" << minErr;
      }
  }

  // The error must decrease quadratically as eps shrinks in the
  // truncation-dominated regime (1e-3 down to ~1e-5). This is the
  // discriminating test: a wrong Hessian gives errors that bottom out
  // at O(1) regardless of eps.
  TEST(Rodin_Adaptation_Barrier, CentralFDErrorDecreasesQuadraticallyInTruncationRegime)
  {
    auto mesh = makeSmallTriangleMesh();
    auto [cellCache, cellToLocal] = precomputeCellGeometry(mesh);
    const auto params = defaultParams(1);

    const auto& cell = cellCache[cellCache.size() / 2];
    const auto uv = deterministicUv(/*seed=*/3, Real(1e-3));
    const auto local0 = evaluateBarrierLocal(cell, uv, kGamma, kBeta, params);
    ASSERT_TRUE(local0.valid);
    const auto p = deterministicDir(7);
    const Math::Vector<Real> Hp = local0.hess * p;

    // Three eps in the truncation regime. Between successive eps the
    // ratio should be ~ (eps_k / eps_{k-1})^2 = 1e-2.
    const std::array<Real, 3> epsList{ Real(1e-3), Real(1e-4), Real(1e-5) };
    std::array<Real, 3> errs;
    for (std::size_t i = 0; i < epsList.size(); ++i)
      errs[i] = relativeError(centralFD(cell, uv, params, p, epsList[i]), Hp);

    // err(1e-4) / err(1e-3) ~ 1e-2 ; err(1e-5) / err(1e-4) ~ 1e-2.
    // Be generous: accept any decrease at least as fast as 1/eps (so
    // factor better than 0.2 per eps step would be linear; we want
    // closer to 0.01 for quadratic).
    EXPECT_LT(errs[1], errs[0] * Real(0.05))
      << "expected ~1e-2 decrease, got " << errs[1] / errs[0]
      << " (err[1e-3]=" << errs[0] << ", err[1e-4]=" << errs[1] << ")";
    EXPECT_LT(errs[2], errs[1] * Real(0.05))
      << "expected ~1e-2 decrease, got " << errs[2] / errs[1]
      << " (err[1e-4]=" << errs[1] << ", err[1e-5]=" << errs[2] << ")";
  }

  // The error must increase as 1/eps in the roundoff-dominated regime
  // (1e-8, 1e-9). Together with the test above this proves the V-shape.
  TEST(Rodin_Adaptation_Barrier, CentralFDErrorIncreasesInRoundoffRegime)
  {
    auto mesh = makeSmallTriangleMesh();
    auto [cellCache, cellToLocal] = precomputeCellGeometry(mesh);
    const auto params = defaultParams(1);

    const auto& cell = cellCache[cellCache.size() / 2];
    const auto uv = deterministicUv(/*seed=*/4, Real(1e-3));
    const auto local0 = evaluateBarrierLocal(cell, uv, kGamma, kBeta, params);
    ASSERT_TRUE(local0.valid);
    const auto p = deterministicDir(11);
    const Math::Vector<Real> Hp = local0.hess * p;

    const Real e7 = relativeError(centralFD(cell, uv, params, p, Real(1e-7)), Hp);
    const Real e8 = relativeError(centralFD(cell, uv, params, p, Real(1e-8)), Hp);
    const Real e9 = relativeError(centralFD(cell, uv, params, p, Real(1e-9)), Hp);
    // We expect e9 > e8 and e8 >= e7 (within a factor of a few).
    EXPECT_GT(e9, e8 * Real(2))
      << "expected roundoff growth, got e8=" << e8 << " e9=" << e9;
    // e7 may already be in or near the optimum, so only require it
    // smaller than e9.
    EXPECT_LT(e7, e9)
      << "expected e[1e-7] < e[1e-9], got "
      << "e7=" << e7 << ", e9=" << e9;
  }
}
