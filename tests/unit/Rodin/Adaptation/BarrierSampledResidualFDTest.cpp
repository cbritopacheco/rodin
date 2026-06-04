/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//
// Finite-difference verification of the sampled barrier residual:
//
//   computeBarrierSampledCellResidual[i]
//       ~=  ( E(u + eps * e_i) - E(u - eps * e_i) ) / (2 eps)
//
// where E(u) = sum_cells computeBarrierSampledCellEnergy(cell, u, ...).
//
// The relative error |FD - R.p| / max(1, |R.p|) must show the canonical
// V-shape: decreasing as eps^2 in the truncation regime (1e-3..1e-5) and
// increasing for eps < 1e-7 (roundoff regime). The minimum must be < 1e-5.
//

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include <Rodin/Adaptation.h>
#include <Rodin/Geometry.h>
#include <Rodin/Math/SpatialVector.h>
#include <Rodin/Math/Vector.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  using LocalMesh = Mesh<Context::Local>;
  using VectorFES = P1<Math::SpatialVector<Real>, LocalMesh>;
  using VectorGF  = GridFunction<VectorFES, Math::Vector<Real>>;

  LocalMesh makeUnitSquareMesh(std::size_t n)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(2, 0);
    mesh.getConnectivity().compute(2, 1);
    mesh.scale(Real(1) / static_cast<Real>(n - 1));
    return mesh;
  }

  void setDeterministicDisplacement(VectorGF& u, Real scale, std::size_t seed)
  {
    const auto& fes = u.getFiniteElementSpace();
    const auto& mesh = fes.getMesh();
    u.getData().setZero();
    for (auto vIt = mesh.getVertex(); vIt; ++vIt)
    {
      const auto& v = *vIt;
      const Index vi = v.getIndex();
      const auto& dofs = fes.getDOFs(0, vi);
      if (dofs.size() >= 2)
      {
        u.getData()(dofs[0]) = scale * std::sin(Real(1.3) * (vi + 1) + seed);
        u.getData()(dofs[1]) = scale * std::cos(Real(0.9) * (vi + 2) - seed);
      }
    }
  }

  // Total energy = sum over cells.
  Real totalEnergy(
      const LocalMesh& mesh,
      VectorGF& u,
      Real gamma,
      const BarrierParameters& params)
  {
    Real E = Real(0);
    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
      E += computeBarrierSampledCellEnergy(*cellIt, u, gamma, params);
    return E;
  }

  // Global residual assembled from cell-local residuals.
  Math::Vector<Real> assembleGlobalResidual(
      const LocalMesh& mesh,
      VectorGF& u,
      Real gamma,
      const BarrierParameters& params)
  {
    const auto& fes = u.getFiniteElementSpace();
    const std::size_t nDofs = static_cast<std::size_t>(fes.getSize());
    Math::Vector<Real> R =
      Math::Vector<Real>::Zero(static_cast<Eigen::Index>(nDofs));

    for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
    {
      const auto& cell = *cellIt;
      const auto localR =
        computeBarrierSampledCellResidual(cell, u, gamma, params);
      // getDOFs(cellDim, cellIdx) returns global dof indices for the cell
      // in the same order as the FE basis enumeration.
      const auto& cellDofs = fes.getDOFs(
          cell.getDimension(), cell.getIndex());
      for (Eigen::Index i = 0; i < localR.size(); ++i)
        R(cellDofs[static_cast<std::size_t>(i)]) += localR(i);
    }
    return R;
  }

  // Central-FD of dE/du along global dof-space direction p.
  Real centralFDAlongDir(
      const LocalMesh& mesh,
      VectorGF& u,
      Real gamma,
      const BarrierParameters& params,
      const Math::Vector<Real>& p,
      Real eps)
  {
    const Math::Vector<Real> uOrig = u.getData();
    u.getData() = uOrig + eps * p;
    const Real Ep = totalEnergy(mesh, u, gamma, params);
    u.getData() = uOrig - eps * p;
    const Real Em = totalEnergy(mesh, u, gamma, params);
    u.getData() = uOrig;
    return (Ep - Em) / (Real(2) * eps);
  }

  Real relErr(Real fd, Real analytical)
  {
    return std::abs(fd - analytical) / std::max(Real(1), std::abs(analytical));
  }

  Math::Vector<Real> deterministicDir(std::size_t n, std::size_t seed)
  {
    Math::Vector<Real> p(static_cast<Eigen::Index>(n));
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(n); ++i)
      p(i) = std::sin(Real(0.7) * (i + 1) + static_cast<Real>(seed));
    p /= p.norm();
    return p;
  }

  constexpr Real kGamma = Real(1e-1);

  BarrierParameters defaultParams()
  {
    BarrierParameters p;
    p.jMin         = Real(0.1);
    p.domainMeasure = Real(1);
    return p;
  }

  BarrierParameters activeJBarrierParams()
  {
    BarrierParameters p = defaultParams();
    p.jBarrierWeight = Real(0.2);
    p.jBarrierSafeRatio = Real(1.2);
    return p;
  }

  BarrierParameters activeVolumeTetherParams()
  {
    BarrierParameters p = defaultParams();
    p.jVolumeTetherWeight = Real(0.2);
    return p;
  }
}

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Adaptation_BarrierSampled, EnergyAtZeroIsFinite)
  {
    auto mesh = makeUnitSquareMesh(5);
    VectorFES fes(mesh, 2);
    VectorGF u(fes);
    u.getData().setZero();
    const Real E = totalEnergy(mesh, u, kGamma, defaultParams());
    EXPECT_TRUE(std::isfinite(E));
    EXPECT_GE(E, Real(0));
    EXPECT_NEAR(E, Real(0), Real(1e-14));
  }

  // V-shape: minimum FD error < 1e-5.
  TEST(Rodin_Adaptation_BarrierSampled, GlobalResidualVShapeMinError)
  {
    auto mesh = makeUnitSquareMesh(5);
    VectorFES fes(mesh, 2);
    VectorGF u(fes);
    setDeterministicDisplacement(u, Real(1e-2), /*seed=*/0);

    const auto params = defaultParams();
    const Math::Vector<Real> R = assembleGlobalResidual(mesh, u, kGamma, params);

    const std::size_t n = static_cast<std::size_t>(R.size());
    const Math::Vector<Real> dir = deterministicDir(n, /*seed=*/3);
    const Real Rp = R.dot(dir);

    Real minErr = std::numeric_limits<Real>::infinity();
    for (int e = 3; e <= 9; ++e)
    {
      const Real eps = std::pow(Real(10), -Real(e));
      const Real fd  = centralFDAlongDir(mesh, u, kGamma, params, dir, eps);
      minErr = std::min(minErr, relErr(fd, Rp));
    }
    EXPECT_LT(minErr, Real(1e-5))
        << "V-shape minimum relative error: " << minErr;
  }

  // Quadratic decrease in truncation regime.
  TEST(Rodin_Adaptation_BarrierSampled, GlobalResidualQuadraticDecrease)
  {
    auto mesh = makeUnitSquareMesh(5);
    VectorFES fes(mesh, 2);
    VectorGF u(fes);
    setDeterministicDisplacement(u, Real(1e-2), /*seed=*/1);

    const auto params = defaultParams();
    const Math::Vector<Real> R = assembleGlobalResidual(mesh, u, kGamma, params);
    const std::size_t n = static_cast<std::size_t>(R.size());
    const Math::Vector<Real> dir = deterministicDir(n, /*seed=*/5);
    const Real Rp = R.dot(dir);

    const std::array<Real, 3> epsList{ Real(1e-3), Real(1e-4), Real(1e-5) };
    std::array<Real, 3> errs;
    for (std::size_t i = 0; i < 3; ++i)
      errs[i] = relErr(
          centralFDAlongDir(mesh, u, kGamma, params, dir, epsList[i]), Rp);

    EXPECT_LT(errs[1], errs[0] * Real(0.05))
        << "1e-3->1e-4 ratio: " << errs[1] / errs[0];
    EXPECT_LT(errs[2], errs[1] * Real(0.05))
        << "1e-4->1e-5 ratio: " << errs[2] / errs[1];
  }

  // Increase in roundoff regime.
  TEST(Rodin_Adaptation_BarrierSampled, GlobalResidualRoundoffIncrease)
  {
    auto mesh = makeUnitSquareMesh(5);
    VectorFES fes(mesh, 2);
    VectorGF u(fes);
    setDeterministicDisplacement(u, Real(1e-2), /*seed=*/2);

    const auto params = defaultParams();
    const Math::Vector<Real> R = assembleGlobalResidual(mesh, u, kGamma, params);
    const std::size_t n = static_cast<std::size_t>(R.size());
    const Math::Vector<Real> dir = deterministicDir(n, /*seed=*/7);
    const Real Rp = R.dot(dir);

    const Real e8 = relErr(
        centralFDAlongDir(mesh, u, kGamma, params, dir, Real(1e-8)), Rp);
    const Real e9 = relErr(
        centralFDAlongDir(mesh, u, kGamma, params, dir, Real(1e-9)), Rp);
    EXPECT_GT(e9, e8 * Real(2))
        << "Expected roundoff increase e8=" << e8 << " e9=" << e9;
  }

  TEST(Rodin_Adaptation_BarrierSampled, JBarrierResidualVShapeMinError)
  {
    auto mesh = makeUnitSquareMesh(5);
    VectorFES fes(mesh, 2);
    VectorGF u(fes);
    setDeterministicDisplacement(u, Real(1e-2), /*seed=*/21);

    const auto params = activeJBarrierParams();
    const Math::Vector<Real> R =
      assembleGlobalResidual(mesh, u, /*gamma=*/Real(0), params);

    const std::size_t n = static_cast<std::size_t>(R.size());
    const Math::Vector<Real> dir = deterministicDir(n, /*seed=*/23);
    const Real Rp = R.dot(dir);

    Real minErr = std::numeric_limits<Real>::infinity();
    for (int e = 3; e <= 9; ++e)
    {
      const Real eps = std::pow(Real(10), -Real(e));
      const Real fd =
        centralFDAlongDir(mesh, u, Real(0), params, dir, eps);
      minErr = std::min(minErr, relErr(fd, Rp));
    }
    EXPECT_LT(minErr, Real(1e-5))
        << "J-barrier residual V-shape minimum relative error: " << minErr;
  }

  TEST(Rodin_Adaptation_BarrierSampled, VolumeTetherResidualVShapeMinError)
  {
    auto mesh = makeUnitSquareMesh(5);
    VectorFES fes(mesh, 2);
    VectorGF u(fes);
    setDeterministicDisplacement(u, Real(1e-2), /*seed=*/31);

    const auto params = activeVolumeTetherParams();
    const Math::Vector<Real> R =
      assembleGlobalResidual(mesh, u, /*gamma=*/Real(0), params);

    const std::size_t n = static_cast<std::size_t>(R.size());
    const Math::Vector<Real> dir = deterministicDir(n, /*seed=*/33);
    const Real Rp = R.dot(dir);

    Real minErr = std::numeric_limits<Real>::infinity();
    for (int e = 3; e <= 9; ++e)
    {
      const Real eps = std::pow(Real(10), -Real(e));
      const Real fd =
        centralFDAlongDir(mesh, u, Real(0), params, dir, eps);
      minErr = std::min(minErr, relErr(fd, Rp));
    }
    EXPECT_LT(minErr, Real(1e-5))
        << "Volume-tether residual V-shape minimum relative error: "
        << minErr;
  }

  // --------------------------------------------------------------------------
  // Tangent FD tests
  //
  // The tangent K satisfies: (K p)[i] ~= (R(u + eps * p) - R(u - eps * p))[i] / (2 eps)
  // for each component i. We check the V-shape on K.p.

  namespace
  {
    Math::Vector<Real> assembleGlobalResidualVec(
        const LocalMesh& mesh,
        VectorGF& u,
        Real gamma,
        const BarrierParameters& params)
    {
      return assembleGlobalResidual(mesh, u, gamma, params);
    }

    // (R(u+eps p) - R(u-eps p)) / (2 eps) — central FD of R w.r.t. u in direction p.
    Math::Vector<Real> centralFDResidual(
        const LocalMesh& mesh,
        VectorGF& u,
        Real gamma,
        const BarrierParameters& params,
        const Math::Vector<Real>& p,
        Real eps)
    {
      const Math::Vector<Real> uOrig = u.getData();
      u.getData() = uOrig + eps * p;
      const Math::Vector<Real> Rp = assembleGlobalResidualVec(mesh, u, gamma, params);
      u.getData() = uOrig - eps * p;
      const Math::Vector<Real> Rm = assembleGlobalResidualVec(mesh, u, gamma, params);
      u.getData() = uOrig;
      return (Rp - Rm) / (Real(2) * eps);
    }

    // Assemble global tangent matrix (n x n) from cell-local tangents.
    Math::Matrix<Real> assembleGlobalTangent(
        const LocalMesh& mesh,
        VectorGF& u,
        Real gamma,
        const BarrierParameters& params)
    {
      const auto& fes = u.getFiniteElementSpace();
      const std::size_t nDofs = static_cast<std::size_t>(fes.getSize());
      Math::Matrix<Real> K =
        Math::Matrix<Real>::Zero(
            static_cast<Eigen::Index>(nDofs),
            static_cast<Eigen::Index>(nDofs));

      for (auto cellIt = mesh.getCell(); cellIt; ++cellIt)
      {
        const auto& cell = *cellIt;
        const auto localK =
          computeBarrierSampledCellTangent(cell, u, gamma, params);
        const auto& cellDofs = fes.getDOFs(
            cell.getDimension(), cell.getIndex());
        const Eigen::Index nc = localK.rows();
        for (Eigen::Index i = 0; i < nc; ++i)
          for (Eigen::Index j = 0; j < nc; ++j)
            K(cellDofs[static_cast<std::size_t>(i)],
              cellDofs[static_cast<std::size_t>(j)]) += localK(i, j);
      }
      return K;
    }
  }

  TEST(Rodin_Adaptation_BarrierSampled, TangentVShapeMinError)
  {
    auto mesh = makeUnitSquareMesh(5);
    VectorFES fes(mesh, 2);
    VectorGF u(fes);
    setDeterministicDisplacement(u, Real(1e-2), /*seed=*/10);

    const auto params = defaultParams();
    const Math::Matrix<Real> K = assembleGlobalTangent(mesh, u, kGamma, params);

    const std::size_t n = static_cast<std::size_t>(K.rows());
    const Math::Vector<Real> dir = deterministicDir(n, /*seed=*/13);
    const Math::Vector<Real> Kp  = K * dir;

    Real minErr = std::numeric_limits<Real>::infinity();
    for (int e = 3; e <= 8; ++e)
    {
      const Real eps = std::pow(Real(10), -Real(e));
      const Math::Vector<Real> fd =
        centralFDResidual(mesh, u, kGamma, params, dir, eps);
      const Real err =
        (fd - Kp).norm() / std::max(Real(1), Kp.norm());
      minErr = std::min(minErr, err);
    }
    EXPECT_LT(minErr, Real(1e-5))
        << "Tangent V-shape minimum error: " << minErr;
  }

  TEST(Rodin_Adaptation_BarrierSampled, TangentQuadraticDecreaseInTruncation)
  {
    auto mesh = makeUnitSquareMesh(5);
    VectorFES fes(mesh, 2);
    VectorGF u(fes);
    setDeterministicDisplacement(u, Real(1e-2), /*seed=*/11);

    const auto params = defaultParams();
    const Math::Matrix<Real> K = assembleGlobalTangent(mesh, u, kGamma, params);
    const std::size_t n = static_cast<std::size_t>(K.rows());
    const Math::Vector<Real> dir = deterministicDir(n, /*seed=*/15);
    const Math::Vector<Real> Kp  = K * dir;

    const std::array<Real, 3> epsList{ Real(1e-3), Real(1e-4), Real(1e-5) };
    std::array<Real, 3> errs;
    for (std::size_t i = 0; i < 3; ++i)
    {
      const Math::Vector<Real> fd =
        centralFDResidual(mesh, u, kGamma, params, dir, epsList[i]);
      errs[i] = (fd - Kp).norm() / std::max(Real(1), Kp.norm());
    }

    EXPECT_LT(errs[1], errs[0] * Real(0.05))
        << "1e-3->1e-4 ratio: " << errs[1] / errs[0];
    EXPECT_LT(errs[2], errs[1] * Real(0.05))
        << "1e-4->1e-5 ratio: " << errs[2] / errs[1];
  }

  TEST(Rodin_Adaptation_BarrierSampled, JBarrierTangentVShapeMinError)
  {
    auto mesh = makeUnitSquareMesh(5);
    VectorFES fes(mesh, 2);
    VectorGF u(fes);
    setDeterministicDisplacement(u, Real(1e-2), /*seed=*/25);

    const auto params = activeJBarrierParams();
    const Math::Matrix<Real> K =
      assembleGlobalTangent(mesh, u, /*gamma=*/Real(0), params);

    const std::size_t n = static_cast<std::size_t>(K.rows());
    const Math::Vector<Real> dir = deterministicDir(n, /*seed=*/27);
    const Math::Vector<Real> Kp  = K * dir;

    Real minErr = std::numeric_limits<Real>::infinity();
    for (int e = 3; e <= 8; ++e)
    {
      const Real eps = std::pow(Real(10), -Real(e));
      const Math::Vector<Real> fd =
        centralFDResidual(mesh, u, Real(0), params, dir, eps);
      const Real err =
        (fd - Kp).norm() / std::max(Real(1), Kp.norm());
      minErr = std::min(minErr, err);
    }
    EXPECT_LT(minErr, Real(1e-5))
        << "J-barrier tangent V-shape minimum error: " << minErr;
  }

  TEST(Rodin_Adaptation_BarrierSampled, VolumeTetherTangentVShapeMinError)
  {
    auto mesh = makeUnitSquareMesh(5);
    VectorFES fes(mesh, 2);
    VectorGF u(fes);
    setDeterministicDisplacement(u, Real(1e-2), /*seed=*/35);

    const auto params = activeVolumeTetherParams();
    const Math::Matrix<Real> K =
      assembleGlobalTangent(mesh, u, /*gamma=*/Real(0), params);

    const std::size_t n = static_cast<std::size_t>(K.rows());
    const Math::Vector<Real> dir = deterministicDir(n, /*seed=*/37);
    const Math::Vector<Real> Kp  = K * dir;

    Real minErr = std::numeric_limits<Real>::infinity();
    for (int e = 3; e <= 8; ++e)
    {
      const Real eps = std::pow(Real(10), -Real(e));
      const Math::Vector<Real> fd =
        centralFDResidual(mesh, u, Real(0), params, dir, eps);
      const Real err =
        (fd - Kp).norm() / std::max(Real(1), Kp.norm());
      minErr = std::min(minErr, err);
    }
    EXPECT_LT(minErr, Real(1e-5))
        << "Volume-tether tangent V-shape minimum error: " << minErr;
  }

}
