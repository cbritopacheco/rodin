/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file YeohLawTest.cpp
 * @brief Verification test for the compressible Yeoh hyperelastic law.
 *
 * Mesh-free unit test (KinematicState + ConstitutivePoint built directly
 * from a displacement gradient).  Three checks, in 2D and 3D, at a
 * small-strain and a finite-strain state:
 *
 *   1. STRESS consistency:   P == dW/dF        (central finite differences)
 *   2. TANGENT consistency:  dP[delta F]       (central finite differences
 *                                               of P in random directions)
 *   3. DEGENERATION oracle:  Yeoh(c1, 0, 0, k) has the SAME energy density
 *      as MooneyRivlin(c1, 0, k), so W, P and dP must agree to machine
 *      precision.
 *
 * Exit code 0 = all PASS; 1 = some FAIL.  Run:  ./RodinSolidYeohLawTest
 */

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>

#include <Rodin/Types.h>
#include <Rodin/Math/SpatialMatrix.h>
#include <Rodin/Solid/Kinematics/KinematicState.h>
#include <Rodin/Solid/Local/ConstitutivePoint.h>
#include <Rodin/Solid/Constitutive/Yeoh.h>
#include <Rodin/Solid/Constitutive/MooneyRivlin.h>

using namespace Rodin;
using namespace Rodin::Solid;

namespace
{
  bool g_allPass = true;

  void report(const std::string& name, Real err, Real tol)
  {
    const bool pass = std::isfinite(err) && err < tol;
    g_allPass = g_allPass && pass;
    std::cout << "  " << (pass ? "[PASS] " : "[FAIL] ") << name
              << "  max rel. error = " << std::scientific
              << std::setprecision(3) << err << "  (tol " << tol << ")\n";
  }

  // Random displacement gradient with guaranteed J > 0.
  Math::SpatialMatrix<Real> randomH(size_t d, Real scale, std::mt19937& gen)
  {
    std::uniform_real_distribution<Real> dist(-1.0, 1.0);
    Math::SpatialMatrix<Real> H(d, d);
    while (true)
    {
      for (size_t i = 0; i < d; ++i)
        for (size_t j = 0; j < d; ++j)
          H(i, j) = scale * dist(gen);
      Math::SpatialMatrix<Real> F = H;
      for (size_t i = 0; i < d; ++i)
        F(i, i) += 1.0;
      if (F.determinant() > 0.25)
        return H;
    }
  }

  template <class Law>
  Real energyAt(const Law& law, size_t d, const Math::SpatialMatrix<Real>& H)
  {
    KinematicState state(d);
    state.setDisplacementGradient(H);
    ConstitutivePoint cp(state);
    typename Law::Cache cache;
    law.setCache(cache, cp);
    return law.getStrainEnergyDensity(cache, cp);
  }

  template <class Law>
  Math::SpatialMatrix<Real> stressAt(const Law& law, size_t d,
                                     const Math::SpatialMatrix<Real>& H)
  {
    KinematicState state(d);
    state.setDisplacementGradient(H);
    ConstitutivePoint cp(state);
    typename Law::Cache cache;
    law.setCache(cache, cp);
    Math::SpatialMatrix<Real> P(d, d);
    law.getFirstPiolaKirchhoffStress(P, cache, cp);
    return P;
  }

  template <class Law>
  Math::SpatialMatrix<Real> tangentAt(const Law& law, size_t d,
                                      const Math::SpatialMatrix<Real>& H,
                                      const Math::SpatialMatrix<Real>& dF)
  {
    KinematicState state(d);
    state.setDisplacementGradient(H);
    ConstitutivePoint cp(state);
    typename Law::Cache cache;
    law.setCache(cache, cp);
    Math::SpatialMatrix<Real> dP(d, d);
    law.getMaterialTangent(dP, cache, cp, dF);
    return dP;
  }

  // 1. P vs central finite differences of W (component by component).
  template <class Law>
  Real stressFDError(const Law& law, size_t d,
                     const Math::SpatialMatrix<Real>& H, Real eps)
  {
    const Math::SpatialMatrix<Real> P = stressAt(law, d, H);
    Real maxErr = 0.0;
    for (size_t i = 0; i < d; ++i)
    {
      for (size_t j = 0; j < d; ++j)
      {
        Math::SpatialMatrix<Real> Hp = H, Hm = H;
        Hp(i, j) += eps;
        Hm(i, j) -= eps;
        const Real fd =
          (energyAt(law, d, Hp) - energyAt(law, d, Hm)) / (2.0 * eps);
        const Real err =
          std::abs(fd - P(i, j)) / (1.0 + std::abs(P(i, j)));
        maxErr = std::max(maxErr, err);
      }
    }
    return maxErr;
  }

  // 2. dP[dF] vs central finite differences of P along dF.
  template <class Law>
  Real tangentFDError(const Law& law, size_t d,
                      const Math::SpatialMatrix<Real>& H,
                      const Math::SpatialMatrix<Real>& dF, Real eps)
  {
    const Math::SpatialMatrix<Real> dP = tangentAt(law, d, H, dF);
    const Math::SpatialMatrix<Real> Pp = stressAt(
        law, d, Math::SpatialMatrix<Real>(H + eps * dF));
    const Math::SpatialMatrix<Real> Pm = stressAt(
        law, d, Math::SpatialMatrix<Real>(H + (-eps) * dF));
    Real maxErr = 0.0;
    for (size_t i = 0; i < d; ++i)
    {
      for (size_t j = 0; j < d; ++j)
      {
        const Real fd = (Pp(i, j) - Pm(i, j)) / (2.0 * eps);
        const Real err =
          std::abs(fd - dP(i, j)) / (1.0 + std::abs(dP(i, j)));
        maxErr = std::max(maxErr, err);
      }
    }
    return maxErr;
  }

  Real matRelDiff(const Math::SpatialMatrix<Real>& A,
                  const Math::SpatialMatrix<Real>& B)
  {
    Real maxErr = 0.0;
    for (int i = 0; i < A.rows(); ++i)
      for (int j = 0; j < A.cols(); ++j)
        maxErr = std::max(maxErr,
            std::abs(A(i, j) - B(i, j)) / (1.0 + std::abs(B(i, j))));
    return maxErr;
  }
}

int main()
{
  std::mt19937 gen(20260610); // deterministic

  // Coronary-wall-like parameters (small-strain mu = 2 c1).
  const Real c1 = 8.9e5;
  const Real c2 = 4.0e5;
  const Real c3 = 2.0e5;
  const Real kappa = 8.3e6;

  const Yeoh yeoh(c1, c2, c3, kappa);
  const Yeoh yeohReduced(c1, 0.0, 0.0, kappa);
  const MooneyRivlin mooneyReduced(c1, 0.0, kappa);

  const Real eps = 1.0e-6;
  const Real fdTol = 5.0e-5;    // central FD truncation/cancellation budget
  const Real exactTol = 1.0e-12;

  for (size_t d : {size_t(2), size_t(3)})
  {
    for (Real scale : {0.02, 0.2}) // small-strain and finite-strain states
    {
      std::cout << "== d = " << d << ", |H| scale = " << scale << " ==\n";
      const Math::SpatialMatrix<Real> H = randomH(d, scale, gen);

      // 1. Stress vs FD of the energy.
      report("P = dW/dF      (FD)", stressFDError(yeoh, d, H, eps), fdTol);

      // 2. Tangent vs FD of the stress, several random directions.
      Real tErr = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        const Math::SpatialMatrix<Real> dF = randomH(d, 1.0, gen);
        tErr = std::max(tErr, tangentFDError(yeoh, d, H, dF, eps));
      }
      report("dP[dF]         (FD)", tErr, fdTol);

      // 3. Exact degeneration: Yeoh(c1,0,0,k) == MooneyRivlin(c1,0,k).
      const Real wErr =
        std::abs(energyAt(yeohReduced, d, H) - energyAt(mooneyReduced, d, H))
        / (1.0 + std::abs(energyAt(mooneyReduced, d, H)));
      report("W  == MR(c2=0)      ", wErr, exactTol);
      report("P  == MR(c2=0)      ",
             matRelDiff(stressAt(yeohReduced, d, H),
                        stressAt(mooneyReduced, d, H)), exactTol);
      const Math::SpatialMatrix<Real> dF = randomH(d, 1.0, gen);
      report("dP == MR(c2=0)      ",
             matRelDiff(tangentAt(yeohReduced, d, H, dF),
                        tangentAt(mooneyReduced, d, H, dF)), exactTol);
    }
  }

  std::cout << (g_allPass ? "\nALL TESTS PASSED\n" : "\nSOME TESTS FAILED\n");
  return g_allPass ? EXIT_SUCCESS : EXIT_FAILURE;
}
