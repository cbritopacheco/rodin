/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/* Kelvin-ball resistance metrics on a chamber or on the reconstructed domain. */
#ifndef KELVIN_BALL_METRICS_H
#define KELVIN_BALL_METRICS_H

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include <Rodin/Assembly.h>
#include <Rodin/Alert/Info.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include "Common.h"
#include "RotatedNitscheIntegrator.h"
#include "SewedOutput.h"

namespace KelvinBall
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  struct Parameters
  {
      Real h;
      Real nitschePenalty = DefaultNitschePenalty;
      Real stabilizationFactor = DefaultStabilizationFactor;

      Real pressureStabilization() const
      {
        return stabilizationFactor * h * h / Mu;
      }
  };

  struct Values
  {
      Real k = 0;
      Real c = 0;
      Real q = 0;
      Real rho = 0;
      Real nitscheJump = 0;
      Real couplingSymmetry = 0;
  };

  /**
   * @brief Evaluates the resistance tensor and Kelvin-ball objective.
   *
   * The object owns the stabilized Stokes parameters and the complete
   * evaluation procedure.  Its private stages assemble the state families,
   * verify the rotated traces, and contract the resulting resistance tensor.
   */
  class Metrics
  {
    public:
      explicit Metrics(const Parameters& parameters)
        : m_parameters(parameters)
      {}

      Values evaluateChamber(Mesh& fluid) const;

      Values evaluateSewed(Mesh& fluid) const;

    private:
      void checkPairs(const Mesh& mesh, const RotatedNitscheIntegrator& coupling) const
      {
        const auto& locator = coupling.getLocator();
        for (const RotationPair& pair : RotationPairs)
        {
          size_t found = 0;
          size_t missed = 0;
          for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
          {
            if (face->getAttribute() != pair.slave)
              continue;
            const auto& qf = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
            const auto& quadrature = face->getQuadrature(qf);
            for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
            {
              if (locator.locate(
                    pair.master, pair.rotation * quadrature.getPoint(qp).vector()))
                ++found;
              else
                ++missed;
            }
          }
          if (found == 0 || missed != 0)
            throw std::runtime_error(
              "The rotated chamber faces do not cover each other.");
        }
      }

      template <class Space, class VelocityA, class VelocityB>
      Real viscousProduct(const Space& space, const VelocityA& a, const VelocityB& b,
        Real multiplicity) const
      {
        TrialFunction u(space);
        TestFunction v(space);
        const auto Du = Real(0.5) * (Jacobian(u) + Jacobian(u).T());
        const auto Dv = Real(0.5) * (Jacobian(v) + Jacobian(v).T());
        BilinearForm form(u, v);
        form = Integral(Real(2) * Mu * Du, Dv);
        form.assemble();
        return multiplicity * form(a, b);
      }

      template <class Rigid0, class Rigid1, class Rigid2, class Velocity0,
        class Velocity1, class Velocity2, class Pressure0, class Pressure1,
        class Pressure2>
      Real solveChamberFamily(const VelocitySpace& Vh, const PressureSpace& Qh,
        const RotatedNitscheIntegrator& coupling, const Rigid0& rigid0,
        const Rigid1& rigid1, const Rigid2& rigid2, Velocity0& velocity0,
        Velocity1& velocity1, Velocity2& velocity2, Pressure0& pressure0,
        Pressure1& pressure1, Pressure2& pressure2) const
      {
        TrialFunction u0(Vh), u1(Vh), u2(Vh);
        TrialFunction p0(Qh), p1(Qh), p2(Qh);
        TestFunction v0(Vh), v1(Vh), v2(Vh);
        TestFunction q0(Qh), q1(Qh), q2(Qh);
        const auto Du0 = Real(0.5) * (Jacobian(u0) + Jacobian(u0).T());
        const auto Du1 = Real(0.5) * (Jacobian(u1) + Jacobian(u1).T());
        const auto Du2 = Real(0.5) * (Jacobian(u2) + Jacobian(u2).T());
        const auto Dv0 = Real(0.5) * (Jacobian(v0) + Jacobian(v0).T());
        const auto Dv1 = Real(0.5) * (Jacobian(v1) + Jacobian(v1).T());
        const auto Dv2 = Real(0.5) * (Jacobian(v2) + Jacobian(v2).T());
        const Real stabilization = m_parameters.pressureStabilization();

        Problem stokes(u0, p0, u1, p1, u2, p2, v0, q0, v1, q1, v2, q2);
        stokes = Integral(Real(2) * Mu * Du0, Dv0) - Integral(p0, Div(v0)) -
          Integral(Div(u0), q0) - stabilization * Integral(Grad(p0), Grad(q0)) +
          Integral(Real(2) * Mu * Du1, Dv1) - Integral(p1, Div(v1)) -
          Integral(Div(u1), q1) - stabilization * Integral(Grad(p1), Grad(q1)) +
          Integral(Real(2) * Mu * Du2, Dv2) - Integral(p2, Div(v2)) -
          Integral(Div(u2), q2) - stabilization * Integral(Grad(p2), Grad(q2)) +
          DirichletBC(u0, rigid0).on(Gamma) + DirichletBC(u1, rigid1).on(Gamma) +
          DirichletBC(u2, rigid2).on(Gamma) +
          DirichletBC(u0, VectorFunction{0, 0, 0}).on(Outer) +
          DirichletBC(u1, VectorFunction{0, 0, 0}).on(Outer) +
          DirichletBC(u2, VectorFunction{0, 0, 0}).on(Outer);
        stokes.assemble();
        auto& system = stokes.getLinearSystem();
        coupling.assembleStokes<1>(Vh, Qh, stokes.getTrialOffsets(), system, Mu,
          m_parameters.nitschePenalty, -stabilization, FlatSet<Attribute>{Gamma, Outer});
        solveDirect(stokes);
        const Real residual =
          (system.getOperator() * system.getSolution() - system.getVector()).norm() /
          std::max(system.getVector().norm(), Real(1));
        if (!std::isfinite(residual) || residual > LinearResidualTolerance)
          throw std::runtime_error("The chamber Stokes system did not converge.");
        velocity0 = u0.getSolution();
        velocity1 = u1.getSolution();
        velocity2 = u2.getSolution();
        pressure0 = p0.getSolution();
        pressure1 = p1.getSolution();
        pressure2 = p2.getSolution();
        return coupling.familyJump(velocity0, velocity1, velocity2);
      }

      template <class Rigid0, class Rigid1, class Rigid2, class Velocity0,
        class Velocity1, class Velocity2, class Pressure0, class Pressure1,
        class Pressure2>
      Real solveFullFamily(const VelocitySpace& Vh, const PressureSpace& Qh,
        const Rigid0& rigid0, const Rigid1& rigid1, const Rigid2& rigid2,
        Velocity0& velocity0, Velocity1& velocity1, Velocity2& velocity2,
        Pressure0& pressure0, Pressure1& pressure1, Pressure2& pressure2) const
      {
        TrialFunction u0(Vh), u1(Vh), u2(Vh);
        TrialFunction p0(Qh), p1(Qh), p2(Qh);
        TestFunction v0(Vh), v1(Vh), v2(Vh);
        TestFunction q0(Qh), q1(Qh), q2(Qh);
        const auto Du0 = Real(0.5) * (Jacobian(u0) + Jacobian(u0).T());
        const auto Du1 = Real(0.5) * (Jacobian(u1) + Jacobian(u1).T());
        const auto Du2 = Real(0.5) * (Jacobian(u2) + Jacobian(u2).T());
        const auto Dv0 = Real(0.5) * (Jacobian(v0) + Jacobian(v0).T());
        const auto Dv1 = Real(0.5) * (Jacobian(v1) + Jacobian(v1).T());
        const auto Dv2 = Real(0.5) * (Jacobian(v2) + Jacobian(v2).T());
        const Real stabilization = m_parameters.pressureStabilization();
        P0g gauge(Vh.getMesh());
        TrialFunction lambda0(gauge), lambda1(gauge), lambda2(gauge);
        TestFunction eta0(gauge), eta1(gauge), eta2(gauge);
        Problem stokes(u0, p0, u1, p1, u2, p2, lambda0, lambda1, lambda2, v0, q0, v1, q1,
          v2, q2, eta0, eta1, eta2);
        stokes = Integral(Real(2) * Mu * Du0, Dv0) - Integral(p0, Div(v0)) -
          Integral(Div(u0), q0) - stabilization * Integral(Grad(p0), Grad(q0)) +
          Integral(lambda0, q0) + Integral(p0, eta0) + Integral(Real(2) * Mu * Du1, Dv1) -
          Integral(p1, Div(v1)) - Integral(Div(u1), q1) -
          stabilization * Integral(Grad(p1), Grad(q1)) + Integral(lambda1, q1) +
          Integral(p1, eta1) + Integral(Real(2) * Mu * Du2, Dv2) - Integral(p2, Div(v2)) -
          Integral(Div(u2), q2) - stabilization * Integral(Grad(p2), Grad(q2)) +
          Integral(lambda2, q2) + Integral(p2, eta2) + DirichletBC(u0, rigid0).on(Gamma) +
          DirichletBC(u1, rigid1).on(Gamma) + DirichletBC(u2, rigid2).on(Gamma) +
          DirichletBC(u0, VectorFunction{0, 0, 0}).on(Outer) +
          DirichletBC(u1, VectorFunction{0, 0, 0}).on(Outer) +
          DirichletBC(u2, VectorFunction{0, 0, 0}).on(Outer);
        stokes.assemble();
        auto& system = stokes.getLinearSystem();
        solveDirect(stokes);
        const Real residual =
          (system.getOperator() * system.getSolution() - system.getVector()).norm() /
          std::max(system.getVector().norm(), Real(1));
        if (!std::isfinite(residual) || residual > LinearResidualTolerance)
          throw std::runtime_error("The full-domain Stokes system did not converge.");
        velocity0 = u0.getSolution();
        velocity1 = u1.getSolution();
        velocity2 = u2.getSolution();
        pressure0 = p0.getSolution();
        pressure1 = p1.getSolution();
        pressure2 = p2.getSolution();
        return residual;
      }

      template <class Velocity0, class Velocity1, class Velocity2, class Rotation0,
        class Rotation1, class Rotation2>
      Values resistance(const VelocitySpace& Vh, const Velocity0& uT0,
        const Velocity1& uT1, const Velocity2& uT2, const Rotation0& uR0,
        const Rotation1& uR1, const Rotation2& uR2, Real multiplicity) const
      {
        const std::array<const Velocity0*, 3> translations{&uT0, &uT1, &uT2};
        const std::array<const Rotation0*, 3> rotations{&uR0, &uR1, &uR2};
        Math::SpatialMatrix<Real> coupling(3, 3);
        Real k = 0;
        Real q = 0;
        for (size_t i = 0; i < 3; ++i)
        {
          k += viscousProduct(Vh, *translations[i], *translations[i], multiplicity);
          q += viscousProduct(Vh, *rotations[i], *rotations[i], multiplicity);
          for (size_t j = 0; j < 3; ++j)
            coupling(i, j) = viscousProduct(Vh, *translations[i], *rotations[j], 1);
        }
        k /= 3;
        q /= 3;
        const Real c = multiplicity * coupling.trace() / 3;
        if (!(k > 0) || !(q > 0))
          throw std::runtime_error("The resistance scale is not positive.");
        Math::SpatialMatrix<Real> reconstructed(3, 3);
        reconstructed.setZero();
        for (const auto& R : SewedOutput::getCubeRotations())
          reconstructed += R * coupling * R.transpose();
        Math::SpatialMatrix<Real> isotropic(3, 3);
        isotropic.setIdentity();
        isotropic *= c;
        return {k, c, q, std::abs(c) / std::sqrt(k * q), 0,
          (reconstructed - isotropic).norm() / std::max(reconstructed.norm(), Real(1))};
      }

    public:
      template <class Velocity0, class Velocity1, class Velocity2, class Rotation0,
        class Rotation1, class Rotation2, class Pressure0, class Pressure1,
        class Pressure2, class RotationPressure0, class RotationPressure1,
        class RotationPressure2>
      Values evaluateChamber(const VelocitySpace& Vh, const PressureSpace& Qh,
        const RotatedNitscheIntegrator& coupling, Velocity0& uT0, Velocity1& uT1,
        Velocity2& uT2, Rotation0& uR0, Rotation1& uR1, Rotation2& uR2, Pressure0& pT0,
        Pressure1& pT1, Pressure2& pT2, RotationPressure0& pR0, RotationPressure1& pR1,
        RotationPressure2& pR2) const
      {
        checkPairs(Vh.getMesh(), coupling);
        Alert::Info() << Alert::Text<Alert::YellowT>(
                           Alert::Yellow, "Translational Stokes family")
                           .setBold()
                      << Alert::Raise;
        const Real translationJump = solveChamberFamily(Vh, Qh, coupling,
          VectorFunction{1, 0, 0}, VectorFunction{0, 1, 0}, VectorFunction{0, 0, 1}, uT0,
          uT1, uT2, pT0, pT1, pT2);
        Alert::Info() << Alert::Text<Alert::YellowT>(
                           Alert::Yellow, "Rotational Stokes family")
                           .setBold()
                      << Alert::Raise;
        const Real rotationJump = solveChamberFamily(Vh, Qh, coupling,
          VectorFunction{0, -F::z, F::y}, VectorFunction{F::z, 0, -F::x},
          VectorFunction{-F::y, F::x, 0}, uR0, uR1, uR2, pR0, pR1, pR2);
        Values result = resistance(
          Vh, uT0, uT1, uT2, uR0, uR1, uR2, static_cast<Real>(ChamberMultiplicity));
        result.nitscheJump = std::max(translationJump, rotationJump);
        return result;
      }

    private:
      Parameters m_parameters;
  };
}

#endif
