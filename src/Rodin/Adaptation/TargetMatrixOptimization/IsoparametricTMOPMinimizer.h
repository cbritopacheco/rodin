/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_ISOPARAMETRICTMOPMINIMIZER_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_ISOPARAMETRICTMOPMINIMIZER_H

/**
 * @file
 * @brief Residual-only minimization driver for isoparametric TMOP energies.
 *
 * This production driver minimizes a displacement energy E(u) with:
 * - energy evaluation,
 * - first variation / residual assembly R(u; v) = DE(u)[v],
 * - a fixed elliptic preconditioner P = M + ell^2 K,
 * - determinant-safe Armijo backtracking.
 *
 * It deliberately does not assemble or use the TMOP tangent/Hessian.
 */

#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>

#include "Rodin/Assembly/Default.h"
#include "Rodin/Geometry/ParametricTransformation.h"
#include "Rodin/Math.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Types.h"
#include "Rodin/Variational.h"

#include "IsoparametricGeometry.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  struct IsoparametricTMOPMinimizerParameters
  {
    Index maxIterations = 50;
    Real gradientTolerance = Real(1e-10);
    Real stepTolerance = Real(1e-12);
    Real energyTolerance = Real(1e-14);
    Real armijo = Real(1e-4);
    Real alphaMin = Real(1e-8);
    Real alphaShrink = Real(0.5);
    Real preconditionerLength = Real(0.05);
    Real minDetFloor = Real(0);
    bool printIterations = false;
    std::function<void(Real)> acceptedEnergyMonitor;
  };

  struct IsoparametricTMOPMinimizerReport
  {
    Index iterations = 0;
    Index acceptedSteps = 0;
    Index rejectedSteps = 0;
    Index energyEvaluations = 0;
    bool converged = false;
    bool acceptedAnyStep = false;
    Real initialEnergy = 0;
    Real finalEnergy = 0;
    Real initialGradientNorm = 0;
    Real finalGradientNorm = 0;
    Real finalStepNorm = 0;
    Real finalAlpha = 0;
    Real finalDirectionalDerivative = 0;
  };

  template <class Element, class FES, class Data, class Target>
  bool isTargetAdmissible(
      const Geometry::LocalMesh& mesh,
      const Variational::GridFunction<FES, Data>& u,
      const Element& fe,
      const Target& target,
      Real detFloor,
      Geometry::Polytope::Type geometry = Geometry::Polytope::Type::Triangle,
      Index quadratureOrder = 4)
  {
    const auto& qf = QF::PolytopeQuadratureFormula::get(
        quadratureOrder, geometry);
    const auto& ref = Element::getNodes(geometry);
    const auto& conn = mesh.getConnectivity();
    const Index nc = static_cast<Index>(mesh.getCellCount());
    for (Index ci = 0; ci < nc; ++ci)
    {
      if (conn.getGeometry(2, ci) != geometry)
        continue;

      auto cellIt = mesh.getPolytope(2, ci);
      Geometry::PointCloud pc(2, ref.size());
      for (size_t j = 0; j < ref.size(); ++j)
      {
        const Geometry::Point p(*cellIt, ref[j]);
        const auto x = p.getPhysicalCoordinates();
        const auto dx = u.getValue(p);
        pc(0, j) = x[0] + dx[0];
        pc(1, j) = x[1] + dx[1];
      }

      Geometry::ParametricTransformation<Element> trial(pc, fe);
      for (size_t q = 0; q < qf.getSize(); ++q)
      {
        const auto& rc = qf.getPoint(q);
        Math::SpatialMatrix<Real> A;
        trial.jacobian(A, rc);
        if (A.rows() != 2 || A.cols() != 2)
          return false;

        const auto W = target.evaluate(*cellIt, rc);
        if (W.rows() != 2 || W.cols() != 2)
          return false;

        const Real detA = A.determinant();
        const Real detW = W.determinant();
        if (!std::isfinite(detA) || !std::isfinite(detW)
            || !(detA > detFloor) || !(detW > detFloor))
          return false;

        const Real detT = detA / detW;
        if (!std::isfinite(detT) || !(detT > detFloor))
          return false;
      }
    }
    return true;
  }

  template <class Element,
            class FES, class Data,
            class TrialF, class TestF,
            class MakeResidual, class Energy, class IsAdmissible>
  IsoparametricTMOPMinimizerReport minimizeIsoparametricTMOP(
      Geometry::LocalMesh& mesh,
      const Element& fe,
      Variational::GridFunction<FES, Data>& u,
      TrialF& p,
      TestF& v,
      MakeResidual makeResidual,
      Energy energy,
      IsAdmissible isAdmissible,
      Geometry::Attribute interfaceAttribute,
      const IsoparametricTMOPMinimizerParameters& params = {})
  {
    using VectorType = Data;

    IsoparametricTMOPMinimizerReport rep;
    rep.initialEnergy = energy();
    ++rep.energyEvaluations;
    rep.finalEnergy = rep.initialEnergy;

    const Real ell = std::max(params.preconditionerLength, Real(0));
    Variational::BilinearForm P(p, v);
    P =
        Variational::Integral(Variational::Dot(p, v))
      + Variational::Integral((ell * ell) * Variational::Jacobian(p),
          Variational::Jacobian(v));
    P.assemble();
    Eigen::SparseLU<std::decay_t<decltype(P.getOperator())>> lu;
    lu.compute(P.getOperator());
    if (lu.info() != Eigen::Success)
      return rep;

    for (Index it = 0; it < params.maxIterations; ++it)
    {
      Variational::LinearForm R(v);
      R = makeResidual();
      R.assemble();
      const VectorType g = R.getVector();
      if (!g.allFinite())
        break;

      const Real gNorm = g.norm();
      if (it == 0)
        rep.initialGradientNorm = gNorm;
      rep.finalGradientNorm = gNorm;
      if (gNorm <= params.gradientTolerance)
      {
        rep.converged = true;
        break;
      }

      VectorType direction = lu.solve(-g);
      if (lu.info() != Eigen::Success)
        break;
      if (!direction.allFinite())
        break;

      Real directionalDerivative = g.dot(direction);
      if (!(directionalDerivative < Real(0)))
      {
        direction *= Real(-1);
        directionalDerivative = g.dot(direction);
      }
      if (!(directionalDerivative < Real(0)))
      {
        direction = -g;
        directionalDerivative = -g.squaredNorm();
      }
      if (!(directionalDerivative < Real(0)))
        break;

      const VectorType u0 = u.getData();
      const Real e0 = rep.finalEnergy;
      Real alpha = Real(1);
      bool accepted = false;
      Real acceptedEnergy = e0;
      VectorType acceptedU = u0;

      while (alpha >= params.alphaMin)
      {
        u.getData() = u0 + alpha * direction;
        const bool finiteU = u.getData().allFinite();
        const bool admissible = finiteU && isAdmissible();
        const Real eTrial = admissible
          ? energy()
          : std::numeric_limits<Real>::infinity();
        ++rep.energyEvaluations;

        const bool sufficientDecrease =
          std::isfinite(eTrial)
          && eTrial <= e0 + params.armijo * alpha * directionalDerivative;
        if (admissible && sufficientDecrease)
        {
          accepted = true;
          acceptedEnergy = eTrial;
          acceptedU = u.getData();
          break;
        }
        alpha *= params.alphaShrink;
        ++rep.rejectedSteps;
      }

      if (!accepted)
      {
        u.getData() = u0;
        break;
      }

      u.getData() = acceptedU;
      rep.acceptedAnyStep = true;
      ++rep.acceptedSteps;
      ++rep.iterations;
      rep.finalAlpha = alpha;
      rep.finalDirectionalDerivative = directionalDerivative;
      rep.finalStepNorm = (alpha * direction).norm();
      rep.finalEnergy = acceptedEnergy;
      if (params.acceptedEnergyMonitor)
        params.acceptedEnergyMonitor(rep.finalEnergy);

      if (params.printIterations)
      {
        std::cout
          << "tmop_minimize it=" << it
          << " energy=" << rep.finalEnergy
          << " gradient=" << gNorm
          << " alpha=" << alpha
          << " gTp=" << directionalDerivative
          << " step_norm=" << rep.finalStepNorm
          << "\n";
      }

      if (rep.finalStepNorm <= params.stepTolerance)
      {
        rep.converged = true;
        break;
      }
      if (std::abs(e0 - rep.finalEnergy)
          <= params.energyTolerance * std::max(Real(1), std::abs(e0)))
      {
        rep.converged = true;
        break;
      }
    }

    moveMesh(mesh, u, fe, interfaceAttribute);
    return rep;
  }
}

#endif
