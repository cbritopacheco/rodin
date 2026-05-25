/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_ISOPARAMETRICTMOPSOLVER_H
#define RODIN_ADAPTATION_TARGETMATRIXOPTIMIZATION_ISOPARAMETRICTMOPSOLVER_H

/**
 * @file
 * @brief Free-function Newton driver for the isoparametric TMOP problem.
 *
 * Composes with Rodin's existing variational form language:
 *
 *   auto makeResidual = [&] { return quality.residual(u, v) + ... ; };
 *   auto makeTangent  = [&] { return quality.tangent(u, du, v) + ... ; };
 *   auto energy       = [&] { return quality.energy(u) + ... ; };
 *
 *   solveIsoparametricTMOP<Element>(
 *       mesh, fe, u, du, v,
 *       makeResidual, makeTangent, energy,
 *       interfaceAttribute);
 *
 * The driver does full Newton steps plus a curved-validity gate; it
 * does NOT enforce a hard interface manifold constraint. The interface
 * fit term in the residual is responsible for pulling interface vertices
 * onto phi=0 *smoothly* -- Newton stays smooth, and there is no line
 * search branch with alpha-dependent interface updates.
 *
 * If you need exact-on-manifold interface vertices, add a smooth penalty
 * term to the composition (e.g. a vertex-based phi^2 weight), do NOT
 * overwrite u inside a step-rejection branch.
 *
 * On success the driver commits the displacement to the mesh via
 * moveMesh(mesh, u). No state, no class.
 */

#include <Eigen/SparseLU>

#include <cmath>
#include <iostream>
#include <limits>
#include <utility>

#include "Rodin/Assembly/Default.h"
#include "Rodin/Math.h"
#include "Rodin/Types.h"
#include "Rodin/Variational.h"

#include "IsoparametricGeometry.h"

namespace Rodin::Adaptation::TargetMatrixOptimization
{
  /// Outcome of one isoparametric TMOP solve.
  struct IsoparametricTMOPSolverReport
  {
    Index iterations = 0;        ///< Newton steps accepted.
    bool converged = false;      ///< Residual fell below the tolerance.
    Real initialResidual = 0;
    Real finalResidual = 0;
    Real lastQuadraticRate = -1;
    Real lastAcceptedAlpha = 0;
    bool acceptedStep = false;
    Index rejectedSteps = 0;
  };

  /// User-tunable Newton parameters.
  struct IsoparametricTMOPParameters
  {
    Index maxIterations = 8;
    Real residualTolerance = Real(1e-10);
    /// Step-norm cutoff to declare convergence.
    Real stepTolerance = Real(1e-10);
    /// Strict curved-det floor (absolute).
    Real minDetFloor = Real(0);
    /// Print per-iteration Newton diagnostics.
    bool printIterations = false;
  };

  /// Full Newton plus curved-validity gate, for the isoparametric TMOP problem.
  ///
  /// The user composes the residual / tangent / energy with Rodin's existing
  /// form language INSIDE the callables `makeResidual` / `makeTangent` /
  /// `energy`. The driver invokes them once per Newton iteration so that
  /// `LinearForm` / `BilinearForm` are rebound to a fresh expression each
  /// time -- the Rodin-idiomatic assembly pattern.
  ///
  /// @tparam Element  Lagrange-nodal mesh element used to represent the
  ///                  parametric geometry (e.g. RealH1Element<2>).
  ///
  /// @param mesh                 Parametric mesh (already upgraded).
  /// @param fe                   The mesh element instance.
  /// @param u                    Displacement GridFunction.
  /// @param v                    Test function.
  /// @param makeResidual()       Returns the residual expression bound to
  ///                             the current `u`.
  /// @param makeTangent()        Returns the tangent expression bound to
  ///                             the current `u`.
  /// @param energy()             Returns the merit energy at the current u.
  /// @param interfaceAttribute   Attribute identifying interface edges
  ///                             (used only to refresh curved edge
  ///                             transformations in the final moveMesh).
  /// @param params               Solver parameters.
  ///
  /// Smooth method: the interface fit is the level-set penalty term that
  /// the caller put in `makeResidual` / `makeTangent` (Knupp, Kolev, Mittal,
  /// Tomov 2021). No geometric surface projection, no hard manifold
  /// constraint, and no step damping: a rejected full step is reported as a
  /// failed Newton step immediately.
  /// Tangential relaxation is provided automatically by the quality metric.
  template <class Element,
            class FES, class Data,
            class TrialF, class TestF,
            class MakeResidual, class MakeTangent,
            class Energy>
  IsoparametricTMOPSolverReport solveIsoparametricTMOP(
      Geometry::LocalMesh& mesh,
      const Element& fe,
      Variational::GridFunction<FES, Data>& u,
      TrialF& du,
      TestF& v,
      MakeResidual makeResidual,
      MakeTangent makeTangent,
      Energy energy,
      Geometry::Attribute interfaceAttribute,
      const IsoparametricTMOPParameters& params = {})
  {
    using GeoData = Data;

    IsoparametricTMOPSolverReport rep;
    Real previousResidual = -1;
    for (Index it = 0; it < params.maxIterations; ++it)
    {
      Variational::LinearForm R(v);
      R = makeResidual();
      R.assemble();
      const auto r = R.getVector();
      if (!r.allFinite()) break;
      const Real rNorm = r.norm();
      const Real quadraticRate =
        (previousResidual > Real(0) && std::isfinite(previousResidual))
        ? (rNorm / (previousResidual * previousResidual))
        : Real(-1);
      rep.lastQuadraticRate = quadraticRate;
      if (it == 0) rep.initialResidual = rNorm;
      rep.finalResidual = rNorm;
      if (rNorm <= params.residualTolerance)
      { rep.converged = true; break; }

      Variational::BilinearForm J(du, v);
      J = makeTangent();
      J.assemble();
      Eigen::SparseLU<std::decay_t<decltype(J.getOperator())>> lu;
      lu.compute(J.getOperator());
      if (lu.info() != Eigen::Success) break;
      const GeoData dx = lu.solve(-r);
      if (lu.info() != Eigen::Success || !dx.allFinite()) break;

      const GeoData u0 = u.getData();
      const Real eBefore = energy();
      u.getData() = u0 + dx;
      const Real e = energy();
      const bool accepted =
        std::isfinite(e)
        && isCurvedMoveValid(mesh, u, fe, params.minDetFloor);
      rep.lastAcceptedAlpha = Real(1);
      rep.acceptedStep = accepted;
      if (params.printIterations)
      {
        std::cout
          << "tmop_newton it=" << it
          << " residual=" << rNorm
          << " qrate=" << quadraticRate
          << " step_norm=" << dx.norm()
          << " alpha=1"
          << " energy_before=" << eBefore
          << " energy_after=" << e
          << " full_step_valid=" << (accepted ? 1 : 0)
          << "\n";
      }
      if (!accepted)
      {
        u.getData() = u0;
        ++rep.rejectedSteps;
        rep.acceptedStep = false;
        break;
      }
      previousResidual = rNorm;
      ++rep.iterations;
      if (dx.norm() <= params.stepTolerance) break;
    }

    // Commit: move the parametric mesh by the solved displacement.
    moveMesh(mesh, u, fe, interfaceAttribute);
    return rep;
  }
}

#endif
