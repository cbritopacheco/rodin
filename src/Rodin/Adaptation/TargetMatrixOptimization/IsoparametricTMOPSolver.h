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
 * The driver does Newton + Armijo line search + curved-validity gate; it
 * does NOT enforce a hard interface manifold constraint. The interface
 * fit term in the residual is responsible for pulling interface vertices
 * onto phi=0 *smoothly* -- Newton stays smooth, the line search has no
 * alpha-invariant jump on interface dofs, and convergence is well-behaved.
 *
 * If you need exact-on-manifold interface vertices, add a smooth penalty
 * term to the composition (e.g. a vertex-based phi^2 weight), do NOT
 * overwrite u inside the line search.
 *
 * On success the driver commits the displacement to the mesh via
 * moveMesh(mesh, u). No state, no class.
 */

#include <Eigen/SparseLU>

#include <cmath>
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
  };

  /// User-tunable Newton parameters.
  struct IsoparametricTMOPParameters
  {
    Index maxIterations = 8;
    Index lineSearchBudget = 30;
    Real residualTolerance = Real(1e-10);
    /// Step-norm cutoff to declare convergence.
    Real stepTolerance = Real(1e-10);
    /// Strict curved-det floor (absolute).
    Real minDetFloor = Real(0);
    /// Relative slack on the merit non-increase line-search gate. With the
    /// smooth (projection-free) method a near-monotone tolerance is enough.
    Real meritRelativeTolerance = Real(1e-12);
  };

  /// Newton + Armijo line search + interface tangential constraint +
  /// curved-validity gate, for the isoparametric TMOP problem.
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
  /// constraint -- the line search is a plain Armijo backtracking with a
  /// curved-validity gate, so Newton stays smooth. Tangential relaxation is
  /// provided automatically by the quality metric.
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
    for (Index it = 0; it < params.maxIterations; ++it)
    {
      Variational::LinearForm R(v);
      R = makeResidual();
      R.assemble();
      const auto r = R.getVector();
      if (!r.allFinite()) break;
      const Real rNorm = r.norm();
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

      const Real e0 = energy();
      const GeoData u0 = u.getData();
      Real alpha = Real(1);
      bool accepted = false;
      for (Index ls = 0; ls < params.lineSearchBudget; ++ls)
      {
        u.getData() = u0 + alpha * dx;
        const Real e = energy();
        if (std::isfinite(e)
            && e <= e0 * (Real(1) + params.meritRelativeTolerance)
            && isCurvedMoveValid(mesh, u, fe, params.minDetFloor))
        {
          accepted = true;
          break;
        }
        alpha *= Real(0.5);
      }
      if (!accepted)
      {
        u.getData() = u0;
        break;
      }
      ++rep.iterations;
      if (alpha * dx.norm() <= params.stepTolerance) break;
    }

    // Commit: move the parametric mesh by the solved displacement.
    moveMesh(mesh, u, fe, interfaceAttribute);
    return rep;
  }
}

#endif
