/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_WNGIR_PARAMETERS_H
#define RODIN_ADAPTATION_WNGIR_PARAMETERS_H

#include <algorithm>
#include <cstddef>

#include "Rodin/Configure.h"
#include "Rodin/Geometry/Types.h"
#include "Rodin/Types.h"

namespace Rodin::Adaptation::WNGIR
{
  /// @brief Runtime parameters controlling WNGIR assembly and iteration.
  struct Parameters
  {
      /// @brief Coefficients, reference scale and admissible geometry of the model.
      struct Model
      {
          Real h = 0; ///< Fixed background reference size (required).
          Real fit = 1; ///< @f$\kappa_F@f$, target-normal fitting stiffness.
          /// @brief Global strain-variation weights, each assembled with h.
          struct Distribution
          {
              Real deviatoric = Real(1e-4); ///< @f$\kappa_{\rm dev}@f$.
              Real divergence = Real(1e-2); ///< @f$\kappa_{\rm div}@f$, with 1/d normalization.
          };
          Distribution distribution; ///< Centered current-strain distribution.
          Real distortion = 10; ///< @f$Q_{\max}@f$, relative-distortion budget.
          Real jacobian =
            Real(1e-2); ///< @f$j_{\mathrm{safe}}@f$, relative Jacobian floor.
          Real hinge = Real(10); ///< @f$\widehat\mu@f$, hinge/model-decrease ratio.
          Real qualityGuard = Real(0.1); ///< Guard fraction of the identity margins.
          Real jacobianWeight = 1; ///< Relative Jacobian-hinge row weight.
          Real distortionWeight = 1; ///< Relative distortion-hinge row weight.
          Real robustScale =
            0; ///< Positive fixes Welsch scale; zero selects it automatically.
      };

      /// @brief Accuracy requirements and independent work budgets.
      struct Convergence
      {
          /// @brief Geometric success, stationarity and best-effort exit tolerances.
          struct Tolerance
          {
              Real geometric = 0; ///< Zero selects @f$h^{p+1}@f$; positive overrides it.
              Real innerRelative =
                Real(1e-3); ///< Stationarity residual relative to force.
              Real innerAbsolute = Real(1e-12); ///< Absolute stationarity allowance.
              Real linearRelative =
                Real(1e-6); ///< Linear residual tolerance for all backends.
              Real energy = Real(1e-8); ///< Relative energy-change stagnation threshold.
              Real step = 0; ///< Absolute accepted-motion stagnation threshold.
              Real stepOverH = Real(5e-4); ///< Accepted-motion/reference-size threshold.
          };

          /// @brief Independent limits for fitting, Newton, linear solves and searches.
          struct Iterations
          {
              std::size_t outer = 30; ///< Maximum outer fitting iterations.
              std::size_t inner = 15; ///< Maximum Newton corrections per outer iteration.
              std::size_t linear = 1000; ///< Maximum iterations per CG solve only.
              std::size_t backtracks = 32; ///< Maximum outer trial halvings.
              std::size_t stagnation = 5; ///< Consecutive small steps or energy changes.
          };

          Tolerance tolerance; ///< Accuracy and stagnation thresholds.
          Iterations iterations; ///< Work budgets and stagnation persistence.
      };

      /// @brief Directional scaling and actual trial acceptance.
      struct Globalization
      {
          bool directionalNewton = true; ///< Scale the model without a level-set Hessian.
          Real maxStepOverH =
            0; ///< Positive caps predictor motion/h; zero is unrestricted.
          Real armijo = Real(1e-4); ///< Armijo sufficient-decrease coefficient.
      };

      /// @brief Available local linear-solver backends.
      enum class LinearSolver
      {
        CG,
        SparseLU,
        MUMPS
      };
      /// @brief Linear backend selection and thread policy.
      struct Linear
      {
#ifdef RODIN_USE_MUMPS
          LinearSolver solver = LinearSolver::MUMPS; ///< Default direct backend.
#else
          LinearSolver solver = LinearSolver::SparseLU; ///< Default direct backend.
#endif
          Index threads = 0; ///< Zero leaves the backend thread policy unchanged.
      };

      /// @brief Integration and independent geometry-validation orders.
      struct Quadrature
      {
          /**
           * @brief Calibrated P1 or provisional higher-order volume integration policy.
           * @param feOrder Polynomial order of the displacement finite element.
           * @param transformationOrder Polynomial order of the geometry transformation.
           * @param simplex Whether the entity is a simplex; tensor-product degree one is not P1.
           * @returns Two for affine P1; otherwise at least eight. Not an exactness guarantee.
           */
          static size_t getCellOrder(size_t feOrder, size_t transformationOrder = 1,
            bool simplex = true)
          {
            constexpr size_t affineOrder = 2, nonlinearMinimum = 8;
            return feOrder <= 1 && transformationOrder == 1 && simplex
              ? affineOrder : std::max(nonlinearMinimum, 2 * feOrder);
          }

          /**
           * @brief Surface order including the non-polynomial composed level set.
           * Eight is calibrated for affine P1; the higher-order minimum is provisional.
           * Neither is an exactness guarantee for arbitrary analytic coefficients.
           * @param feOrder Polynomial order of the displacement finite element.
           * @param transformationOrder Polynomial order of the geometry transformation.
           * @param simplex Whether the entity is a simplex.
           * @returns Automatic interface integration order.
           */
          static size_t getInterfaceOrder(size_t feOrder, size_t transformationOrder = 1,
            bool simplex = true)
          {
            constexpr size_t affineOrder = 8, nonlinearMinimum = 12;
            return feOrder <= 1 && transformationOrder == 1 && simplex
              ? affineOrder : std::max(nonlinearMinimum, 2 * feOrder + 2);
          }

          /**
           * @brief Independent geometric-validation order for an FE order.
           * @param feOrder Polynomial order of the displacement finite element.
           * @returns Automatic geometric-validation sampling order.
           */
          static size_t getValidationOrder(size_t feOrder)
          {
            constexpr size_t minimumOrder = 32;
            return std::max(minimumOrder, 2 * feOrder + 4);
          }

          /// @brief Surface integration order; specific settings override the common order.
          /// @param feOrder Displacement finite-element order.
          /// @param transformationOrder Geometry transformation order.
          /// @param simplex Whether the entity is a simplex.
          /// @returns Selected surface integration order.
          size_t getSurfaceOrder(size_t feOrder, size_t transformationOrder = 1,
            bool simplex = true) const
          {
            return surface > 0 ? surface : order > 0 ? order
              : getInterfaceOrder(feOrder, transformationOrder, simplex);
          }

          /// @brief Volume integration order; specific settings override the common order.
          /// @param feOrder Displacement finite-element order.
          /// @param transformationOrder Geometry transformation order.
          /// @param simplex Whether the entity is a simplex.
          /// @returns Selected volume integration order.
          size_t getVolumeOrder(size_t feOrder, size_t transformationOrder = 1,
            bool simplex = true) const
          {
            return volume > 0 ? volume : order > 0 ? order
              : getCellOrder(feOrder, transformationOrder, simplex);
          }

          /// @brief Independent actual-quality sampling order.
          /// @param feOrder Displacement finite-element order.
          /// @param transformationOrder Geometry transformation order.
          /// @param simplex Whether the entity is a simplex.
          /// @returns Independent quality order: two for affine P1, otherwise at least sixteen.
          size_t getQualityOrder(size_t feOrder, size_t transformationOrder = 1,
            bool simplex = true) const
          {
            constexpr size_t affineOrder = 2, nonlinearMinimum = 16;
            return quality > 0 ? quality : feOrder <= 1 && transformationOrder == 1 && simplex
              ? affineOrder : std::max(nonlinearMinimum, 2 * feOrder + 4);
          }

          std::size_t order = 0; ///< Zero selects automatic integration orders.
          std::size_t surface = 0; ///< Zero uses the common or automatic surface order.
          std::size_t volume = 0; ///< Zero uses the common or automatic volume order.
          std::size_t quality = 0; ///< Zero selects independent automatic quality sampling; vertices are always checked.
          std::size_t validation = 0; ///< Zero selects an independent validation order.
      };

      Model model; ///< Fitting, distribution and quality model.
      Convergence convergence; ///< Tolerances and work budgets.
      Globalization globalization; ///< Predictor scaling and outer acceptance.
      Linear linear; ///< Linear backend and thread policy.
      Quadrature quadrature; ///< Integration and validation orders.
      Optional<Geometry::Attribute> interfaceAttribute; ///< Marked facets to fit.
      bool trace = false; ///< Print diagnostics; accepted geometry is always validated.
      bool traceQualityWitness =
        false; ///< Record same-point quality predictions at each trial's limiting cell.
  };
}

#endif
