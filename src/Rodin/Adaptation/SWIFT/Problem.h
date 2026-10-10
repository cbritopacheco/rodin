/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_SWIFT_PROBLEM_H
#define RODIN_ADAPTATION_SWIFT_PROBLEM_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "Rodin/Solver/NewtonSolver.h"
#include "HingeProblem.h"
#include "LinearSolver.h"
#include "Rodin/Variational/Integrator.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/Trace.h"
#include "Rodin/Variational/IdentityMatrix.h"
#include "Rodin/Variational/VectorFunction.h"

#include "Rodin/Location.h"

#include "../CellDeformation.h"
#include "ForwardDecls.h"
#include "Loss.h"
#include "Parameters.h"
#include "DirectionalNewton.h"
#include "Distribution.h"
#include "HingeForce.h"
#include "HingeMetric.h"
#include "Report.h"
#include "FittingTensor.h"
#include "FittingForce.h"

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Sampled physical displacement norm, independent of FE coefficients.
   * @param[in] mesh Reference mesh.
   * @param[in] fes Displacement finite element space.
   * @param[in] cells Cells on which to sample physical displacement.
   * @param[in] step Displacement to evaluate.
   * @param[in] validationOrder Optional quadrature-order override.
   * @returns Largest sampled absolute displacement component.
   */
  template <class Mesh, class FES, class Displacement>
  Real getPhysicalDisplacementNorm(const Mesh& mesh, const FES& fes,
    const std::vector<Index>& cells, const Displacement& step,
    std::size_t validationOrder = 0)
  {
    if constexpr (requires { step.acquire(); })
      step.acquire();
    Real maximum = 0;
    for (const Index index : cells)
    {
      const auto cell = mesh.getCell(index);
      const Geometry::Polytope::Traits traits(cell->getGeometry());
      const auto evaluate = [&](const Geometry::Point& point) {
        const auto value = step.getValue(point);
        Real norm = 0;
        for (std::size_t component = 0; component < mesh.getDimension(); ++component)
        {
          if (!std::isfinite(value(component)))
            return std::numeric_limits<Real>::infinity();
          norm = std::max(norm, std::abs(value(component)));
        }
        return norm;
      };
      for (std::size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
      {
        maximum =
          std::max(maximum, evaluate(Geometry::Point(*cell, traits.getVertex(vertex))));
      }
      const auto& fe = fes.getFiniteElement(cell->getDimension(), index);
      // This is a polynomial FE field, not the composed level-set residual.
      const auto& qf = QF::PolytopeQuadratureFormula::get(validationOrder > 0
          ? validationOrder
          : std::max<size_t>(6, 2 * fe.getOrder() + 4),
        cell->getGeometry());
      const auto& quadrature = cell->getQuadrature(qf);
      for (std::size_t q = 0; q < quadrature.getSize(); ++q)
        maximum = std::max(maximum, evaluate(quadrature.getPoint(q)));
    }
    return maximum;
  }

  /**
   * @brief Strain-distributed Welsch Implicit-interface Fitting Technique (SWIFT).
   *
   * Fits marked facets to the zero level set of @f$\phi@f$ using the metric
   * @f$M=F+D@f$ (Fitting, Distribution), affine quadratic quality hinges,
   * directional Newton scaling and an outer Armijo line search. The trial
   * function owns the accumulated displacement; the solver does not move mesh
   * vertices, change transformations or remesh. See @ref guides-swift for setup,
   * parameters, monitoring and metric extension examples.
   *
   * @par Architecture
   * A native @ref Rodin::Variational::Problem assembles the predictor metric
   * and fitting force. @ref SWIFT::HingeProblem assembles the inner tangent
   * and stationarity residual; @ref Rodin::Solver::NewtonSolver drives its
   * corrections through a merit-backtracking policy.
   * One retained @ref SWIFT::LinearSolver handles both problems and their global
   * mean-strain subtraction. The outer loop owns directional scaling, actual
   * geometry checks, Armijo acceptance and response reporting.
   * Configuration and results are exposed as Parameters and Report. Their
   * standalone types remain available to configure a solve before its
   * finite-element types are known. Component types share the SWIFT namespace.
   * Supplied trial/test storage selects the backend; currently only local Eigen
   * assembly is implemented. PETSc storage is rejected at compile time.
   *
   * @section swift-usage Usage
   * @subsection swift-usage-setup Prepare the mesh and displacement space
   * The following two-dimensional example fits the internal line
   * @f$x_1=0.5@f$ to @f$x_1=0.55@f$. No boundary is fixed, so a uniform
   * translation remains available. The application marks the interface before
   * constructing the solver; SWIFT does not classify it automatically.
   * @code{.cpp}
   * using namespace Rodin;
   * using namespace Rodin::Geometry;
   * using namespace Rodin::Variational;
   *
   * // Build a unit-square triangular mesh and its required facet connectivity.
   * constexpr std::size_t n = 5;
   * constexpr Real referenceSpacing = Real(1) / Real(n - 1);
   * constexpr Geometry::Attribute Interface = 2;
   * auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
   * mesh.scale(referenceSpacing);
   * mesh.getConnectivity().compute(2, 1);
   * mesh.getConnectivity().compute(1, 0);
   * mesh.getConnectivity().compute(1, 2);
   *
   * // Mark facets whose vertices lie on the initial interface.
   * for (auto face = mesh.getFace(); face; ++face)
   * {
   *   bool onInterface = true;
   *   for (const Index vertex : face->getVertices())
   *     onInterface &= std::abs(mesh.getVertexCoordinates(vertex)(0) - Real(0.5))
   *                    < Real(1e-12);
   *   if (onInterface)
   *     mesh.setAttribute({1, face->getIndex()}, Interface);
   * }
   *
   * // u owns the accumulated displacement; v tests increments on the same space.
   * P1<Math::SpatialVector<Real>, LocalMesh> space(mesh, 2);
   * TrialFunction u(space);
   * TestFunction v(space);
   * Adaptation::SWIFT::Problem fitting(u, v);
   * @endcode
   * These snippets use the public headers `Rodin/Adaptation.h`,
   * `Rodin/Geometry.h` and `Rodin/Variational.h`, with `<cmath>` for the
   * coordinate test and `<iostream>` for the monitor below.
   *
   * @subsection swift-usage-solve Configure, monitor and solve
   * @code{.cpp}
   * decltype(fitting)::Parameters parameters;
   * // Use the background spacing, held fixed as the displacement changes.
   * parameters.model.h = referenceSpacing;
   * // Set motion stiffness separately from quality recovery and admissibility.
   * parameters.model.fit = 1;
   * parameters.model.distribution.deviatoric = Real(1e-4);
   * parameters.model.distribution.divergence = Real(1e-2);
   * parameters.model.hinge = 10;
   * parameters.model.distortion = 10;
   * parameters.model.jacobian = Real(0.01);
   * // These are work budgets; zero geometric tolerance selects the automatic target.
   * parameters.convergence.iterations.outer = 30;
   * parameters.convergence.iterations.inner = 15;
   * parameters.convergence.tolerance.geometric = 0;
   * fitting.setParameters(parameters).setInterfaceAttribute(Interface);
   *
   * // Observe accepted steps and the final report without modifying the problem.
   * fitting.setMonitor([](const decltype(fitting)::Report& state) {
   *   std::cout << state.iterations << ' ' << state.energy << ' '
   *             << state.geometricSup << ' ' << state.minJ << ' '
   *             << state.maxQRel << '\n';
   * });
   *
   * // Supply the target level set and its gradient; no target Hessian is needed.
   * RealFunction phi([](const Geometry::Point& point) {
   *   return point.x() - Real(0.55);
   * });
   * VectorFunction gradient(Real(1), Real(0));
   * const auto report = fitting.solve(phi, gradient);
   * std::cout << report.getReasonString() << '\n';
   * const bool targetReached =
   *   report.reason == decltype(fitting)::Report::Reason::GeometricTarget;
   * // solve() updates u, but leaves mesh vertices and transformations unchanged.
   * const auto& displacement = u.getSolution();
   * @endcode
   * A best-effort result may remain quality-admissible without reaching the
   * target. Interpretation requires both the stopping reason and the geometric
   * and quality fields; successful linear or inner solves are not geometric
   * convergence. The final monitor snapshot can repeat the last accepted step.
   *
   * @subsection swift-usage-boundaries Boundary conditions and metric extensions
   * If a boundary has been marked with an application-defined `FixedBoundary`
   * attribute, its current position can be held fixed by adding a homogeneous
   * increment condition before solving:
   * @code{.cpp}
   * const auto fixed =
   *   DirichletBC(u, VectorFunction(Real(0), Real(0))).on(FixedBoundary);
   * fitting += fixed;
   * // Remove the same condition when this boundary is allowed to move again.
   * fitting -= fixed;
   * @endcode
   * Nonzero increment values and identification conditions are not supported.
   * Additional bilinear metric terms can be added and removed through
   * getMetric(); they alter the direction model, not the fitting objective or
   * force. For example, adding a mass term also penalizes translations and is
   * therefore a different motion model. See @ref guides-swift-extension.
   *
   * @subsection swift-usage-lifetime Ownership, reuse and geometry application
   * The mesh is obtained from the trial function's finite element space.
   * Boundary attributes and @f$\phi,\nabla\phi@f$ are supplied by the
   * application. The mesh, space and trial/test functions must outlive the
   * solver. Homogeneous Dirichlet conditions constrain increments, not the
   * accumulated displacement. Applying the displacement to mesh geometry
   * remains a separate application step. A repeated solve takes the existing
   * displacement as its initial guess; replacing the mesh or space requires
   * a new solver. See @ref guides-swift for a complete setup example.
   *
   * @subsection swift-realization Architecture
   * Native bilinear and linear forms assemble fitting and distribution;
   * native variational problems own the predictor and inner systems.
   * The inner model is solved through @ref Rodin::Solver::NewtonSolver.
   * The local Eigen backend supports CG, SparseLU and optional MUMPS.
   * Added bilinear forms from getMetric() change @f$M_k@f$, not @f$E@f$ or
   * @f$f_k@f$. Custom additions must preserve the required solve properties;
   * an arbitrary bilinear form need not be positive semidefinite.
   *
   * @section swift-model Model
   * @subsection swift-energy Setting and fitting energy
   * Let @f$\Omega_0\subset\mathbb R^d@f$ be the fixed background domain,
   * @f$\Gamma_0@f$ its marked facet skeleton, and @f$V_h@f$ the continuous
   * vector-valued finite element increment space, with homogeneous essential
   * conditions. For the current displacement @f$u_k@f$, define
   * @f[
   * T_k(x)=x+u_k(x),\qquad A_k=I+\nabla u_k,\qquad
   * j_k=\det A_k,\qquad r_k=\phi\circ T_k,\qquad
   * g_k=(\nabla\phi)\circ T_k.
   * @f]
   * The robust objective and its negative first variation are
   * @f[
   * E(u)=G^{-2}\int_{\Gamma_0}\rho_\sigma(\phi(x+u(x)))\,dS_x,
   * \qquad \rho_\sigma(r)=\frac{\sigma^2}{2}
   *       \left(1-e^{-r^2/\sigma^2}\right),
   * @f]
   * @f[
   * f_k[z]=-DE(u_k)[z]
   *       =-G^{-2}\int_{\Gamma_0}e^{-r_k^2/\sigma^2}r_k\,g_k\cdot z\,dS_x.
   * @f]
   * Both the robust scale @f$\sigma@f$ and gradient normalization @f$G@f$
   * remain fixed throughout a solve. @f$G@f$ is the maximum target-gradient
   * magnitude sampled on the initial interface. Integration uses the reference
   * facet measure, not the deformed surface measure. The supplied gradient must
   * equal @f$\nabla\phi@f$ for @f$f_k@f$ to be the negative first variation.
   *
   * @subsection swift-metric Fitting and distribution metric
   * The frozen bilinear form @f$M_k=F_k+D_k@f$ is assembled as
   * @f[
   * F_k[v,z]=\kappa_F G^{-2}\int_{\Gamma_0}
   *          (g_k\cdot v)(g_k\cdot z)\,dS_x,
   * @f]
   * @f[
   * D_k[v,z]=h\int_{\Omega_k}\left[
   *   \kappa_{\rm dev}(\operatorname{dev}\varepsilon(v)
   *        -\overline{\operatorname{dev}\varepsilon(v)}):
   *       (\operatorname{dev}\varepsilon(z)
   *        -\overline{\operatorname{dev}\varepsilon(z)})
   *   +\frac{\kappa_{\rm div}}{d}
   *       (\operatorname{div}v-\overline{\operatorname{div}v})
   *       (\operatorname{div}z-\overline{\operatorname{div}z})
   *   \right]dy.
   * @f]
   * The symmetric and deviatoric parts are defined by
   * @f[
   * \operatorname{sym}L=\frac{L+L^T}{2},\qquad
   * \operatorname{dev}L=L-\frac{\operatorname{tr}L}{d}I.
   * @f]
   * Here @f$h@f$ is the fixed reference mesh size. @f$F_k@f$ is the Hessian of the half-squared
   * fitting residual with the level-set Hessian omitted, scaled by
   * @f$\kappa_F@f$. It is not robust-weighted. @f$D_k@f$ is a pullback of
   * centered strain variation in the current configuration. Here
   * @f$\varepsilon(v)=\operatorname{sym}\nabla_yv@f$ and overbars are global
   * current-volume averages. Coherent affine current-coordinate motions cost
   * zero; independently varying element strains generally do not. Continuity
   * and the global solve distribute interface motion through shared DOFs.
   *
   * @subsection swift-predictor Predictor and directional Newton scaling
   * First solve @f$M_k[p_k,z]=f_k[z]@f$ on the constrained increment space.
   * Directional Newton rescales this physical predictor using
   * @f[
   * a_k=\begin{cases}
   * \dfrac{f_k[p_k]}{c_k}, & \lambda_{\max}=0,\\
   * \min\left\{\dfrac{f_k[p_k]}{c_k},\,
   * \dfrac{\ell_{\max}}{\|p_k\|_{\infty,\mathrm{samp}}}\right\},
   * & \lambda_{\max}>0,
   * \end{cases}
   * \qquad \bar p_k=a_kp_k,\qquad \bar M_k=M_k/a_k.
   * @f]
   * By default, @f$a_k=f_k[p_k]/c_k@f$ is unrestricted. A positive optional
   * motion bound is @f$\ell_{\max}=h\lambda_{\max}@f$, where
   * @f$\lambda_{\max}>0@f$ is a dimensionless motion limit; the norm is the sampled
   * componentwise maximum of the physical field, not its coefficient vector.
   * The scalar curvature is selected from
   * @f[
   * c_k^{\mathrm N}=G^{-2}\int_{\Gamma_0}
   *   \rho_\sigma''(r_k)(g_k\cdot p_k)^2\,dS_x,
   * \qquad
   * c_k^{\mathrm{GN}}=G^{-2}\int_{\Gamma_0}
   *   e^{-r_k^2/\sigma^2}(g_k\cdot p_k)^2\,dS_x.
   * @f]
   * The Newton value is used when positive; otherwise the positive weighted
   * Gauss--Newton value is used.
   * No level-set Hessian is required. This is a scalar line model, not a full
   * outer Newton iteration. Without directional scaling, @f$a_k=1@f$.
   *
   * @subsection swift-inner Affine quality model and inner Newton problem
   * Relative distortion is defined, for @f$\det A>0@f$, by
   * @f[
   * Q(A)=\frac{|A|_F^2}{d(\det A)^{2/d}}.
   * @f]
   * Its minimum is one; rotations and positive uniform scaling do not change
   * it. The Jacobian floor separately guards against collapse and inversion.
   * At the frozen outer state, the slacks for an increment @f$v@f$ are
   * @f[
   * s_J(v)=j_k-j_{\mathrm{safe}}+Dj(A_k)[\nabla v],\qquad
   * s_Q(v)=Q_{\max}-Q(A_k)-DQ(A_k)[\nabla v].
   * @f]
   * The guard widths are
   * @f[
   * \delta_J=\gamma(1-j_{\mathrm{safe}}),\qquad
   * \delta_Q=\gamma(Q_{\max}-1).
   * @f]
   * With @f$(t)_+=\max(t,0)@f$, the affine squared-hinge penalty is
   * @f[
   * B_k(v)=\frac{\mu_k}{2}\sum_K\sum_{s\in\mathcal S_K}
   *   \left[w^J_{K,s}\kappa_J\left(1-\frac{s_J(v)}{\delta_J}\right)_+^2
   *        +w^Q_{K,s}\kappa_Q\left(1-\frac{s_Q(v)}{\delta_Q}\right)_+^2\right]_{K,s},
   * \qquad \mu_k=\widehat\mu\frac{f_k[\bar p_k]}{2|\Omega_0|}.
   * @f]
   * The set @f$\mathcal S_K@f$ contains the same closed-form reference covering
   * points as the actual-quality checks. For each constraint, @ref QualitySamples
   * computes capped nonlinear guard penetration at the current geometry and
   * the full physical predictor, retaining the larger risk @f$r^c_{K,s}@f$.
   * With mapped cell mass @f$W_K@f$ and witness count @f$N_K@f$, it sets
   * @f[
   * w^c_{K,s}=\frac{W_K}{N_K}
   * \left(\frac12+\frac{N_Kr^c_{K,s}}{2\sum_t r^c_{K,t}}\right),
   * \qquad c\in\{J,Q\}.
   * @f]
   * If all risks vanish, weights are equal. Risks are capped at 100; an
   * inverted predictor receives maximal distortion risk without evaluating
   * distortion there. Both measures preserve cell mass and are strictly
   * positive. They are frozen throughout the inner solve, not differentiated.
   * Validation remains unweighted. This adaptive discrete penalty is not a
   * quadrature approximation of a fixed volume penalty and does not certify
   * quality between witnesses. Here @f$\gamma\in(0,1)@f$ specifies
   * the guard fraction. The inner problem is
   * @f[
   * \min_{v\in V_h}\Psi_k(v),\qquad
   * \Psi_k(v)=\tfrac12\bar M_k[v,v]-f_k[v]+B_k(v).
   * @f]
   * Starting from @f$\bar p_k@f$, Newton corrections solve the stationarity equation
   * @f[
   * (\bar M_k+D^2B_k(v_m))\,\delta v_m
   *       =f_k-\bar M_kv_m-DB_k(v_m).
   * @f]
   * Merit backtracking decreases @f$\Psi_k@f$. The model is piecewise quadratic:
   * within a fixed active hinge set its Hessian is constant. Affine slacks
   * contribute rank-one row products, not @f$D^2j@f$ or @f$D^2Q@f$.
   * These finite penalties are not hard constraints; actual nonlinear quality
   * is checked again by the outer line search. Inactive hinges require no
   * correction. Inner stopping uses the Euclidean stationarity residual on
   * the free solve coordinates:
   * @f[
   * \|\nabla\Psi_k(v_m)\|_2\le\tau_{\mathrm{abs}}
   *       +\tau_{\mathrm{rel}}\|f_k\|_2.
   * @f]
   *
   * @subsection swift-outer Outer globalization and geometric stopping
   * A descending inner direction @f$v_k@f$ is accepted with
   * @f$\alpha_k\in\{1,1/2,1/4,\ldots\}@f$ only if
   * @f[
   * E(u_k+\alpha_kv_k)\le E(u_k)-c_A\alpha_k f_k[v_k],\qquad
   * j(u_k+\alpha_kv_k)>j_{\mathrm{safe}},\qquad
   * Q(u_k+\alpha_kv_k)<Q_{\max}.
   * @f]
   * Quality inequalities are evaluated on the geometric validation samples.
   * The update is @f$u_{k+1}=u_k+\alpha_kv_k@f$. Geometric success requires
   * @f[
   * D_{\infty,\mathrm{samp}}(u)=\max_{x\in\mathcal S_{\Gamma_0}}
   *   \frac{|\phi(x+u(x))|}{|\nabla\phi(x+u(x))|}
   *   \le\varepsilon_{\mathrm{geom}}.
   * @f]
   * The sample set includes interface vertices and validation quadrature points.
   * The default target is @f$h^{p+1}@f$, where @f$p@f$ is the interface finite
   * element order. This is a sampled local distance estimate, not a certified
   * Hausdorff bound. Small accepted steps, small relative energy changes and
   * exhausted iteration budgets are best-effort exits, not target convergence.
   *
   * @section swift-parameters Parameters
   * The values below describe the default model and stopping policy. The
   * reference size must be supplied by the application. Automatically selected
   * scales are fixed for the duration of a solve; they are not recomputed from
   * the shrinking or expanding mesh. The corresponding C++ fields are listed
   * in @ref guides-swift-controls and @ref Parameters.
   *
   * | Parameter | C++ | Description | Value |
   * |-----------|-----|-------------|-------|
   * | @f$h@f$ | `model.h` | Background reference size; scales distribution, motion limit and geometric target | Required |
   * | @f$\kappa_F@f$ | `model.fit` | Target-normal motion stiffness, not a force multiplier | @f$1@f$ |
   * | @f$\kappa_{\rm dev}@f$ | `model.distribution.deviatoric` | Cost of deviatoric-strain variation about its global mean | @f$10^{-4}@f$ |
   * | @f$\kappa_{\rm div}@f$ | `model.distribution.divergence` | Cost of divergence variation about its global mean, normalized by @f$d@f$ | @f$10^{-2}@f$ |
   * | @f$\widehat\mu@f$ | `model.hinge` | Hinge strength relative to predicted fitting improvement | @f$10@f$ |
   * | @f$\kappa_J@f$ | `model.jacobianWeight` | Relative Jacobian-hinge row weight | @f$1@f$ |
   * | @f$\kappa_Q@f$ | `model.distortionWeight` | Relative distortion-hinge row weight | @f$1@f$ |
   * | @f$\gamma@f$ | `model.qualityGuard` | Guard fraction defining hinge activation margins | @f$0.1@f$ |
   * | @f$j_{\mathrm{safe}}@f$ | `model.jacobian` | Relative Jacobian floor in the hinge and actual acceptance checks | @f$10^{-2}@f$ |
   * | @f$Q_{\max}@f$ | `model.distortion` | Maximum admissible relative distortion | @f$10@f$ |
   * | @f$\sigma@f$ | `model.robustScale` | Welsch robustness scale in level-set units | Automatic (`0`) |
   * | @f$\lambda_{\max}@f$ | `globalization.maxStepOverH` | Optional predictor-motion limit divided by @f$h@f$; zero disables it | @f$0@f$ |
   * | @f$c_A@f$ | `globalization.armijo` | Armijo sufficient-decrease coefficient | @f$10^{-4}@f$ |
   * | @f$\varepsilon_{\mathrm{geom}}@f$ | `convergence.tolerance.geometric` | Maximum sampled normalized-residual target | @f$h^{p+1}@f$, automatic (`0`) |
   * | @f$\tau_{\mathrm{rel}}@f$ | `convergence.tolerance.innerRelative` | Inner stationarity tolerance relative to the fitting force norm | @f$10^{-3}@f$ |
   * | @f$\tau_{\mathrm{abs}}@f$ | `convergence.tolerance.innerAbsolute` | Absolute allowance in the inner stationarity test | @f$10^{-12}@f$ |
   * | @f$\tau_{\mathrm{lin}}@f$ | `convergence.tolerance.linearRelative` | Relative linear residual tolerance | @f$10^{-6}@f$ |
   * | @f$N_{\mathrm{outer}}@f$ | `convergence.iterations.outer` | Maximum outer iterations | @f$30@f$ |
   * | @f$N_{\mathrm{inner}}@f$ | `convergence.iterations.inner` | Maximum Newton corrections per outer iteration | @f$15@f$ |
   * | @f$N_{\mathrm{CG}}@f$ | `convergence.iterations.linear` | Maximum iterations per CG solve; not a direct-solver limit | @f$1000@f$ |
   * | @f$N_{\mathrm{bt}}@f$ | `convergence.iterations.backtracks` | Maximum outer trial halvings | @f$32@f$ |
   * | @f$\varepsilon_{\mathrm{step}}/h@f$ | `convergence.tolerance.stepOverH` | Small accepted-motion threshold for a best-effort exit | @f$5\times10^{-4}@f$ |
   * | @f$\varepsilon_E@f$ | `convergence.tolerance.energy` | Relative energy-change threshold for a best-effort exit | @f$10^{-8}@f$ |
   * | @f$N_{\mathrm{stag}}@f$ | `convergence.iterations.stagnation` | Consecutive small steps or energy changes before stagnation exit | @f$5@f$ |
   *
   * @subsection swift-motion-controls Motion and robustness
   * @ref Parameters::Model::fit controls normal-motion stiffness in the
   * metric, not the magnitude of the fitting force. Increasing @f$\kappa_F@f$
   * therefore does not mean stronger attraction to the interface.
   * @ref Parameters::Model::Distribution::deviatoric and
   * @ref Parameters::Model::Distribution::divergence control variations
   * of shape-changing strain and local volume-change rate about their global
   * means. Uniform global shear, stretching and scaling are not penalized.
   * Larger weights distribute motion more coherently but can impede fitting.
   * Relative weights select the predictor direction; directional Newton
   * selects its physical length. Neither term enforces quality by itself.
   *
   * @ref Parameters::Model::h is the fixed background reference size, not the
   * size of the deformed elements. It scales distribution, the automatic
   * geometric target and the predictor-motion cap.
   * @ref Parameters::Globalization::maxStepOverH bounds the scaled
   * predictor in units of @f$h@f$ before quality recovery when positive; zero
   * leaves directional Newton unrestricted. The hinge solve can
   * alter that predictor, and outer backtracking checks the resulting motion.
   * @ref Parameters::Model::robustScale sets @f$\sigma@f$ in level-set units.
   * Residuals much larger than @f$\sigma@f$ have reduced influence on the
   * force. This protects against poorly classified observations, but excessive
   * downweighting can also weaken useful fitting motion. Zero selects an
   * automatic scale, which remains fixed during the solve.
   *
   * @subsection swift-quality-controls Quality recovery and admissibility
   * @ref Parameters::Model::hinge scales the finite hinge penalty relative to
   * the predicted fitting improvement. Increasing @f$\widehat\mu@f$ strengthens
   * recovery when hinges are active; it has no effect on inactive hinges and
   * does not guarantee nonlinear feasibility.
   * @ref Parameters::Model::qualityGuard sets their activation margins relative
   * to the identity margins. A wider guard reacts earlier, but also changes
   * the slack normalization and hence the penalty curvature.
   * @ref Parameters::Model::jacobianWeight and @ref Parameters::Model::distortionWeight set the
   * relative importance of the Jacobian and distortion hinge rows.
   *
   * @ref Parameters::Model::distortion is a spendable distortion budget, not a
   * quantity that must be minimized. @ref Parameters::Model::jacobian protects
   * against relative volume collapse through the hinge and actual quality
   * checks. The same Jacobian floor is used throughout. Increasing a Jacobian floor restricts compression,
   * whereas @f$Q@f$ is insensitive to positive isotropic scaling.
   *
   * @subsection swift-stopping-controls Accuracy and work limits
   * @ref Parameters::Convergence::Tolerance::geometric specifies the sampled
   * distance target; zero selects @f$h^{p+1}@f$.
   * @ref Parameters::Convergence::Tolerance::innerRelative and
   * @ref Parameters::Convergence::Tolerance::innerAbsolute control stationarity of the
   * direction problem, not geometric accuracy. Linear residual accuracy is
   * controlled separately by @ref Parameters::Convergence::Tolerance::linearRelative.
   * @ref Parameters::Convergence::Iterations::outer and
   * @ref Parameters::Convergence::Iterations::inner are work limits, not convergence
   * certificates. Persistent small accepted motion or small energy changes
   * produce best-effort exits. The report distinguishes these exits from a
   * geometric target hit and records fitting energy, maximum sampled error,
   * quality and iteration counts. Defaults and diagnostic fields are listed
   * in @ref guides-swift-controls and @ref guides-swift-results.
   *
   * @section swift-theory Theory
   * @subsection swift-structure Metric structure and solvability
   * For positive metric weights and @f$j_k>0@f$, both terms are integrals of
   * squares. With nonnegative quadrature weights this property also holds for
   * the assembled form:
   * @f[
   * M_k[v,v]=F_k[v,v]+D_k[v,v]\ge0.
   * @f]
   * The unresolved space consists of motions invisible to both observations
   * and distribution. Formally, its continuous characterization is
   * @f[
   * \mathcal N_k=\left\{v\in V_h:\;
   *   g_k\cdot v=0\text{ on }\Gamma_0,\quad
   *   v\circ T_k^{-1}\text{ is globally affine in }\Omega_k\right\}.
   * @f]
   * The discrete kernel instead depends on the quadrature samples. A
   * compatible predictor is unique only modulo this kernel. Selecting a
   * representative does not require penalizing translations, rotations or
   * isotropic strain in the metric. The present solver imposes no nullspace
   * gauge; singular systems remain subject to the selected linear backend.
   * With both centered weights positive, Korn's inequality controls H1 modulo
   * global affine motions on a fixed connected bounded Lipschitz domain.
   * If the divergence weight is zero, the kernel also includes conformal
   * motions and is infinite-dimensional in the continuous two-dimensional setting.
   *
   * On a fixed finite-dimensional quotient, a symmetric positive-semidefinite
   * form is positive definite after its complete kernel has been removed.
   * A mesh-uniform @f$H^1@f$ estimate is a stronger, separate statement:
   * @f[
   * M_k[v,v]\ge c\inf_{w\in\mathcal N_k}
   *             \|v-w\|_{H^1(\Omega_0)}^2.
   * @f]
   * Such an estimate requires control of the deformation, observation geometry,
   * discrete space and quadrature, with @f$c>0@f$ independent of refinement.
   * It does not follow from positive semidefiniteness alone.
   *
   * @subsection swift-convergence Scope of convergence diagnostics
   * Affine squared hinges are convex. Their active-region Hessian therefore
   * preserves positive semidefiniteness:
   * @f[
   * D^2B_k(v)[z,z]\ge0,\qquad
   * D^2\Psi_k(v)[z,z]=\bar M_k[z,z]+D^2B_k(v)[z,z]\ge0.
   * @f]
   * On a fixed active region with an invertible reduced Hessian, the inner
   * stationarity equation is linear and an exact full Newton correction solves
   * that local quadratic problem. Active-region changes and backtracking
   * prevent an unconditional one-correction claim.
   *
   * For a descending direction, accepted outer steps satisfy
   * @f[
   * E(u_k)-E(u_{k+1})\ge c_A\alpha_k f_k[v_k]>0.
   * @f]
   * This is sufficient decrease, not a proof that the complete outer sequence
   * reaches the geometric target. Likewise, @f$D_{\infty,\mathrm{samp}}@f$
   * measures a finite sample set, not the Hausdorff distance. Sampled quality
   * checks do not certify every point of a curved element, and a target hit
   * does not establish the refinement law
   * @f[
   * D_\infty(h)\le C h^{p+1}
   * @f]
   * with a refinement-independent constant @f$C@f$.
   */
  template <class TrialFunctionType, class TestFunctionType>
  class Problem
  {
    public:
      /// @brief Model, convergence and backend configuration.
      using Parameters = SWIFT::Parameters;
      /// @brief Accepted-geometry and iteration diagnostics.
      using Report = SWIFT::Report;

    private:
      using Displacement = std::remove_reference_t<
        decltype(std::declval<TrialFunctionType&>().getSolution())>;
      using ProblemType = std::decay_t<decltype(Variational::Problem(
        std::declval<TrialFunctionType&>(), std::declval<TestFunctionType&>()))>;
      using LinearSystemType = typename ProblemType::LinearSystemType;
      static_assert(
        std::is_same_v<typename FormLanguage::Traits<LinearSystemType>::OperatorType,
          Math::SparseMatrix<Real>> &&
          std::is_same_v<typename FormLanguage::Traits<LinearSystemType>::VectorType,
            Math::Vector<Real>>,
        "SWIFT::Problem currently supports local Eigen storage only; "
        "PETSc/MPI assembly is not implemented.");
      using BilinearFormType = std::decay_t<decltype(Variational::BilinearForm(
        std::declval<TrialFunctionType&>(), std::declval<TestFunctionType&>()))>;
      using LinearFormType = std::decay_t<decltype(Variational::LinearForm(
        std::declval<TestFunctionType&>()))>;

      using SpatialVec = Math::SpatialVector<Real>;
      using SpatialMat = Math::SpatialMatrix<Real>;
      using Loss = SWIFT::Loss;
      using Hinge = HingeState;
      using HingeProblem = SWIFT::HingeProblem<TrialFunctionType, TestFunctionType>;
      using LinearSolver = SWIFT::LinearSolver<LinearSystemType>;
      using Distribution =
        SWIFT::Distribution<TrialFunctionType, TestFunctionType, Displacement>;
      using HingeMetric =
        SWIFT::HingeMetric<TrialFunctionType, TestFunctionType, Displacement>;
      using HingeForce = SWIFT::HingeForce<TestFunctionType, Displacement>;

    public:
      /// @brief Observer of accepted-step and final diagnostics.
      using Monitor = std::function<void(const Report&)>;

      /**
       * @brief Observes accepted steps and the final report without changing the solve.
       * @param[in] monitor Observer, or an empty optional to remove it.
       * @returns This problem.
       */
      Problem& setMonitor(Optional<Monitor> monitor)
      {
        m_monitor = std::move(monitor);
        return *this;
      }

    private:
      /// @brief Mesh-validity summary of a displacement.
      struct AdmissibilityState
      {
          /// @brief Smallest Jacobian over all quality witnesses.
          Real minJ = std::numeric_limits<Real>::infinity();

          /// @brief Largest Jacobian over all quality witnesses.
          Real maxJ = -std::numeric_limits<Real>::infinity();

          /// @brief Largest relative distortion over the non-inverted points.
          Real maxQ = 0;

          /// @brief Number of quality witnesses at or below the Jacobian floor.
          std::size_t inadmissibleCount = 0;
          Real affineMinJ = std::numeric_limits<Real>::infinity();
          Real affineMaxQ = 0;
          Index qualityCell = 0;
          Real qualityCurrent = 0;
          Real qualityLinearChange = 0;
          Real qualityQuadraticChange = 0;
      };

      /// @brief Interface-fit summary of a displacement.
      struct SurfaceState
      {
          /// @brief Fixed-scale energy used by the current Armijo search.
          Real energy = 0;

          /// @brief Measure of the whole interface.
          Real totalLen = 0;

          /// @brief Root-mean-square residual over the complete interface.
          Real residualRMS = 0;

          /// @brief Supremum residual over the complete interface.
          Real residualSup = 0;
      };

      /**
       * @brief Kink of the discrete interface across facet boundaries.
       *
       * These statistics describe the classified interface only; they do not
       * alter the geometric target or the solver tolerances.
       */
      struct NormalJump
      {
          /// @brief Measure-weighted RMS angle between adjacent facet normals.
          Real rms = 0;

          /// @brief Largest angle between adjacent facet normals.
          Real max = 0;
      };

    public:
      /**
       * @brief Constructs the fitting problem from trial and test functions.
       * @param[in,out] du Trial function owning accumulated displacement.
       * @param[in] v Matching test function.
       */
      Problem(TrialFunctionType& du, TestFunctionType& v)
        : m_u(&du.getSolution()),
          m_trialUUID(du.getUUID()),
          m_duStep(du.getFiniteElementSpace()),
          m_vStep(v.getFiniteElementSpace()),
          m_metricProblem(m_duStep, m_vStep),
          m_hingeProblem(m_duStep, m_vStep),
          m_linearSolver(m_hingeProblem),
          m_distributionForm(m_duStep, m_vStep),
          m_fittingMetric(m_duStep, m_vStep),
          m_fittingForce(m_vStep),
          m_additionalMetric(du, v)
      {}

      /**
       * @brief Copying a problem with owned solver bindings is disabled.
       * @param other Problem that cannot be copied.
       */
      Problem(const Problem& other) = delete;
      /**
       * @brief Copy assignment of a problem with owned solver bindings is disabled.
       * @param other Problem that cannot be assigned.
       * @returns No value, since assignment is deleted.
       */
      Problem& operator=(const Problem& other) = delete;

      /**
       * @brief Additional metric terms, assembled afresh at each outer iteration.
       * These augment @f$F+D@f$; they do not change the fitting energy or force.
       * @returns Modifiable additional bilinear form.
       */
      BilinearFormType& getMetric()
      {
        return m_additionalMetric;
      }

      /**
       * @brief Inspects additional metric terms.
       * @returns Read-only additional bilinear form.
       */
      const BilinearFormType& getMetric() const
      {
        return m_additionalMetric;
      }

      /**
       * @brief Selects the marked interface on the displacement mesh.
       * @param[in] attribute Interface facet attribute.
       * @returns This problem.
       */
      Problem& setInterfaceAttribute(Geometry::Attribute attribute)
      {
        m_parameters.interfaceAttribute = attribute;
        return *this;
      }

      /**
       * @brief Holds the current boundary displacement fixed during fitting.
       *
       * Only value-prescribing, homogeneous conditions on the constructor's
       * trial function are supported. They constrain increments, not accumulated
       * displacement. Nonzero values and identification constraints are rejected.
       * @param[in] condition Homogeneous condition on the constructor's trial function.
       * @returns This problem.
       */
      Problem& operator+=(const Variational::DirichletBCBase<Real>& condition)
      {
        if (condition.getOperand().getUUID() != m_trialUUID)
          Alert::Exception() << "SWIFT boundary conditions must use its trial function."
                             << Alert::Raise;
        m_boundaryConditions.add(condition);
        return *this;
      }

      /**
       * @brief Sets SWIFT runtime parameters.
       * @param[in] parameters Model, convergence and backend controls.
       * @returns This problem.
       */
      Problem& setParameters(const Parameters& parameters)
      {
#ifndef RODIN_USE_MUMPS
        if (parameters.linear.solver == Parameters::LinearSolver::MUMPS)
          Alert::Exception() << "SWIFT MUMPS solves require RODIN_USE_MUMPS."
                             << Alert::Raise;
#endif
        if (!std::isfinite(parameters.model.fit) || !(parameters.model.fit > Real(0)))
          Alert::Exception() << "SWIFT requires a finite positive fitting weight."
                             << Alert::Raise;
        for (const Real value :
          {parameters.model.jacobianWeight, parameters.model.distortionWeight,
            parameters.model.distribution.deviatoric,
            parameters.model.distribution.divergence,
            parameters.convergence.tolerance.innerAbsolute,
            parameters.convergence.tolerance.energy,
            parameters.convergence.tolerance.step,
            parameters.convergence.tolerance.stepOverH, parameters.model.robustScale})
        {
          if (!std::isfinite(value) || value < Real(0))
            Alert::Exception()
              << "SWIFT weights and tolerances must be finite and nonnegative."
              << Alert::Raise;
        }
        if (!std::isfinite(parameters.convergence.tolerance.innerRelative) ||
          !(parameters.convergence.tolerance.innerRelative > Real(0)) ||
          !std::isfinite(parameters.convergence.tolerance.linearRelative) ||
          !(parameters.convergence.tolerance.linearRelative > Real(0)) ||
          !(parameters.globalization.armijo > Real(0) &&
            parameters.globalization.armijo < Real(1)) ||
          parameters.convergence.iterations.inner == 0 ||
          parameters.convergence.iterations.stagnation == 0)
          Alert::Exception()
            << "SWIFT requires valid residual tolerances, budgets and line search."
            << Alert::Raise;
        if (!std::isfinite(parameters.globalization.maxStepOverH) ||
          parameters.globalization.maxStepOverH < Real(0))
          Alert::Exception()
            << "Directional Newton requires a finite nonnegative step/h bound."
            << Alert::Raise;
        if (!(parameters.model.qualityGuard > Real(0) &&
              parameters.model.qualityGuard < Real(1)) ||
          !(parameters.model.jacobian > Real(0) && parameters.model.jacobian < Real(1)) ||
          !(parameters.model.distortion > Real(1)) ||
          !std::isfinite(parameters.model.distortion) ||
          !std::isfinite(parameters.model.hinge) || parameters.model.hinge < Real(0) ||
          !std::isfinite(parameters.convergence.tolerance.geometric) ||
          parameters.convergence.tolerance.geometric < Real(0))
          Alert::Exception()
            << "SWIFT requires 0 < guard,jacobian < 1, finite distortion > 1 and "
               "nonnegative penalty/fit tolerance."
            << Alert::Raise;
        m_parameters = parameters;
        return *this;
      }

      /**
       * @brief Inspects current controls.
       * @returns The current SWIFT parameters.
       */
      const Parameters& getParameters() const
      {
        return m_parameters;
      }

      /**
       * @brief Inspects the most recent solve.
       * @returns The last completed fitting report.
       */
      const Report& getReport() const
      {
        return m_report;
      }

      /**
       * @brief Solves the SWIFT fitting problem on marked interface facets.
       *
       * The supplied sensitivity @p grad must equal the derivative of @p phi
       * at the moved quadrature points for the assembled force to be the exact
       * first variation of the line-search energy. An independently supplied
       * sensitivity is supported, but then defines a pseudo-gradient.
       * @param[in] phi Target level-set function.
       * @param[in] grad Target gradient in physical coordinates.
       * @returns Accepted-geometry and iteration diagnostics.
       */
      template <class PhiDerived, class GradDerived>
      Report solve(const Variational::RealFunctionBase<PhiDerived>& phi,
        const Variational::VectorFunctionBase<Real, GradDerived>& grad)
      {
        using Rodin::Index;
        const Parameters& p = m_parameters;
        Displacement& u = *m_u;
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        using Mesh = std::remove_cvref_t<decltype(mesh)>;
        const std::size_t meshDim = mesh.getDimension();

        Report rep;
        const Real h = p.model.h;
        if (!std::isfinite(h) || !(h > Real(0)))
          Alert::Exception() << "SWIFT requires a finite positive reference mesh size."
                             << Alert::Raise;
        assembleBoundaryConditions();
        const Real acceptedJacobianFloor = p.model.jacobian;
        const Real stepTol = p.convergence.tolerance.step;
        using Clock = std::chrono::steady_clock;
        auto secondsSince = [](Clock::time_point t0) -> Real {
          return std::chrono::duration<Real>(Clock::now() - t0).count();
        };
        auto setupTic = Clock::now();

        // ============================================================
        // Per-frame geometry tabulation (only u changes per iteration).
        // Field values come from GridFunction::getValue (cached), so the
        // tables store quadrature geometry only.
        // ============================================================
        std::vector<Index> validationCells;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const Index cellIndex = cell->getIndex();
          validationCells.push_back(cellIndex);
        }

        if (!p.interfaceAttribute)
        {
          rep.reason = Report::Reason::MissingInterface;
          return finish(rep);
        }

        // One bounding-volume locator per solve: the background mesh is fixed
        // for the frame, and the index builds lazily on first query.
        const Location::AABB<Mesh> locator(mesh);
        std::vector<Index> interfaceFacets;
        for (auto face = mesh.getFace(); face; ++face)
        {
          if (face->getAttribute() == *p.interfaceAttribute)
            interfaceFacets.push_back(face->getIndex());
        }
        if (interfaceFacets.empty())
        {
          rep.reason = Report::Reason::EmptyInterface;
          return finish(rep);
        }
        const auto normalJump = getNormalJump(mesh, fes, interfaceFacets, meshDim);
        rep.normalJumpRMS = normalJump.rms;
        rep.normalJumpMax = normalJump.max;
        const Real geometricScale = std::pow(
          h, fes.getFiniteElement(meshDim - 1, interfaceFacets.front()).getOrder() + 1);
        const Real geometricTarget = p.convergence.tolerance.geometric > Real(0)
          ? p.convergence.tolerance.geometric
          : geometricScale;
        if (!std::isfinite(geometricScale) || !(geometricScale > Real(0)) ||
          !std::isfinite(geometricTarget) || !(geometricTarget > Real(0)))
          Alert::Exception() << "SWIFT requires a finite positive geometric target."
                             << Alert::Raise;
        rep.geometricSupTarget = geometricTarget;
        const auto robustScale = getRobustScale(mesh, fes, phi, grad, interfaceFacets, h);
        if (!(robustScale.gradientScale > Real(0)) ||
          !std::isfinite(robustScale.gradientScale))
        {
          rep.reason = Report::Reason::DegenerateGradient;
          return finish(rep);
        }
        const Real sigma = robustScale.sigma;
        const Real levelSetMeshScale = h * robustScale.gradientScale;
        const Real dataNormalization =
          Real(1) / (robustScale.gradientScale * robustScale.gradientScale);
        const Loss loss(sigma);
        rep.sigma = sigma;
        rep.levelSetGradientScale = robustScale.gradientScale;
        const Real domainMeasure = getDomainMeasure(mesh, fes, validationCells);

        // ============================================================
        // Per-iteration field evaluations through the GridFunction.
        // ============================================================
        auto admissibility = [&](const Displacement& gf) {
          return getAdmissibilityState(mesh, fes, validationCells, gf, meshDim);
        };
        auto surfaceState = [&](const Displacement& gf) {
          return getSurfaceState(
            mesh, fes, gf, phi, interfaceFacets, loss, dataNormalization, locator);
        };
        auto recordSurfaceState = [&rep](const SurfaceState& state) {
          rep.energy = state.energy;
          rep.residualRMS = state.residualRMS;
          rep.residualSup = state.residualSup;
          rep.interfaceMeasure = state.totalLen;
        };

        rep.tSetup = secondsSince(setupTic);

        Displacement vK(fes);
        Displacement predictor(fes);
        Displacement scratch(fes);
        Displacement previousU(fes);
        Displacement uTrial(fes);

        // The interface coefficients are not polynomial, so the quadrature
        // order is given to the integrator explicitly, per face.
        const Variational::Integrator::OrderType surfaceOrder =
          [&](const Geometry::Polytope& face) -> std::size_t {
          const auto& fe = fes.getFiniteElement(face.getDimension(), face.getIndex());
          return p.quadrature.getSurfaceOrder(fe.getOrder(),
            face.getTransformation().getOrder(),
            Geometry::Polytope::Traits(face.getGeometry()).getVertexCount() ==
              face.getDimension() + 1);
        };

        AdmissibilityState initialAdm{};
        {
          initialAdm = admissibility(u);
          rep.minJ = initialAdm.minJ;
          rep.maxJ = initialAdm.maxJ;
          rep.maxQRel = initialAdm.maxQ;
        }
        if (initialAdm.inadmissibleCount > 0 ||
          initialAdm.minJ <= acceptedJacobianFloor ||
          initialAdm.maxQ >= p.model.distortion)
        {
          rep.reason = Report::Reason::InvalidInitialGeometry;
          return finish(rep);
        }

        const SurfaceState initialSurface = surfaceState(u);
        recordSurfaceState(initialSurface);
        if (!(initialSurface.totalLen > Real(0)))
        {
          rep.reason = Report::Reason::EmptyInterface;
          return finish(rep);
        }

        Real ePrev = initialSurface.energy;
        const auto traceStart = setupTic;
        // Validation concerns the accepted outer geometry, not the inner QP iterate.
        const auto traceGeometry = [&](const auto& geometry, std::size_t accepted,
                                     const char* phase) {
          const auto flags = std::cout.flags();
          const auto precision = std::cout.precision();
          std::cout << "      swift geometry: outer=" << accepted << "  phase=" << phase
                    << std::scientific
                    << std::setprecision(std::numeric_limits<Real>::max_digits10)
                    << "  geom_rms=" << geometry.rms << "  geom_sup=" << geometry.sup
                    << "  geom_c=" << geometry.sup / geometricScale
                    << "  normal_rms=" << geometry.normalRMS << "  h=" << h
                    << "  min_j=" << rep.minJ << "  max_qrel=" << rep.maxQRel
                    << "  inner_total=" << rep.innerIterations
                    << "  inner_last=" << rep.lastInnerIterations
                    << "  inner_converged=" << rep.innerConverged
                    << "  inner_residual=" << rep.innerResidual
                    << "  inner_relative_residual=" << rep.innerRelativeResidual
                    << "  inner_residual_tolerance=" << rep.innerResidualTolerance
                    << "  geom_sup_target=" << geometricTarget
                    << "  inactive_hinge_skips=" << rep.inactiveHingeSkips
                    << "  direct_analyses=" << rep.directAnalyses
                    << "  direct_factorizations=" << rep.directFactorizations
                    << "  fit_energy=" << rep.energy
                    << "  interface_measure=" << rep.interfaceMeasure
                    << "  inner_assembly_seconds=" << rep.tInnerAssembly
                    << "  inner_solve_seconds=" << rep.tInnerSolve
                    << "  inner_ls_seconds=" << rep.tInnerLineSearch
                    << "  seconds=" << secondsSince(traceStart) << std::endl;
          std::cout.flags(flags);
          std::cout.precision(precision);
        };
        {
          const auto geometry = getInterfaceGeometryState(
            mesh, fes, u, phi, grad, interfaceFacets, meshDim, locator);
          if (p.trace)
            traceGeometry(geometry, 0, "initial");
          rep.geometricRMS = geometry.rms;
          rep.geometricSup = geometry.sup;
          rep.geometricConstant = geometry.sup / geometricScale;
          rep.normalRMS = geometry.normalRMS;
          rep.qualityBudgetSatisfied = true;
          if (!std::isfinite(geometry.sup))
          {
            rep.reason = Report::Reason::InvalidGeometry;
            return finish(rep);
          }
          if (geometry.sup <= geometricTarget)
          {
            rep.geometricTargetReached = true;
            rep.reason = Report::Reason::GeometricTarget;
            return finish(rep);
          }
        }

        // ============================================================
        // Nonlinear iteration.
        // ============================================================
        std::size_t consecutiveSmallAcceptedSteps = 0;
        std::size_t consecutiveSmallEnergyChanges = 0;
        const auto recordLinearSolve = [&rep](std::size_t iterations, Real error) {
          rep.linearIterations += iterations;
          ++rep.linearSolveCount;
          rep.maxLinearIterations = std::max(rep.maxLinearIterations, iterations);
          rep.linearError = error;
        };
        for (; rep.iterations < p.convergence.iterations.outer; ++rep.iterations)
        {
          auto tic = Clock::now();
          FittingTensor tensor(grad, u, locator, p, dataNormalization, meshDim);
          const auto fittingIntegrand = Variational::Dot(tensor * m_duStep, m_vStep);
          auto fittingMetric = Variational::FaceIntegral(fittingIntegrand);
          // Native order propagation is valid for a known tensor on an affine
          // entity. Curved surface measures and unknown compositions retain
          // the non-polynomial policy. Explicit orders always take precedence.
          fittingMetric.setOrder(
            [&, integrand = fittingIntegrand](const Geometry::Polytope& face) {
              if (p.quadrature.surface > 0 || p.quadrature.order > 0)
                return surfaceOrder(face);
              const auto order = integrand.getOrder(face);
              const auto& transformation =
                mesh.getPolytopeTransformation(face.getDimension(), face.getIndex());
              if (order && transformation.getOrder() == 1)
                return *order;
              return surfaceOrder(face);
            });
          fittingMetric.over(*p.interfaceAttribute);
          FittingForce forceCoeff(
            phi, grad, u, locator, loss, dataNormalization, meshDim);
          auto surfaceForce = Variational::FaceIntegral(forceCoeff, m_vStep);
          surfaceForce.setOrder(surfaceOrder);
          surfaceForce.over(*p.interfaceAttribute);
          // The observation metric and the fitting force depend on the outer
          // displacement, not on the hinge increment, so they are assembled
          // here and reused by every correction below.
          m_fittingMetric = fittingMetric;
          m_fittingMetric.assemble();
          if (!m_additionalMetric.getLocalIntegrators().empty() ||
            !m_additionalMetric.getGlobalIntegrators().empty())
          {
            m_additionalMetric.assemble();
            m_fittingMetric.getOperator() += m_additionalMetric.getOperator();
          }
          const Real deviatoric = p.model.h * p.model.distribution.deviatoric;
          const Real divergence = p.model.h * p.model.distribution.divergence;
          size_t order = 0;
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            order = std::max(order,
              p.quadrature.getVolumeOrder(
                fes.getFiniteElement(meshDim, cell->getIndex()).getOrder(),
                cell->getTransformation().getOrder(),
                Geometry::Polytope::Traits(cell->getGeometry()).getVertexCount() ==
                  meshDim + 1));
          }
          const Distribution distribution(
            m_duStep, m_vStep, u, deviatoric, divergence, order);
          m_distributionForm = distribution;
          m_distributionForm.assemble();
          m_centering = distribution.getCentering();
          for (const auto& [dof, value] : m_fixedIncrementDOFs)
            m_centering.row(dof).setZero();
          m_linearSolver.setCentering(m_centering);
          m_fittingForce = surfaceForce;
          m_fittingForce.assemble();
          bool solveOk = true;
          Real predictorAction = Real(0);
          typename ProblemType::ProblemBodyType predictorBody(m_distributionForm);
          predictorBody = predictorBody + m_fittingMetric - m_fittingForce;
          m_metricProblem = predictorBody;
          assembleStep();
          rep.tAssembly += secondsSince(tic);

          Math::SparseMatrix<Real> fixedMetric;
          Math::Vector<Real> fixedForce;
          if constexpr (std::is_same_v<
                          typename FormLanguage::Traits<LinearSystemType>::OperatorType,
                          Math::SparseMatrix<Real>>)
          {
            {
              const auto& system = m_metricProblem.getLinearSystem();
              fixedMetric = system.getOperator();
              fixedForce = system.getVector();
            }
          }
          tic = Clock::now();
          solveOk = solvePredictor(vK, rep);
          recordLinearSolve(m_linearSolver.getIterations(), m_linearSolver.getError());
          rep.tSolve += secondsSince(tic);
          if (!solveOk)
          {
            rep.reason = Report::Reason::PredictorFailure;
            break;
          }

          predictor = vK;
          predictorAction = getSurfaceForceAction(mesh, fes, u, predictor, phi, grad,
            interfaceFacets, loss, dataNormalization, meshDim, locator);
          if (!(predictorAction > Real(0)) || !std::isfinite(predictorAction))
          {
            rep.reason = Report::Reason::NonDescentPredictor;
            break;
          }
          rep.predictorScale = Real(1);
          if (p.globalization.directionalNewton)
          {
            const auto curvatures = getSurfaceDirectionalCurvature(mesh, fes, u,
              predictor, phi, grad, interfaceFacets, loss, dataNormalization, locator);
            const Real norm = getPhysicalDisplacementNorm(
              mesh, fes, validationCells, predictor, p.quadrature.validation);
            rep.predictorScale =
              getDirectionalNewtonStep(predictorAction, curvatures.first,
                curvatures.second, norm, h * p.globalization.maxStepOverH);
            if (!(rep.predictorScale > Real(0)) || !std::isfinite(rep.predictorScale))
            {
              rep.reason = Report::Reason::InvalidScaling;
              break;
            }
            // The predictor, frozen quadratic form and hinge scale describe the same motion.
            predictor *= rep.predictorScale;
            vK = predictor;
            predictorAction *= rep.predictorScale;
            fixedMetric *= Real(1) / rep.predictorScale;
            m_distributionForm.getOperator() *= Real(1) / rep.predictorScale;
            m_fittingMetric.getOperator() *= Real(1) / rep.predictorScale;
            m_centering /= std::sqrt(rep.predictorScale);
            m_linearSolver.setCentering(m_centering);
            if (p.trace)
            {
              const auto precision = std::cout.precision();
              std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
                        << "      swift directional: outer=" << rep.iterations
                        << "  curvature=" << curvatures.first
                        << "  fitting_curvature=" << curvatures.second
                        << "  scale=" << rep.predictorScale
                        << "  step_over_h=" << norm * rep.predictorScale / h << '\n';
              std::cout.precision(precision);
            }
          }
          rep.predictorAction = predictorAction;

          {
            const Real modelDecrease = Real(0.5) * std::max(Real(0), predictorAction);
            const Real hingeCoefficient = domainMeasure > Real(0)
              ? p.model.hinge * modelDecrease / domainMeasure
              : Real(0);
            rep.hingeCoefficient = hingeCoefficient;
            const size_t innerIterations = p.convergence.iterations.inner;
            const Real predictorNorm =
              std::max(std::abs(predictor.max()), std::abs(predictor.min()));
            const Real residualScale = fixedForce.norm();
            const Real residualTolerance = p.convergence.tolerance.innerAbsolute +
              p.convergence.tolerance.innerRelative * residualScale;
            bool innerConverged = false;
            rep.lastInnerIterations = 0;
            rep.lastInnerAlpha = 0;
            rep.minInnerAlpha = 1;
            rep.fullInnerSteps = 0;
            if (!hasActiveHinges(
                  mesh, fes, validationCells, u, vK, meshDim, hingeCoefficient))
            {
              Math::Vector<Real> image = fixedMetric * vK.getData();
              image -= m_centering * (m_centering.transpose() * vK.getData());
              rep.innerResidual = (image - fixedForce).norm();
              rep.innerResidualTolerance = residualTolerance;
              rep.innerRelativeResidual = residualScale > Real(0)
                ? rep.innerResidual / residualScale
                : rep.innerResidual;
              innerConverged = rep.innerResidual <= residualTolerance;
              rep.innerConverged = innerConverged;
              rep.innerRelativeCorrection = Real(0);
              if (p.trace)
                std::cout << "        hinge skip: outer=" << rep.iterations
                          << "  reason=inactive-hinges  residual=" << rep.innerResidual
                          << "  rel=" << rep.innerRelativeResidual
                          << "  converged=" << innerConverged << '\n';
              if (innerConverged)
                ++rep.inactiveHingeSkips;
            }
            if (!innerConverged)
            {
              HingeMetric hingeMetric(m_duStep, m_vStep, u, vK, predictor, p, hingeCoefficient);
              HingeForce hingeForce(m_vStep, u, vK, predictor, p, hingeCoefficient);
              typename ProblemType::ProblemBodyType body(m_distributionForm);
              body = body + m_fittingMetric + hingeMetric - m_fittingForce - hingeForce;
              m_hingeProblem = body;
              m_hingeProblem.setState(vK)
                .setBoundaryDOFs(m_fixedIncrementDOFs)
                .setCentering(m_centering);
              using Newton = Solver::NewtonSolver<LinearSolver>;
              Newton newton(m_linearSolver);
              newton.setMaxIterations(innerIterations)
                .setAbsoluteTolerance(residualTolerance)
                .setRelativeTolerance(Real(0))
                .setStepTolerance(Real(0));
              const auto recordResidual = [&](Real residual, size_t inner) {
                const Real assemblySeconds = m_hingeProblem.getAssemblySeconds();
                rep.tAssembly += assemblySeconds;
                rep.tInnerAssembly += assemblySeconds;
                rep.innerResidual = residual;
                rep.innerResidualTolerance = residualTolerance;
                rep.innerRelativeResidual =
                  residualScale > Real(0) ? residual / residualScale : residual;
                innerConverged = std::isfinite(residual) && residual <= residualTolerance;
                if (p.trace)
                {
                  const auto precision = std::cout.precision();
                  std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
                            << "        hinge residual: outer=" << rep.iterations
                            << "  inner=" << inner << "  residual=" << residual
                            << "  force_norm=" << residualScale
                            << "  rel=" << rep.innerRelativeResidual
                            << "  tolerance=" << residualTolerance
                            << "  converged=" << innerConverged << '\n';
                  std::cout.precision(precision);
                }
              };
              newton.setMonitor([&](const typename Newton::Report& report) {
                if (report.converged || !std::isfinite(report.finalResidual))
                  recordResidual(report.finalResidual, report.iterations);
              });
              newton.setStepPolicy([&](auto&, auto& system, auto& newtonReport) {
                recordResidual(newtonReport.finalResidual, newtonReport.iterations);
                const size_t linearIterations = m_linearSolver.getIterations();
                const Real linearError = m_linearSolver.getError();
                recordLinearSolve(linearIterations, linearError);
                const Real innerSolveSeconds = m_linearSolver.getSolveSeconds();
                rep.tSolve += innerSolveSeconds;
                rep.tInnerSolve += innerSolveSeconds;
                solveOk = m_linearSolver.success();
                if (!solveOk)
                {
                  if (p.trace)
                    std::cout << "        hinge inner=" << (newtonReport.iterations + 1)
                              << "  outer=" << rep.iterations << "  linear_ok=0"
                              << "  cg_it=" << linearIterations
                              << "  cg_err=" << linearError << std::endl;
                  rep.reason = Report::Reason::LinearFailure;
                  return typename Newton::StepResult{false, false, Real(0)};
                }
                scratch.setData(system.getSolution());
                const Real correctionNorm =
                  std::max(std::abs(scratch.max()), std::abs(scratch.min()));
                const Real iterateNorm = std::max(std::abs(vK.max()), std::abs(vK.min()));
                const Real correctionScale = std::max(iterateNorm, predictorNorm);
                const Real relativeCorrection = correctionScale > Real(0)
                  ? correctionNorm / correctionScale
                  : (correctionNorm == Real(0) ? Real(0)
                                               : std::numeric_limits<Real>::infinity());
                rep.innerRelativeCorrection = relativeCorrection;
                ++rep.innerIterations;
                ++rep.lastInnerIterations;
                rep.maxInnerIterations =
                  std::max(rep.maxInnerIterations, rep.lastInnerIterations);
                Real innerAlpha = Real(1);
                {
                  const auto meritStart = Clock::now();
                  const auto merit = [&](const Displacement& increment) {
                    Math::Vector<Real> image;
                    image = fixedMetric * increment.getData();
                    image -=
                      m_centering * (m_centering.transpose() * increment.getData());
                    const Real quadratic = Real(0.5) * increment.getData().dot(image);
                    const Real force = fixedForce.dot(increment.getData());
                    const Real quality = getHingeEnergy(mesh, fes, validationCells, u,
                      increment, predictor, meshDim, hingeCoefficient);
                    return std::pair<Real, Real>{quadratic - force + quality,
                      std::abs(quadratic) + std::abs(force) + std::abs(quality)};
                  };
                  const auto before = merit(vK);
                  const Real slope = -system.getVector().dot(scratch.getData());
                  bool meritAccepted = false;
                  std::size_t backtracks = 0;
                  Real after = std::numeric_limits<Real>::infinity();
                  for (; backtracks < InnerMaxBacktracks; ++backtracks)
                  {
                    uTrial = scratch;
                    uTrial *= innerAlpha;
                    uTrial += vK;
                    const auto trial = merit(uTrial);
                    after = trial.first;
                    const Real rounding = MeritRoundoffFactor *
                      std::numeric_limits<Real>::epsilon() *
                      std::max(before.second, trial.second);
                    if (std::isfinite(before.first) && std::isfinite(after) &&
                      std::isfinite(slope) && slope < Real(0) &&
                      after <= before.first +
                          InnerArmijoCoefficient * innerAlpha * std::min(Real(0), slope) +
                          rounding)
                    {
                      meritAccepted = true;
                      break;
                    }
                    innerAlpha *= BacktrackingReduction;
                  }
                  rep.innerBacktracks += backtracks;
                  rep.tInnerLineSearch += secondsSince(meritStart);
                  if (p.trace)
                  {
                    const auto precision = std::cout.precision();
                    std::cout << std::setprecision(
                                   std::numeric_limits<Real>::max_digits10)
                              << "        hinge merit: outer=" << rep.iterations
                              << "  inner=" << (newtonReport.iterations + 1)
                              << "  before=" << before.first << "  after=" << after
                              << "  slope=" << slope << "  alpha=" << innerAlpha
                              << "  backtracks=" << backtracks
                              << "  accepted=" << meritAccepted << std::endl;
                    std::cout.precision(precision);
                  }
                  if (!meritAccepted)
                  {
                    if (p.trace)
                      std::cout << "        hinge inner=" << (newtonReport.iterations + 1)
                                << "  outer=" << rep.iterations << "  linear_ok=1"
                                << "  merit_ok=0  converged=0  cg_it=" << linearIterations
                                << "  cg_err=" << linearError << std::endl;
                    rep.reason = Report::Reason::InnerLineSearchFailure;
                    solveOk = false;
                    return typename Newton::StepResult{false, false, Real(0)};
                  }
                }
                rep.lastInnerAlpha = innerAlpha;
                rep.minInnerAlpha = std::min(rep.minInnerAlpha, innerAlpha);
                if (innerAlpha >= Real(1) - FullStepTolerance)
                  ++rep.fullInnerSteps;
                if (p.trace)
                {
                  const auto precision = std::cout.precision();
                  std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
                            << "        hinge inner=" << (newtonReport.iterations + 1)
                            << "  outer=" << rep.iterations << "  corr=" << correctionNorm
                            << "  iterate=" << iterateNorm
                            << "  rel=" << relativeCorrection << "  alpha=" << innerAlpha
                            << "  linear_ok=1  cg_it=" << linearIterations
                            << "  cg_err=" << linearError
                            << "  mu_eff=" << hingeCoefficient << std::endl;
                  std::cout.precision(precision);
                }
                scratch *= innerAlpha;
                vK += scratch;

                return typename Newton::StepResult{true, false, scratch.getData().norm()};
              });
              newton.solve(vK);
              innerConverged = newton.converged();
              // Certify the last accepted iterate without spending another correction.
              if (solveOk &&
                newton.getReport().reason == Newton::ConvergedReason::MaxIterations)
              {
                m_hingeProblem.assemble();
                recordResidual(
                  m_hingeProblem.getLinearSystem().getVector().norm(), innerIterations);
              }
            }
            rep.innerConverged = innerConverged;
            if (solveOk && !innerConverged)
            {
              rep.reason = Report::Reason::InnerIterationLimit;
              solveOk = false;
            }
          }
          if (!solveOk)
            break;

          if (!(std::isfinite(vK.max()) && std::isfinite(vK.min())))
          {
            rep.reason = Report::Reason::NonfiniteDirection;
            break;
          }

          Real directionAction = getSurfaceForceAction(mesh, fes, u, vK, phi, grad,
            interfaceFacets, loss, dataNormalization, meshDim, locator);
          const Real predictorNorm =
            std::max(std::abs(predictor.max()), std::abs(predictor.min()));
          Real directionNorm = std::max(std::abs(vK.max()), std::abs(vK.min()));
          rep.directionAction = directionAction;
          rep.descentRatio =
            predictorAction > Real(0) ? directionAction / predictorAction : Real(0);
          rep.directionNormRatio =
            predictorNorm > Real(0) ? directionNorm / predictorNorm : Real(0);
          if (!(std::isfinite(directionAction) && directionAction > Real(0)))
          {
            rep.reason = Report::Reason::NonDescentDirection;
            break;
          }

          // ---- Nonlinear line search on TRUE geometry ----
          tic = Clock::now();
          Real alpha = Real(1);
          bool accepted = false;
          std::size_t backtracks = 0;
          AdmissibilityState adm{};
          Real eTrial = std::numeric_limits<Real>::infinity();
          SurfaceState trialSurface{};
          {
            previousU = u;
            for (; backtracks <= p.convergence.iterations.backtracks; ++backtracks)
            {
              // uTrial = previousU + alpha * vK
              uTrial = vK;
              uTrial *= alpha;
              uTrial += previousU;
              if ((uTrial.getData().array() == previousU.getData().array()).all())
                break;
              bool jOK = true;
              bool qOK = true;
              {
                adm = p.trace ? getAdmissibilityState(
                                  mesh, fes, validationCells, uTrial, meshDim, &previousU)
                              : admissibility(uTrial);
                jOK = adm.inadmissibleCount == 0 && adm.minJ > acceptedJacobianFloor;
                qOK = adm.maxQ < p.model.distortion;
              }
              bool eOK = true;
              if (jOK && qOK)
              {
                trialSurface = surfaceState(uTrial);
                eTrial = trialSurface.energy;
                const Real sufficientDecrease =
                  p.globalization.armijo * alpha * directionAction;
                eOK = std::isfinite(eTrial) && eTrial <= ePrev - sufficientDecrease;
              }
              if (p.trace)
              {
                const auto precision = std::cout.precision();
                std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
                          << "      swift trial: outer=" << rep.iterations
                          << "  alpha=" << alpha << "  affine_min_j=" << adm.affineMinJ
                          << "  affine_max_qrel=" << adm.affineMaxQ
                          << "  min_j=" << adm.minJ << "  max_qrel=" << adm.maxQ
                          << "  j_ok=" << jOK << "  q_ok=" << qOK
                          << "  energy_checked=" << (jOK && qOK)
                          << "  e_ok=" << (jOK && qOK && eOK) << std::endl;
                if (p.traceQualityWitness)
                  std::cout << "      swift quality witness: outer=" << rep.iterations
                            << "  alpha=" << alpha << "  cell=" << adm.qualityCell
                            << "  current_q=" << adm.qualityCurrent
                            << "  actual_q=" << adm.maxQ
                            << "  linear_change=" << adm.qualityLinearChange
                            << "  quadratic_change=" << adm.qualityQuadraticChange
                            << "  remainder="
                            << adm.maxQ - adm.qualityCurrent - adm.qualityLinearChange
                            << "  current_margin="
                            << p.model.distortion - adm.qualityCurrent << std::endl;
                std::cout.precision(precision);
              }
              if (jOK && qOK && eOK)
              {
                u = uTrial;
                accepted = true;
                break;
              }
              if (!jOK)
                ++rep.jacobianRejections;
              if (!qOK)
                ++rep.distortionRejections;
              if (jOK && qOK && !eOK)
                ++rep.energyRejections;
              alpha *= BacktrackingReduction;
            }
          }
          rep.tLineSearch += secondsSince(tic);
          if (!accepted)
          {
            u = previousU;
            rep.reason = Report::Reason::LineSearchFailure;
            break;
          }

          rep.lastAlpha = alpha;
          {
            // Measure the accepted FE field, not its coefficient vector.
            scratch = u;
            scratch -= previousU;
            rep.acceptedStep = getPhysicalDisplacementNorm(
              mesh, fes, validationCells, scratch, p.quadrature.validation);
          }
          rep.minJ = adm.minJ;
          rep.maxJ = adm.maxJ;
          rep.maxQRel = adm.maxQ;

          const auto& surf = trialSurface;
          const Real eNow = eTrial;
          rep.backtracks += backtracks;
          rep.actualPredictedDecrease = alpha * directionAction > Real(0)
            ? (ePrev - eNow) / (alpha * directionAction)
            : Real(0);
          recordSurfaceState(surf);
          bool geometricTargetReached = false;
          {
            const auto geometry = getInterfaceGeometryState(
              mesh, fes, u, phi, grad, interfaceFacets, meshDim, locator);
            rep.geometricRMS = geometry.rms;
            rep.geometricSup = geometry.sup;
            rep.geometricConstant = geometry.sup / geometricScale;
            rep.normalRMS = geometry.normalRMS;
            rep.qualityBudgetSatisfied =
              rep.minJ > acceptedJacobianFloor && rep.maxQRel < p.model.distortion;
            rep.geometricTargetReached =
              rep.qualityBudgetSatisfied && geometry.sup <= geometricTarget;
            if (m_monitor)
            {
              auto snapshot = rep;
              ++snapshot.iterations;
              (*m_monitor)(snapshot);
            }
            if (p.trace)
              traceGeometry(geometry, rep.iterations + 1, "accepted");
            if (!std::isfinite(geometry.sup))
            {
              rep.reason = Report::Reason::InvalidGeometry;
              ++rep.iterations;
              break;
            }
            geometricTargetReached = geometry.sup <= geometricTarget;
          }
          if (!(surf.totalLen > Real(0)))
          {
            rep.reason = Report::Reason::EmptyInterface;
            ++rep.iterations;
            break;
          }
          if (p.trace)
            std::cout << "      swift it=" << std::setw(3) << rep.iterations
                      << "  E=" << std::scientific << std::setprecision(3) << eNow
                      << "  residualRMS=" << surf.residualRMS << "  residualRMS/(hG)="
                      << (levelSetMeshScale > Real(0)
                             ? surf.residualRMS / levelSetMeshScale
                             : Real(0))
                      << "  residualSup=" << surf.residualSup
                      << "  step/h=" << (h > Real(0) ? rep.acceptedStep / h : Real(0))
                      << "  linIt=" << rep.linearIterations << "  alpha=" << alpha
                      << "  muEff=" << rep.hingeCoefficient
                      << "  innerIt=" << rep.lastInnerIterations
                      << "  innerRel=" << rep.innerRelativeCorrection
                      << "  predictor_action=" << rep.predictorAction
                      << "  inner_residual=" << rep.innerResidual
                      << "  inner_relative_residual=" << rep.innerRelativeResidual
                      << "  desc=" << rep.descentRatio
                      << "  dirNorm=" << rep.directionNormRatio
                      << "  ared/pred=" << rep.actualPredictedDecrease
                      << "  bt=" << backtracks << "  rejJ=" << rep.jacobianRejections
                      << "  rejQ=" << rep.distortionRejections
                      << "  rejE=" << rep.energyRejections << "  min_j=" << rep.minJ
                      << "  max_j=" << rep.maxJ << "  max_Q=" << rep.maxQRel << '\n';

          if (geometricTargetReached)
          {
            rep.reason = Report::Reason::GeometricTarget;
            ++rep.iterations;
            break;
          }
          const Real stepThreshold = stepTol + h * p.convergence.tolerance.stepOverH;
          consecutiveSmallAcceptedSteps =
            rep.acceptedStep <= stepThreshold ? consecutiveSmallAcceptedSteps + 1 : 0;
          const Real energyScale =
            std::max({std::abs(ePrev), std::abs(eNow), std::numeric_limits<Real>::min()});
          const Real energyChange = std::abs(ePrev - eNow) / energyScale;
          consecutiveSmallEnergyChanges = p.convergence.tolerance.energy > Real(0) &&
              energyChange <= p.convergence.tolerance.energy
            ? consecutiveSmallEnergyChanges + 1
            : 0;
          if (consecutiveSmallAcceptedSteps >= p.convergence.iterations.stagnation ||
            consecutiveSmallEnergyChanges >= p.convergence.iterations.stagnation)
          {
            rep.reason =
              consecutiveSmallAcceptedSteps >= p.convergence.iterations.stagnation
              ? Report::Reason::SmallAcceptedSteps
              : Report::Reason::SmallEnergyChanges;
            ++rep.iterations;
            break;
          }
          ePrev = eNow;
        }

        const InterfaceGeometryState geometry = getInterfaceGeometryState(
          mesh, fes, u, phi, grad, interfaceFacets, meshDim, locator);
        rep.geometricRMS = geometry.rms;
        rep.geometricSup = geometry.sup;
        rep.geometricConstant = geometry.sup / geometricScale;
        rep.normalRMS = geometry.normalRMS;
        rep.qualityBudgetSatisfied =
          rep.minJ > acceptedJacobianFloor && rep.maxQRel < p.model.distortion;
        rep.geometricTargetReached =
          rep.qualityBudgetSatisfied && geometry.sup <= geometricTarget;
        if (p.trace)
          traceGeometry(geometry, rep.iterations, "final");

        return finish(rep);
      }

    private:
      /**
       * @brief Publishes final diagnostics to the optional monitor.
       * @param report Completed solve diagnostics.
       * @returns The stored final report.
       */
      Report finish(const Report& report)
      {
        m_report = report;
        if (m_monitor)
          (*m_monitor)(m_report);
        return m_report;
      }

      /// @brief Geometric discrepancy of the complete fitted interface.
      struct InterfaceGeometryState
      {
          Real rms = std::numeric_limits<Real>::infinity();
          Real sup = std::numeric_limits<Real>::infinity();
          Real normalRMS = std::numeric_limits<Real>::infinity();
      };

      /**
       * @brief Evaluates directional robust curvature without the level-set Hessian.
       * @param mesh Fixed reference mesh.
       * @param fes Displacement finite element space defining quadrature orders.
       * @param current Frozen outer displacement.
       * @param direction Proposed incremental displacement.
       * @param phi Target level-set function.
       * @param grad Target level-set gradient.
       * @param interfaceFacets Indices of the classified interface facets.
       * @param loss Fixed-scale robust loss.
       * @param normalization Energy and force normalization factor.
       * @param locator Target-evaluation point locator.
       * @returns Robust scalar curvature and its nonnegative influence-weighted counterpart.
       */
      template <class Mesh, class FES, class PhiType, class GradType, class LocatorType>
      std::pair<Real, Real> getSurfaceDirectionalCurvature(const Mesh& mesh,
        const FES& fes, const Displacement& current, const Displacement& direction,
        const PhiType& phi, const GradType& grad,
        const std::vector<Index>& interfaceFacets, const Loss& loss, Real normalization,
        const LocatorType& locator) const
      {
        if constexpr (requires { current.acquire(); })
          current.acquire();
        if constexpr (requires { direction.acquire(); })
          direction.acquire();
        if constexpr (requires { phi.acquire(); })
          phi.acquire();
        if constexpr (requires { grad.acquire(); })
          grad.acquire();
        std::vector<std::pair<Real, Real>> facetCurvatures(interfaceFacets.size());
#ifdef RODIN_USE_OPENMP
#pragma omp parallel
#endif
        {
          DeformationMap<Displacement, LocatorType> deformation(current, locator);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(interfaceFacets.size()); ++i)
          {
            const auto face = mesh.getFace(interfaceFacets[static_cast<size_t>(i)]);
            const auto& qf = getQuadrature(*face, fes);
            const auto& quadrature = face->getQuadrature(qf);
            for (size_t q = 0; q < quadrature.getSize(); ++q)
            {
              const auto& point = quadrature.getPoint(q);
              const Variational::IntegrationPoint ip(point, &qf, q);
              const auto value = direction.getValue(point);
              const ResidualState state(phi, grad, deformation, ip, loss);
              const Real residual = state.getResidual();
              const Real action = state.getGradient().dot(value);
              const Real weight =
                qf.getWeight(q) * point.getDistortion() * normalization * action * action;
              facetCurvatures[i].first += weight * loss.getCurvature(residual);
              facetCurvatures[i].second += weight * state.getWeight();
            }
          }
        }
        std::pair<Real, Real> curvature{};
        for (const auto& value : facetCurvatures)
        {
          curvature.first += value.first;
          curvature.second += value.second;
        }
        return curvature;
      }

      /**
       * @brief Selects assembly quadrature for the entity and displacement order.
       * Affine P1 uses calibrated surface/volume orders; higher-order or
       * curved entities use provisional non-polynomial policies.
       * Specific orders override the common integration order.
       * @param polytope Entity on which quadrature is requested.
       * @param fes Displacement finite element space defining quadrature orders.
       * @returns The cached surface or volume quadrature formula.
       */
      template <class FES>
      const QF::QuadratureFormulaBase& getQuadrature(
        const Geometry::Polytope& polytope, const FES& fes) const
      {
        const auto& fe =
          fes.getFiniteElement(polytope.getDimension(), polytope.getIndex());
        const bool isInterface = polytope.getDimension() < fes.getMesh().getDimension();
        const bool simplex =
          Geometry::Polytope::Traits(polytope.getGeometry()).getVertexCount() ==
          polytope.getDimension() + 1;
        const std::size_t order = isInterface
          ? m_parameters.quadrature.getSurfaceOrder(
              fe.getOrder(), polytope.getTransformation().getOrder(), simplex)
          : m_parameters.quadrature.getVolumeOrder(
              fe.getOrder(), polytope.getTransformation().getOrder(), simplex);
        return QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
      }

      /**
       * @brief Selects diagnostic quadrature independently of assembly quadrature.
       * @param polytope Entity on which quadrature is requested.
       * @param fes Displacement finite element space defining quadrature orders.
       * @returns The cached geometric-response quadrature formula.
       */
      template <class FES>
      const QF::QuadratureFormulaBase& getGeometricValidationQuadrature(
        const Geometry::Polytope& polytope, const FES& fes) const
      {
        const auto& fe =
          fes.getFiniteElement(polytope.getDimension(), polytope.getIndex());
        const std::size_t order = m_parameters.quadrature.validation > 0
          ? m_parameters.quadrature.validation
          : Parameters::Quadrature::getValidationOrder(fe.getOrder());
        return QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
      }

      /**
       * @brief Integrates the fixed reference domain measure.
       * @param mesh Fixed reference mesh.
       * @param fes Displacement finite element space defining quadrature orders.
       * @param validationCells Reference cell indices included in the check.
       * @returns The reference volume used to normalize the quadratic model.
       */
      template <class Mesh, class FES>
      Real getDomainMeasure(
        const Mesh& mesh, const FES& fes, const std::vector<Index>& validationCells) const
      {
        if constexpr (requires { mesh.getShard(); })
        {
          return mesh.getMeasure(mesh.getDimension());
        }
        else
        {
          Real measure = 0;
          for (const Index cellIndex : validationCells)
          {
            const auto cell = mesh.getCell(cellIndex);
            const auto& qf = getQuadrature(*cell, fes);
            const auto& quadrature = cell->getQuadrature(qf);
            for (std::size_t q = 0; q < quadrature.getSize(); ++q)
              measure += qf.getWeight(q) * quadrature.getPoint(q).getDistortion();
          }
          return measure;
        }
      }

      /**
       * @brief Evaluates the negative robust-energy first variation on a direction.
       * @param mesh Fixed reference mesh.
       * @param fes Displacement finite element space defining quadrature orders.
       * @param current Frozen outer displacement.
       * @param direction Proposed incremental displacement.
       * @param phi Target level-set function.
       * @param grad Target level-set gradient.
       * @param interfaceFacets Indices of the classified interface facets.
       * @param loss Fixed-scale robust loss.
       * @param normalization Energy and force normalization factor.
       * @param dimension Spatial dimension.
       * @param locator Target-evaluation point locator.
       * @returns The normalized fitting-force action.
       */
      template <class Mesh, class FES, class PhiType, class GradType, class LocatorType>
      Real getSurfaceForceAction(const Mesh& mesh, const FES& fes,
        const Displacement& current, const Displacement& direction, const PhiType& phi,
        const GradType& grad, const std::vector<Index>& interfaceFacets, const Loss& loss,
        Real normalization, std::size_t dimension, const LocatorType& locator) const
      {
        if constexpr (requires { current.acquire(); })
          current.acquire();
        if constexpr (requires { direction.acquire(); })
          direction.acquire();
        if constexpr (requires { phi.acquire(); })
          phi.acquire();
        if constexpr (requires { grad.acquire(); })
          grad.acquire();
        std::vector<Real> facetActions(interfaceFacets.size(), Real(0));
#ifdef RODIN_USE_OPENMP
#pragma omp parallel
#endif
        {
          FittingForce force(phi, grad, current, locator, loss, normalization, dimension);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(interfaceFacets.size()); ++i)
          {
            Real& facetAction = facetActions[static_cast<std::size_t>(i)];
            const Index facetIndex = interfaceFacets[static_cast<std::size_t>(i)];
            const auto face = mesh.getFace(facetIndex);
            const auto& qf = getQuadrature(*face, fes);
            const auto& quadrature = face->getQuadrature(qf);
            for (std::size_t q = 0; q < quadrature.getSize(); ++q)
            {
              const auto& point = quadrature.getPoint(q);
              const Variational::IntegrationPoint ip(point, &qf, q);
              facetAction += qf.getWeight(q) * point.getDistortion() *
                force.getValue(ip).dot(direction.getValue(point));
            }
          }
        }
        Real action = 0;
        for (const Real facetAction : facetActions)
          action += facetAction;
        return action;
      }

      /**
       * @brief Detects active affine hinges at the actual-quality witnesses.
       * @param mesh Fixed reference mesh.
       * @param fes Displacement finite element space defining quadrature orders.
       * @param cells Reference cell indices included in the check.
       * @param current Frozen outer displacement.
       * @param inner Current inner increment.
       * @param dimension Spatial dimension.
       * @param coefficient Effective hinge penalty coefficient.
       * @returns Whether any hinge is active or the frozen state is inadmissible.
       */
      template <class Mesh, class FES>
      bool hasActiveHinges(const Mesh& mesh, const FES& fes,
        const std::vector<Index>& cells, const Displacement& current,
        const Displacement& inner, std::size_t dimension, Real coefficient) const
      {
        bool active = false;
#ifdef RODIN_USE_OPENMP
#pragma omp parallel reduction(|| : active)
#endif
        {
          auto currentJacobian = Variational::Jacobian(current);
          auto innerJacobian = Variational::Jacobian(inner);
          CellDeformation deformation(dimension);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(cells.size()); ++i)
          {
            const auto cell = mesh.getCell(cells[static_cast<std::size_t>(i)]);
            const QualitySamples samples(*cell,
              fes.getFiniteElement(dimension, cell->getIndex()).getOrder(), m_parameters);
            samples.forEach(
              [&](const Variational::IntegrationPoint& ip, Real) {
              deformation.setDisplacementGradient(currentJacobian.getValue(ip));
              const Hinge state(
                deformation, innerJacobian.getValue(ip), m_parameters, coefficient);
              active = active || !state.isAdmissible() ||
                state.getJacobianHessian() != Real(0) ||
                state.getDistortionHessian() != Real(0);
            });
          }
        }
        return active;
      }

      /**
       * @brief Integrates the affine squared-hinge energy.
       * @param mesh Fixed reference mesh.
       * @param fes Displacement finite element space defining quadrature orders.
       * @param validationCells Reference cell indices included in the check.
       * @param current Frozen outer displacement.
       * @param inner Current inner increment.
       * @param predictor Frozen directionally scaled predictor defining the weights.
       * @param dimension Spatial dimension.
       * @param coefficient Effective hinge penalty coefficient.
       * @returns The penalty energy of the current inner increment.
       */
      template <class Mesh, class FES>
      Real getHingeEnergy(const Mesh& mesh, const FES& fes,
        const std::vector<Index>& validationCells, const Displacement& current,
        const Displacement& inner, const Displacement& predictor,
        std::size_t dimension, Real coefficient) const
      {
        Real energy = 0;
#ifdef RODIN_USE_OPENMP
#pragma omp parallel reduction(+ : energy)
#endif
        {
          auto currentJacobian = Variational::Jacobian(current);
          auto innerJacobian = Variational::Jacobian(inner);
          CellDeformation deformation(dimension);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(validationCells.size()); ++i)
          {
            const auto cell = mesh.getCell(validationCells[static_cast<std::size_t>(i)]);
            const QualitySamples samples(*cell,
              fes.getFiniteElement(dimension, cell->getIndex()).getOrder(), m_parameters);
            auto predictorJacobian = Variational::Jacobian(predictor);
            samples.forEachHinge(currentJacobian, predictorJacobian,
              [&](const Variational::IntegrationPoint& ip, Real weightJ, Real weightQ) {
              deformation.setDisplacementGradient(currentJacobian.getValue(ip));
              const Hinge state(
                deformation, innerJacobian.getValue(ip), m_parameters, coefficient);
              energy += state.getEnergy(m_parameters, coefficient, weightJ, weightQ);
            });
          }
        }
        return energy;
      }

      /**
       * @brief Checks the actual deformed cell geometry.
       * @param mesh Fixed reference mesh.
       * @param fes Displacement finite element space defining quadrature orders.
       * @param validationCells Reference cell indices included in the check.
       * @param u Displacement whose deformed interface or geometry is evaluated.
       * @param dimension Spatial dimension.
       * @param current Optional frozen displacement for affine predictions and quality witnesses.
       * @returns Sampled Jacobian and distortion extrema and violation counts.
       */
      template <class Mesh, class FES>
      AdmissibilityState getAdmissibilityState(const Mesh& mesh, const FES& fes,
        const std::vector<Index>& validationCells, const Displacement& u,
        std::size_t dimension, const Displacement* current = nullptr) const
      {
        if constexpr (requires { u.acquire(); })
          u.acquire();
        Real minJ = std::numeric_limits<Real>::infinity();
        Real maxJ = -std::numeric_limits<Real>::infinity();
        Real maxQ = Real(0);
        Real affineMinJ = std::numeric_limits<Real>::infinity();
        Real affineMaxQ = Real(0);
        if (current)
        {
          if constexpr (requires { current->acquire(); })
            current->acquire();
        }
        std::size_t inadmissibleCount = 0;
        struct Witness
        {
            Real actual = 0, current = 0, linear = 0, quadratic = 0;
            Index cell = 0;
        };
        const bool recordWitness = current && m_parameters.traceQualityWitness;
        std::vector<Witness> witnesses(recordWitness ? validationCells.size() : 0);
#ifdef RODIN_USE_OPENMP
#pragma omp parallel reduction(min : minJ, affineMinJ)                                   \
  reduction(max : maxJ, maxQ, affineMaxQ) reduction(+ : inadmissibleCount)
#endif
        {
          CellDeformation deformation(dimension);
          CellDeformation currentDeformation(dimension);
          auto displacementJacobian = Variational::Jacobian(u);
          auto currentJacobian = Variational::Jacobian(current ? *current : u);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(validationCells.size()); ++i)
          {
            const Index cellIndex = validationCells[static_cast<std::size_t>(i)];
            const auto cell = mesh.getCell(cellIndex);
            const auto& fe = fes.getFiniteElement(dimension, cellIndex);
            const QualitySamples samples(*cell, fe.getOrder(), m_parameters);
            samples.forEach([&](const Variational::IntegrationPoint& ip, Real) {
              deformation.setDisplacementGradient(displacementJacobian.getValue(ip));
              const Real j = deformation.getJacobian();
              minJ = std::min(minJ, j);
              maxJ = std::max(maxJ, j);
              if (!std::isfinite(j) || j <= m_parameters.model.jacobian ||
                !deformation.isAdmissible())
                ++inadmissibleCount;
              if (deformation.isAdmissible())
              {
                const Real quality = deformation.getRelativeDistortion();
                if (!std::isfinite(quality))
                  ++inadmissibleCount;
                else
                  maxQ = std::max(maxQ, quality);
              }
              if (current)
              {
                const auto gradient = currentJacobian.getValue(ip);
                const SpatialMat increment(displacementJacobian.getValue(ip) - gradient);
                currentDeformation.setDisplacementGradient(gradient);
                affineMinJ = std::min(affineMinJ,
                  currentDeformation.getJacobian() +
                    currentDeformation.getJacobianAction(increment));
                affineMaxQ = std::max(affineMaxQ,
                  currentDeformation.getRelativeDistortion() +
                    currentDeformation.getRelativeDistortionAction(increment));
                if (recordWitness && deformation.isAdmissible())
                {
                  auto& witness = witnesses[static_cast<size_t>(i)];
                  const Real actual = deformation.getRelativeDistortion();
                  if (actual > witness.actual)
                    witness = {actual, currentDeformation.getRelativeDistortion(),
                      currentDeformation.getRelativeDistortionAction(increment),
                      Real(0.5) *
                        currentDeformation.getRelativeDistortionSecondAction(
                          increment, increment),
                      cellIndex};
                }
              }
            });
          }
        }
        AdmissibilityState result{
          minJ, maxJ, maxQ, inadmissibleCount, affineMinJ, affineMaxQ};
        Witness limiting;
        for (const auto& witness : witnesses)
        {
          if (witness.actual > limiting.actual)
            limiting = witness;
        }
        result.qualityCell = limiting.cell;
        result.qualityCurrent = limiting.current;
        result.qualityLinearChange = limiting.linear;
        result.qualityQuadraticChange = limiting.quadratic;
        return result;
      }

      /**
       * @brief Evaluates the interface fit of a displacement.
       *
       * The residual is the level set read at the deformed image of each
       * interface quadrature point. Energy uses the robust loss, but residual
       * diagnostics include the complete interface without an active-weight cutoff.
       */
      /**
       * @brief Computes robust energy and interface residual diagnostics.
       * @param mesh Fixed reference mesh.
       * @param fes Displacement finite element space defining quadrature orders.
       * @param u Displacement whose deformed interface or geometry is evaluated.
       * @param phi Target level-set function.
       * @param interfaceFacets Indices of the classified interface facets.
       * @param loss Fixed-scale robust loss.
       * @param normalization Energy and force normalization factor.
       * @param locator Target-evaluation point locator.
       * @returns The energy, measure and residual norms at assembly quadrature.
       */
      template <class Mesh, class FES, class PhiType, class LocatorType>
      SurfaceState getSurfaceState(const Mesh& mesh, const FES& fes,
        const Displacement& u, const PhiType& phi,
        const std::vector<Index>& interfaceFacets, const Loss& loss, Real normalization,
        const LocatorType& locator) const
      {
        if constexpr (requires { u.acquire(); })
          u.acquire();
        if constexpr (requires { phi.acquire(); })
          phi.acquire();
        struct SurfaceAccumulation
        {
            SurfaceState state;
            Real squaredResidual = 0;
        };
        std::vector<SurfaceAccumulation> facetStates(interfaceFacets.size());
#ifdef RODIN_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (Index i = 0; i < static_cast<Index>(interfaceFacets.size()); ++i)
        {
          auto& accumulation = facetStates[static_cast<std::size_t>(i)];
          SurfaceState& facetState = accumulation.state;
          DeformationMap deformation(u, locator);
          const Index facetIdx = interfaceFacets[static_cast<std::size_t>(i)];
          const auto face = mesh.getFace(facetIdx);
          const auto& qf = getQuadrature(*face, fes);
          const auto& quad = face->getQuadrature(qf);
          for (std::size_t q = 0; q < quad.getSize(); ++q)
          {
            const auto& src = quad.getPoint(q);
            const Variational::IntegrationPoint ip(src, &qf, q);
            const Real w = qf.getWeight(q) * src.getDistortion();
            const Real r = phi.getValue(deformation.getMovedPoint(ip));
            facetState.totalLen += w;
            facetState.energy += w * normalization * loss.getValue(r);
            accumulation.squaredResidual += w * r * r;
            facetState.residualSup = std::max(facetState.residualSup, std::abs(r));
          }
        }

        SurfaceState state;
        Real squared = 0;
        // Combine the per-facet positive quadrature sums before normalizing the
        // Residual over the complete interface.
        for (const auto& accumulation : facetStates)
        {
          const SurfaceState& facetState = accumulation.state;
          state.energy += facetState.energy;
          state.totalLen += facetState.totalLen;
          state.residualSup = std::max(state.residualSup, facetState.residualSup);
          squared += accumulation.squaredResidual;
        }
        state.residualRMS = state.totalLen > Real(0)
          ? std::sqrt(std::max(Real(0), squared) / state.totalLen)
          : std::numeric_limits<Real>::infinity();
        if (!(state.totalLen > Real(0)))
          state.residualSup = std::numeric_limits<Real>::infinity();
        return state;
      }

      /**
       * @brief Evaluates position and normal errors on the fitted interface.
       *
       * At a mapped interface point @f$y=T_h(x)@f$, the normalized residual
       * @f$\phi(y)/|\nabla\phi(y)|@f$ is the first-order signed distance to the
       * target zero set. The fitted normal is obtained from the tangents mapped
       * by @f$I+\nabla u_h@f$; consequently, curved and higher-order displacement
       * fields are measured without replacing their image by an affine facet.
       * Normal orientation is ignored because interior-facet orientation is not
       * intrinsic to the interface skeleton.
       */
      /**
       * @brief Computes geometric discrepancies on the complete fitted interface.
       * @param mesh Fixed reference mesh.
       * @param fes Displacement finite element space defining quadrature orders.
       * @param u Displacement whose deformed interface or geometry is evaluated.
       * @param phi Target level-set function.
       * @param grad Target level-set gradient.
       * @param interfaceFacets Indices of the classified interface facets.
       * @param dimension Spatial dimension.
       * @param locator Target-evaluation point locator.
       * @returns Sampled RMS and maximum geometric discrepancy and normal mismatch.
       */
      template <class Mesh, class FES, class PhiType, class GradType, class LocatorType>
      InterfaceGeometryState getInterfaceGeometryState(const Mesh& mesh, const FES& fes,
        const Displacement& u, const PhiType& phi, const GradType& grad,
        const std::vector<Index>& interfaceFacets, std::size_t dimension,
        const LocatorType& locator) const
      {
        if constexpr (requires { u.acquire(); })
          u.acquire();
        if constexpr (requires { phi.acquire(); })
          phi.acquire();
        if constexpr (requires { grad.acquire(); })
          grad.acquire();

        struct Accumulation
        {
            Real measure = 0;
            Real squaredDistance = 0;
            Real maximumDistance = 0;
            Real squaredNormal = 0;
            std::size_t invalidCount = 0;
        };
        std::vector<Accumulation> facetStates(interfaceFacets.size());
        constexpr Real gradientFloor = Real(1e-14);

#ifdef RODIN_USE_OPENMP
#pragma omp parallel
#endif
        {
          DeformationMap deformation(u, locator);
#ifdef RODIN_USE_OPENMP
#pragma omp for schedule(static)
#endif
          for (Index i = 0; i < static_cast<Index>(interfaceFacets.size()); ++i)
          {
            Accumulation& state = facetStates[static_cast<std::size_t>(i)];
            const Index facetIndex = interfaceFacets[static_cast<std::size_t>(i)];
            const auto face = mesh.getFace(facetIndex);
            const auto& fe = fes.getFiniteElement(face->getDimension(), facetIndex);
            const auto& qf = getGeometricValidationQuadrature(*face, fes);
            const auto& quadrature = face->getQuadrature(qf);
            // Endpoints/vertices contribute to the maximum, not the surface quadrature integral.
            const Geometry::Polytope::Traits traits(face->getGeometry());
            for (std::size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
            {
              const Geometry::Point point(*face, traits.getVertex(vertex));
              const auto& moved =
                deformation.getMovedPoint(Variational::IntegrationPoint(point));
              const Real norm = grad.getValue(moved).norm();
              const Real distance = std::abs(phi.getValue(moved)) / norm;
              if (!(norm > gradientFloor) || !std::isfinite(norm) ||
                !std::isfinite(distance))
                ++state.invalidCount;
              else
                state.maximumDistance = std::max(state.maximumDistance, distance);
            }
            for (std::size_t q = 0; q < quadrature.getSize(); ++q)
            {
              const auto& point = quadrature.getPoint(q);
              const Variational::IntegrationPoint ip(point, &qf, q);
              const auto& J = point.getJacobian();
              SpatialMat mappedJacobian = J;
              const auto& rc = point.getReferenceCoordinates();
              for (std::size_t local = 0; local < fe.getCount(); ++local)
              {
                const Index dof =
                  fes.getGlobalIndex({face->getDimension(), facetIndex}, local);
                for (std::size_t axis = 0; axis < face->getDimension(); ++axis)
                {
                  for (std::size_t component = 0; component < dimension; ++component)
                  {
                    mappedJacobian(static_cast<Eigen::Index>(component),
                      static_cast<Eigen::Index>(axis)) += u[dof] *
                      fe.getBasis(local).template getDerivative<1>(component, axis)(rc);
                  }
                }
              }

              SpatialVec fittedNormal = SpatialVec::Zero(dimension);
              Real mappedMeasure = 0;
              if (dimension == 1)
              {
                mappedMeasure = Real(1);
                fittedNormal(0) = Real(1);
              }
              else if (dimension == 2)
              {
                SpatialVec tangent = SpatialVec::Zero(2);
                tangent(0) = mappedJacobian(0, 0);
                tangent(1) = mappedJacobian(1, 0);
                mappedMeasure = tangent.norm();
                fittedNormal(0) = tangent(1);
                fittedNormal(1) = -tangent(0);
              }
              else
              {
                SpatialVec a = SpatialVec::Zero(3);
                SpatialVec b = SpatialVec::Zero(3);
                for (std::size_t r = 0; r < 3; ++r)
                {
                  a(static_cast<Eigen::Index>(r)) = mappedJacobian(r, 0);
                  b(static_cast<Eigen::Index>(r)) = mappedJacobian(r, 1);
                }
                fittedNormal = a.cross(b);
                mappedMeasure = fittedNormal.norm();
              }

              const Geometry::Point movedPoint = deformation.getMovedPoint(ip);
              const SpatialVec targetGradient = grad.getValue(movedPoint);
              const Real gradientNorm = targetGradient.norm();
              const Real fittedNormalNorm = fittedNormal.norm();
              if (!(mappedMeasure > gradientFloor && gradientNorm > gradientFloor &&
                    fittedNormalNorm > gradientFloor) ||
                !std::isfinite(mappedMeasure) || !std::isfinite(gradientNorm) ||
                !std::isfinite(fittedNormalNorm))
              {
                ++state.invalidCount;
                continue;
              }

              fittedNormal /= fittedNormalNorm;
              const SpatialVec targetNormal = targetGradient / gradientNorm;
              const Real distance = std::abs(phi.getValue(movedPoint)) / gradientNorm;
              if (!std::isfinite(distance))
              {
                ++state.invalidCount;
                continue;
              }
              const Real normalErrorSquared = Real(2) -
                Real(2) *
                  std::clamp(std::abs(fittedNormal.dot(targetNormal)), Real(0), Real(1));
              const Real weight = qf.getWeight(q) * mappedMeasure;
              state.measure += weight;
              state.squaredDistance += weight * distance * distance;
              state.maximumDistance = std::max(state.maximumDistance, distance);
              state.squaredNormal += weight * normalErrorSquared;
            }
          }
        }

        Accumulation total;
        for (const Accumulation& state : facetStates)
        {
          total.measure += state.measure;
          total.squaredDistance += state.squaredDistance;
          total.maximumDistance = std::max(total.maximumDistance, state.maximumDistance);
          total.squaredNormal += state.squaredNormal;
          total.invalidCount += state.invalidCount;
        }
        if (!(total.measure > Real(0)) || total.invalidCount > 0)
          return {};
        return {std::sqrt(std::max(Real(0), total.squaredDistance) / total.measure),
          total.maximumDistance,
          std::sqrt(std::max(Real(0), total.squaredNormal) / total.measure)};
      }

      /**
       * @brief Normal-jump statistics of the discrete interface.
       *
       * Facets are adjacent when they share a vertex in two dimensions or an
       * edge in three, and each adjacent pair contributes the angle between its
       * normals weighted by the mean of the two facet measures. Normals are
       * compared up to sign, since facet orientation is not consistent across
       * the interface.
       */
      /**
       * @brief Computes jumps between averaged adjacent interface normals.
       * @param mesh Fixed reference mesh.
       * @param fes Displacement finite element space defining quadrature orders.
       * @param interfaceFacets Indices of the classified interface facets.
       * @param dimension Spatial dimension.
       * @returns Measure-weighted RMS and maximum normal-angle jumps.
       */
      template <class Mesh, class FES>
      NormalJump getNormalJump(const Mesh& mesh, const FES& fes,
        const std::vector<Index>& interfaceFacets, std::size_t dimension) const
      {
        // Averaged unit normal and measure of each facet. The normal is the
        // quadrature average of the pointwise unit normals, so a curved facet
        // contributes its mean orientation; a facet degenerate throughout
        // falls back to a fixed direction rather than a null vector.
        constexpr Real degenerate = Real(1e-30);
        std::vector<Math::SpatialVector<Real>> normals(interfaceFacets.size());
        std::vector<Real> measures(interfaceFacets.size(), Real(0));

        for (std::size_t i = 0; i < interfaceFacets.size(); ++i)
        {
          const auto face = mesh.getFace(interfaceFacets[i]);
          const auto& qf = getQuadrature(*face, fes);
          const auto& quad = face->getQuadrature(qf);

          Math::SpatialVector<Real> n = Math::SpatialVector<Real>::Zero(dimension);
          Real measure = 0;
          for (std::size_t q = 0; q < quad.getSize(); ++q)
          {
            const auto& pt = quad.getPoint(q);
            const Real w = qf.getWeight(q) * pt.getDistortion();
            const auto& J = pt.getJacobian();
            // The facet normal: in 2D the tangent rotated a quarter turn, in
            // 3D the cross product of the two tangents.
            Math::SpatialVector<Real> nq;
            if (dimension == 2)
            {
              nq.resize(2);
              nq(0) = J(1, 0);
              nq(1) = -J(0, 0);
            }
            else
            {
              Math::SpatialVector<Real> a = Math::SpatialVector<Real>::Zero(dimension);
              Math::SpatialVector<Real> b = Math::SpatialVector<Real>::Zero(dimension);
              for (std::size_t r = 0; r < dimension; ++r)
              {
                a(static_cast<Eigen::Index>(r)) = J(r, 0);
                b(static_cast<Eigen::Index>(r)) = J(r, 1);
              }
              nq = a.cross(b);
            }
            const Real norm = nq.norm();
            if (norm > degenerate)
            {
              n += w * (nq / norm);
              measure += w;
            }
          }

          if (n.norm() <= degenerate)
          {
            n.setZero();
            n(0) = Real(1);
          }
          else
          {
            n /= n.norm();
          }
          measures[i] = measure;
          normals[i] = n;
        }

        const auto edgeKey = [](Index a, Index b) -> std::uint64_t {
          const auto lo = static_cast<std::uint64_t>(std::min(a, b));
          const auto hi = static_cast<std::uint64_t>(std::max(a, b));
          return (lo << 32) ^ hi;
        };

        // Facets sharing an incidence, keyed by the shared entity.
        UnorderedMap<std::uint64_t, std::vector<std::size_t>> incident;
        if (dimension == 2)
        {
          incident.reserve(2 * interfaceFacets.size());
          for (std::size_t i = 0; i < interfaceFacets.size(); ++i)
          {
            for (const Index v : mesh.getFace(interfaceFacets[i])->getVertices())
              incident[static_cast<std::uint64_t>(v)].push_back(i);
          }
        }
        else if (dimension == 3)
        {
          incident.reserve(3 * interfaceFacets.size());
          for (std::size_t i = 0; i < interfaceFacets.size(); ++i)
          {
            const auto& vv = mesh.getFace(interfaceFacets[i])->getVertices();
            incident[edgeKey(vv[0], vv[1])].push_back(i);
            incident[edgeKey(vv[1], vv[2])].push_back(i);
            incident[edgeKey(vv[2], vv[0])].push_back(i);
          }
        }

        NormalJump jump;
        Real squared = 0;
        Real weight = 0;
        for (const auto& [key, faces] : incident)
        {
          (void)key;
          for (std::size_t a = 0; a < faces.size(); ++a)
          {
            for (std::size_t b = a + 1; b < faces.size(); ++b)
            {
              const auto ia = faces[a], ib = faces[b];
              const Real dot = std::abs(normals[ia].dot(normals[ib]));
              const Real theta = std::acos(std::max(Real(-1), std::min(Real(1), dot)));
              const Real w = Real(0.5) * (measures[ia] + measures[ib]);
              squared += w * theta * theta;
              weight += w;
              jump.max = std::max(jump.max, theta);
            }
          }
        }
        jump.rms = weight > Real(0) ? std::sqrt(squared / weight) : Real(0);
        return jump;
      }

      /**
       * @brief Robust scale @f$\sigma@f$ of the interface residual.
       *
       * The loss weight requires a scale separating the residuals to be fitted
       * from those treated as outliers. It is the 90th percentile of the
       * undeformed residual, floored at @f$3hG_\phi@f$, where @f$G_\phi@f$ is
       * the maximum sampled target-gradient norm. The percentile adapts to the
       * initial mismatch, while the floor converts the geometric resolution to
       * level-set units.
       *
       * Fixed once per frame rather than re-estimated per iteration, which
       * would let the scale chase its own progress and never reject anything.
       * A positive @ref Parameters::Model::robustScale overrides the automatic
       * mesh-dependent selection.
       */
      struct RobustScale
      {
          /// @brief Fixed robust-loss scale in level-set units.
          Real sigma = 0;
          /// @brief Maximum sampled target-gradient norm.
          Real gradientScale = 0;
      };

      /**
       * @brief Estimates a fixed robust scale from the initial interface.
       * @param mesh Fixed reference mesh.
       * @param fes Displacement finite element space defining quadrature orders.
       * @param phi Target level-set function.
       * @param grad Target level-set gradient.
       * @param interfaceFacets Indices of the classified interface facets.
       * @param h Reference mesh spacing.
       * @returns The robust residual scale and sampled gradient normalization.
       */
      template <class Mesh, class FES, class PhiType, class GradType>
      RobustScale getRobustScale(const Mesh& mesh, const FES& fes, const PhiType& phi,
        const GradType& grad, const std::vector<Index>& interfaceFacets, Real h) const
      {
        std::vector<Real> residuals;
        Real gradientScale = 0;
        for (const Index facetIdx : interfaceFacets)
        {
          const auto face = mesh.getFace(facetIdx);
          const auto& qf = getQuadrature(*face, fes);
          const auto& quad = face->getQuadrature(qf);
          for (std::size_t q = 0; q < quad.getSize(); ++q)
          {
            const auto& point = quad.getPoint(q);
            residuals.push_back(std::abs(phi.getValue(point)));
            gradientScale = std::max(gradientScale, grad.getValue(point).norm());
          }
        }

        Real sigma = m_parameters.model.robustScale;
        if (!(sigma > Real(0)))
        {
          sigma = RobustScaleMeshFloor * h * gradientScale;
          const std::size_t k90 = static_cast<std::size_t>(
            RobustScaleQuantile * static_cast<Real>(residuals.size() - 1));
          std::nth_element(residuals.begin(), residuals.begin() + k90, residuals.end());
          sigma = std::max(sigma, residuals[k90]);
        }
        sigma = std::max(sigma, std::sqrt(std::numeric_limits<Real>::min()));
        return {sigma, gradientScale};
      }

      void assembleBoundaryConditions()
      {
        m_fixedIncrementDOFs.clear();
        for (auto& condition : m_boundaryConditions)
        {
          condition.assemble();
          const auto* values =
            std::get_if<typename Variational::DirichletBCBase<Real>::ValueDOFs>(
              &condition.getDOFs());
          if (!values)
            Alert::Exception() << "SWIFT supports homogeneous value conditions only."
                               << Alert::Raise;
          for (const auto& [dof, value] : *values)
          {
            if (value != Real(0))
              Alert::Exception()
                << "SWIFT boundary increments must be zero." << Alert::Raise;
            m_fixedIncrementDOFs[dof] = Real(0);
          }
        }
      }

      /// @brief Assemble one symmetric operator for the model, factors and residuals.
      void assembleStep()
      {
        m_metricProblem.assemble();
        using OperatorType =
          typename FormLanguage::Traits<LinearSystemType>::OperatorType;
        if constexpr (std::is_same_v<OperatorType, Math::SparseMatrix<Real>>)
        {
          auto& matrix = m_metricProblem.getLinearSystem().getOperator();
          matrix =
            (Real(0.5) * (matrix + Math::SparseMatrix<Real>(matrix.transpose()))).eval();
          if (!m_fixedIncrementDOFs.empty())
            m_metricProblem.getLinearSystem().eliminate(m_fixedIncrementDOFs);
          matrix.makeCompressed();
        }
      }

      /**
       * @brief Solves the predictor and transfers it to a displacement field.
       * The retained linear adapter also serves the native inner Newton solve.
       */
      /**
       * @brief Solves the assembled centered metric for the unconstrained predictor.
       * @param out Returned predictor displacement, updated only on success.
       * @param report Destination for linear-solver diagnostics.
       * @returns Whether the centered physical residual met the linear tolerance.
       */
      bool solvePredictor(Displacement& out, Report& report)
      {
        m_linearSolver.setParameters(m_parameters).setReport(&report);
        auto& system = m_metricProblem.getLinearSystem();
        m_linearSolver.solve(system);
        if (!m_linearSolver.success())
          return false;
        m_duStep.getSolution().setData(system.getSolution());
        out = m_duStep.getSolution();
        return true;
      }

      /// Dimensionless sufficient-decrease coefficient for the frozen inner merit.
      static constexpr Real InnerArmijoCoefficient = Real(1e-4);
      /// Safeguard on trials per inner correction, independent of the Newton cap.
      static constexpr size_t InnerMaxBacktracks = 32;
      /// Dimensionless reduction shared by the inner and outer Armijo searches.
      static constexpr Real BacktrackingReduction = Real(0.5);
      /// Relative floating-point allowance in the quadratic inner merit.
      static constexpr Real MeritRoundoffFactor = Real(64);
      /// Dimensionless tolerance for reporting an undamped inner correction.
      static constexpr Real FullStepTolerance = Real(1e-12);
      /// Minimum robust scale in units of the reference mesh spacing.
      static constexpr Real RobustScaleMeshFloor = Real(3);
      /// Residual quantile used by automatic robust-scale selection.
      static constexpr Real RobustScaleQuantile = Real(0.9);

      Displacement* m_u;
      Identifiable::UUID m_trialUUID;
      TrialFunctionType m_duStep;
      TestFunctionType m_vStep;
      ProblemType m_metricProblem;
      HingeProblem m_hingeProblem;
      LinearSolver m_linearSolver;
      BilinearFormType m_distributionForm;
      Math::Matrix<Real> m_centering;
      /**
       * @brief Fitting metric and fitting force at the outer displacement.
       * Both depend on the outer displacement only, so they are assembled once
       * per nonlinear iteration and reused by every hinge correction.
       */
      BilinearFormType m_fittingMetric;
      LinearFormType m_fittingForce;
      BilinearFormType m_additionalMetric;
      Variational::EssentialBoundary<Real> m_boundaryConditions;
      IndexMap<Real> m_fixedIncrementDOFs;
      Parameters m_parameters;
      Report m_report;
      Optional<Monitor> m_monitor;
  };
}

#endif
