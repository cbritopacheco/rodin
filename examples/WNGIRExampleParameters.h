/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXAMPLES_WNGIREXAMPLEPARAMETERS_H
#define RODIN_EXAMPLES_WNGIREXAMPLEPARAMETERS_H

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

#include <Rodin/Adaptation/WNGIR/Parameters.h>
#include <Rodin/Adaptation/WNGIR/Report.h>
#include <Rodin/MMG/MeshOptimizer.h>

namespace Rodin::Examples
{
  /// @brief Lossless canonical responses, separate from the human-readable summary.
  inline void printWNGIRResponses(const Adaptation::WNGIRReport& report)
  {
    const auto flags = std::cout.flags();
    const auto precision = std::cout.precision();
    std::cout << std::scientific
              << std::setprecision(std::numeric_limits<Real>::max_digits10)
              << "    wngir responses: energy=" << report.energy
              << " geom_sup=" << report.geometricSup
              << " geom_sup_target=" << report.geometricSupTarget
              << " target_hit=" << report.geometricTargetReached
              << " quality_ok=" << report.qualityBudgetSatisfied
              << " outer=" << report.iterations
              << " inner_total=" << report.innerIterations
              << " inner_max=" << report.maxInnerIterations
              << " inner_last=" << report.lastInnerIterations
              << " inner_converged=" << report.innerConverged
              << " inner_residual=" << report.innerResidual
              << " inner_relative_residual=" << report.innerRelativeResidual
              << " inner_residual_tolerance=" << report.innerResidualTolerance
              << " min_j=" << report.minJ << " max_qrel=" << report.maxQRel
              << " exit=" << report.getReasonString() << '\n';
    std::cout.flags(flags);
    std::cout.precision(precision);
  }

  inline std::string wngirOutput(const std::string& name)
  {
    std::filesystem::create_directories("wngir");
    return "wngir/" + name;
  }

  struct WNGIRExampleDefaults
  {
      std::size_t maxIterations =
        Adaptation::WNGIRParameters{}.convergence.iterations.outer;
      std::size_t quadratureOrder = 0;
      Real kappaF = Adaptation::WNGIRParameters{}.model.fit;
      Real kappaD = Adaptation::WNGIRParameters{}.model.distribution;
      Real kappaJ = Adaptation::WNGIRParameters{}.model.kappaJ;
      Real kappaQ = Adaptation::WNGIRParameters{}.model.kappaQ;
  };

  inline bool findOption(
    int argc, char** argv, const std::string& name, std::string* value)
  {
    const std::string prefix = "--" + name + "=";
    const std::string flag = "--" + name;
    for (int i = 1; i < argc; ++i)
    {
      const std::string arg(argv[i]);
      if (arg.rfind(prefix, 0) == 0)
      {
        if (value)
          *value = arg.substr(prefix.size());
        return true;
      }
      if (arg == flag)
      {
        if (value)
        {
          if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0)
            *value = argv[i + 1];
          else
            value->clear();
        }
        return true;
      }
    }
    return false;
  }

  inline Real realOption(int argc, char** argv, const std::string& name, Real fallback)
  {
    std::string value;
    if (!findOption(argc, argv, name, &value) || value.empty())
      return fallback;
    return static_cast<Real>(std::atof(value.c_str()));
  }

  inline Real realOption(int argc, char** argv, const std::string& name,
    const std::string& legacyName, Real fallback)
  {
    const Real legacy = realOption(argc, argv, legacyName, fallback);
    return realOption(argc, argv, name, legacy);
  }

  inline std::size_t sizeOption(
    int argc, char** argv, const std::string& name, std::size_t fallback)
  {
    std::string value;
    if (!findOption(argc, argv, name, &value) || value.empty())
      return fallback;
    return static_cast<std::size_t>(std::strtoull(value.c_str(), nullptr, 10));
  }

  inline bool boolOption(int argc, char** argv, const std::string& name, bool fallback)
  {
    std::string value;
    if (!findOption(argc, argv, name, &value))
      return fallback;
    if (value.empty())
      return true;
    return std::strtoull(value.c_str(), nullptr, 10) != 0;
  }

  inline std::string stringOption(
    int argc, char** argv, const std::string& name, std::string fallback)
  {
    std::string value;
    if (!findOption(argc, argv, name, &value) || value.empty())
      return fallback;
    return value;
  }

  /// Remesh the linear background before constructing spaces or curved maps.
  /// MMG lengths are physical: hmin = 0.1 h, hmax = h, hausd = 0.05 h.
  inline void remeshWNGIRBackground(
    Geometry::Mesh<Context::Local>& mesh, int argc, char** argv, Real h)
  {
    if (!boolOption(argc, argv, "initial-mmg", false))
      return;
    const Real hmin = realOption(
      argc, argv, "initial-mmg-hmin", realOption(argc, argv, "hmin", Real(0.1) * h));
    const Real hmax = realOption(argc, argv, "initial-mmg-hmax", h);
    const Real hausd = realOption(argc, argv, "initial-mmg-hausd", Real(0.05) * h);
    if (!(std::isfinite(hmin) && std::isfinite(hmax) && std::isfinite(hausd) &&
          hmin > 0 && hmax >= hmin && hausd > 0))
    {
      try
      {
        Alert::Exception() << "Initial MMG requires 0 < hmin <= hmax and hausd > 0."
                           << Alert::Raise;
      }
      catch (const Alert::Exception&)
      {
        std::exit(EXIT_FAILURE);
      }
    }
    std::cout << "  initial MMG remesh: h=" << h << " hmin=" << hmin << " hmax=" << hmax
              << " hausd=" << hausd << '\n';
    MMG::Mesh mmgMesh(std::move(mesh));
    MMG::Optimizer optimizer;
    // Preserve the corners and ridges of the computational box.
    optimizer.setHMin(hmin).setHMax(hmax).setHausdorff(hausd).setAngleDetection(true);
    optimizer.optimize(mmgMesh);
    mesh = std::move(static_cast<MMG::Mesh::Parent&>(mmgMesh));
  }

  inline Adaptation::WNGIRParameters makeWNGIRParameters(int argc, char** argv, Real h,
    Geometry::Attribute interfaceAttribute, const WNGIRExampleDefaults& defaults = {})
  {
    constexpr const char* options[] = {"wngir-kappa-f", "wngir-robust-scale",
      "wngir-kappa-j", "wngir-kappa-q", "wngir-jsafe", "wngir-qmax",
      "wngir-quality-guard", "wngir-kappa-d", "wngir-directional-newton",
      "wngir-directional-newton-max-step-h", "wngir-quality-witness",
      "wngir-direct-solver", "wngir-direct-threads", "wngir-geometric-sup-tol",
      "wngir-primal-barrier-iterations", "wngir-primal-barrier-relative-tol",
      "wngir-primal-barrier-absolute-tol", "wngir-stagnation-iterations", "wngir-mu-hat",
      "wngir-omega-min", "wngir-max-backtracks", "wngir-armijo", "wngir-jls",
      "wngir-energy-stag-tol", "wngir-step-tol", "wngir-step-h-tol", "wngir-steps",
      "wngir-cg-rtol", "wngir-cg-max-iters", "wngir-trace"};
    for (int i = 1; i < argc; ++i)
    {
      const std::string argument(argv[i]);
      if (!argument.starts_with("--wngir-"))
        continue;
      const auto name = argument.substr(2,
        argument.find('=') == std::string::npos ? std::string::npos
                                                : argument.find('=') - 2);
      if (std::none_of(std::begin(options), std::end(options),
            [&](const char* option) { return name == option; }))
        Alert::Exception() << "Unknown or removed WNGIR option: " << name << Alert::Raise;
    }
    Adaptation::WNGIRParameters p;
    p.model.h = h;

    p.model.fit = realOption(argc, argv, "wngir-kappa-f", defaults.kappaF);
    p.model.robustScale =
      realOption(argc, argv, "wngir-robust-scale", p.model.robustScale);

    p.model.kappaJ = realOption(argc, argv, "wngir-kappa-j", defaults.kappaJ);
    p.model.kappaQ = realOption(argc, argv, "wngir-kappa-q", defaults.kappaQ);
    p.model.jacobian = realOption(argc, argv, "wngir-jsafe", "j-safe", p.model.jacobian);
    p.model.distortion = realOption(argc, argv, "wngir-qmax", p.model.distortion);
    p.model.qualityGuard =
      realOption(argc, argv, "wngir-quality-guard", p.model.qualityGuard);
    p.model.distribution = realOption(argc, argv, "wngir-kappa-d", defaults.kappaD);
    p.globalization.directionalNewton = boolOption(
      argc, argv, "wngir-directional-newton", p.globalization.directionalNewton);
    p.globalization.maxStepOverH = realOption(
      argc, argv, "wngir-directional-newton-max-step-h", p.globalization.maxStepOverH);
    p.traceQualityWitness = boolOption(argc, argv, "wngir-quality-witness", false);
    const auto defaultSolver =
      p.linear.solver == Adaptation::WNGIRParameters::LinearSolver::MUMPS ? "mumps"
                                                                          : "sparse-lu";
    const auto linearSolver =
      stringOption(argc, argv, "wngir-direct-solver", defaultSolver);
    if (linearSolver == "mumps")
      p.linear.solver = Adaptation::WNGIRParameters::LinearSolver::MUMPS;
    else if (linearSolver == "sparse-lu")
      p.linear.solver = Adaptation::WNGIRParameters::LinearSolver::SparseLU;
    else if (linearSolver == "cg")
      p.linear.solver = Adaptation::WNGIRParameters::LinearSolver::CG;
    else
      Alert::Exception() << "Unknown WNGIR solver: " << linearSolver << Alert::Raise;
    p.linear.threads = sizeOption(argc, argv, "wngir-direct-threads", 0);
    p.convergence.tolerance.geometric =
      realOption(argc, argv, "wngir-geometric-sup-tol", 0);
    p.convergence.iterations.inner = sizeOption(
      argc, argv, "wngir-primal-barrier-iterations", p.convergence.iterations.inner);
    p.convergence.tolerance.innerRelative = realOption(argc, argv,
      "wngir-primal-barrier-relative-tol", p.convergence.tolerance.innerRelative);
    p.convergence.tolerance.innerAbsolute = realOption(argc, argv,
      "wngir-primal-barrier-absolute-tol", p.convergence.tolerance.innerAbsolute);
    p.convergence.iterations.stagnation = sizeOption(
      argc, argv, "wngir-stagnation-iterations", p.convergence.iterations.stagnation);
    p.model.hinge = realOption(argc, argv, "wngir-mu-hat", p.model.hinge);
    p.model.omegaMin = realOption(argc, argv, "wngir-omega-min", p.model.omegaMin);
    p.convergence.iterations.backtracks =
      sizeOption(argc, argv, "wngir-max-backtracks", p.convergence.iterations.backtracks);
    p.globalization.armijo =
      realOption(argc, argv, "wngir-armijo", p.globalization.armijo);

    p.globalization.jMin = realOption(argc, argv, "j-min", p.globalization.jMin);
    p.globalization.jFloor = realOption(
      argc, argv, "wngir-jls", "j-ls", std::max(p.globalization.jMin, p.model.jacobian));
    p.convergence.tolerance.energy =
      realOption(argc, argv, "wngir-energy-stag-tol", p.convergence.tolerance.energy);
    p.convergence.tolerance.step =
      realOption(argc, argv, "wngir-step-tol", p.convergence.tolerance.step);
    p.convergence.tolerance.stepOverH =
      realOption(argc, argv, "wngir-step-h-tol", p.convergence.tolerance.stepOverH);

    p.quadrature.order = sizeOption(argc, argv, "quad-order", defaults.quadratureOrder);
    p.quadrature.validation =
      sizeOption(argc, argv, "geometric-validation-order", p.quadrature.validation);
    p.convergence.iterations.outer =
      sizeOption(argc, argv, "wngir-steps", defaults.maxIterations);

    p.convergence.tolerance.linearRelative =
      realOption(argc, argv, "wngir-cg-rtol", p.convergence.tolerance.linearRelative);
    p.convergence.iterations.linear =
      sizeOption(argc, argv, "wngir-cg-max-iters", p.convergence.iterations.linear);
    p.interfaceAttribute = interfaceAttribute;
    p.trace =
      boolOption(argc, argv, "trace", boolOption(argc, argv, "wngir-trace", false));
    return p;
  }
}

#endif
