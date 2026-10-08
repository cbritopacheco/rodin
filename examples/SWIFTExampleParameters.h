/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_EXAMPLES_SWIFTEXAMPLEPARAMETERS_H
#define RODIN_EXAMPLES_SWIFTEXAMPLEPARAMETERS_H

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

#include <Rodin/Adaptation/SWIFT/Parameters.h>
#include <Rodin/Adaptation/SWIFT/Report.h>
#include <Rodin/MMG/MeshOptimizer.h>

namespace Rodin::Examples
{
  /// @brief Lossless canonical responses, separate from the human-readable summary.
  inline void printSWIFTResponses(const Adaptation::SWIFT::Report& report)
  {
    const auto flags = std::cout.flags();
    const auto precision = std::cout.precision();
    std::cout << std::scientific
              << std::setprecision(std::numeric_limits<Real>::max_digits10)
              << "    swift responses: energy=" << report.energy
              << " geom_sup=" << report.geometricSup
              << " geom_c=" << report.geometricConstant
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

  inline std::string swiftOutput(const std::string& name)
  {
    std::filesystem::create_directories("swift");
    return "swift/" + name;
  }

  struct SWIFTExampleDefaults
  {
      std::size_t maxIterations =
        Adaptation::SWIFT::Parameters{}.convergence.iterations.outer;
      std::size_t quadratureOrder = 0;
      Real fit = Adaptation::SWIFT::Parameters{}.model.fit;
      Real deviatoric = Adaptation::SWIFT::Parameters{}.model.distribution.deviatoric;
      Real divergence = Adaptation::SWIFT::Parameters{}.model.distribution.divergence;
      Real jacobianWeight = Adaptation::SWIFT::Parameters{}.model.jacobianWeight;
      Real distortionWeight = Adaptation::SWIFT::Parameters{}.model.distortionWeight;
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
  inline void remeshSWIFTBackground(
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

  inline Adaptation::SWIFT::Parameters makeSWIFTParameters(int argc, char** argv, Real h,
    Geometry::Attribute interfaceAttribute, const SWIFTExampleDefaults& defaults = {})
  {
    constexpr const char* options[] = {"swift-fit", "swift-robust-scale",
      "swift-jacobian-weight", "swift-distortion-weight", "swift-jacobian", "swift-distortion",
      "swift-quality-guard", "swift-distribution-deviatoric",
      "swift-distribution-divergence", "swift-directional-newton",
      "swift-max-step-over-h", "swift-quality-witness",
      "swift-linear-solver", "swift-linear-threads",
      "swift-geometric-tolerance", "swift-inner-iterations",
      "swift-inner-relative-tolerance", "swift-inner-absolute-tolerance",
      "swift-stagnation-iterations", "swift-hinge",
      "swift-backtracks", "swift-armijo", "swift-energy-tolerance",
      "swift-step-tolerance", "swift-step-over-h-tolerance", "swift-outer-iterations", "swift-linear-relative-tolerance",
      "swift-linear-iterations", "swift-trace"};
    for (int i = 1; i < argc; ++i)
    {
      const std::string argument(argv[i]);
      if (argument.starts_with("--wngir-"))
        Alert::Exception() << "WNGIR options have been renamed to --swift-*: "
                           << argument << Alert::Raise;
      if (!argument.starts_with("--swift-"))
        continue;
      const auto name = argument.substr(2,
        argument.find('=') == std::string::npos ? std::string::npos
                                                : argument.find('=') - 2);
      if (std::none_of(std::begin(options), std::end(options),
            [&](const char* option) { return name == option; }))
        Alert::Exception() << "Unknown or removed SWIFT option: " << name << Alert::Raise;
    }
    Adaptation::SWIFT::Parameters p;
    p.model.h = h;

    p.model.fit = realOption(argc, argv, "swift-fit", defaults.fit);
    p.model.robustScale =
      realOption(argc, argv, "swift-robust-scale", p.model.robustScale);

    p.model.jacobianWeight =
      realOption(argc, argv, "swift-jacobian-weight", defaults.jacobianWeight);
    p.model.distortionWeight =
      realOption(argc, argv, "swift-distortion-weight", defaults.distortionWeight);
    p.model.jacobian = realOption(argc, argv, "swift-jacobian", p.model.jacobian);
    p.model.distortion = realOption(argc, argv, "swift-distortion", p.model.distortion);
    p.model.qualityGuard =
      realOption(argc, argv, "swift-quality-guard", p.model.qualityGuard);
    p.model.distribution.deviatoric =
      realOption(argc, argv, "swift-distribution-deviatoric", defaults.deviatoric);
    p.model.distribution.divergence =
      realOption(argc, argv, "swift-distribution-divergence", defaults.divergence);
    p.globalization.directionalNewton = boolOption(
      argc, argv, "swift-directional-newton", p.globalization.directionalNewton);
    p.globalization.maxStepOverH = realOption(
      argc, argv, "swift-max-step-over-h", p.globalization.maxStepOverH);
    p.traceQualityWitness = boolOption(argc, argv, "swift-quality-witness", false);
    const auto defaultSolver =
      p.linear.solver == Adaptation::SWIFT::Parameters::LinearSolver::MUMPS ? "mumps"
                                                                          : "sparse-lu";
    const auto linearSolver =
      stringOption(argc, argv, "swift-linear-solver", defaultSolver);
    if (linearSolver == "mumps")
      p.linear.solver = Adaptation::SWIFT::Parameters::LinearSolver::MUMPS;
    else if (linearSolver == "sparse-lu")
      p.linear.solver = Adaptation::SWIFT::Parameters::LinearSolver::SparseLU;
    else if (linearSolver == "cg")
      p.linear.solver = Adaptation::SWIFT::Parameters::LinearSolver::CG;
    else
      Alert::Exception() << "Unknown SWIFT solver: " << linearSolver << Alert::Raise;
    p.linear.threads = sizeOption(argc, argv, "swift-linear-threads", 0);
    p.convergence.tolerance.geometric =
      realOption(argc, argv, "swift-geometric-tolerance", 0);
    p.convergence.iterations.inner = sizeOption(
      argc, argv, "swift-inner-iterations", p.convergence.iterations.inner);
    p.convergence.tolerance.innerRelative = realOption(argc, argv,
      "swift-inner-relative-tolerance", p.convergence.tolerance.innerRelative);
    p.convergence.tolerance.innerAbsolute = realOption(argc, argv,
      "swift-inner-absolute-tolerance", p.convergence.tolerance.innerAbsolute);
    p.convergence.iterations.stagnation = sizeOption(
      argc, argv, "swift-stagnation-iterations", p.convergence.iterations.stagnation);
    p.model.hinge = realOption(argc, argv, "swift-hinge", p.model.hinge);
    p.convergence.iterations.backtracks =
      sizeOption(argc, argv, "swift-backtracks", p.convergence.iterations.backtracks);
    p.globalization.armijo =
      realOption(argc, argv, "swift-armijo", p.globalization.armijo);

    p.convergence.tolerance.energy =
      realOption(argc, argv, "swift-energy-tolerance", p.convergence.tolerance.energy);
    p.convergence.tolerance.step =
      realOption(argc, argv, "swift-step-tolerance", p.convergence.tolerance.step);
    p.convergence.tolerance.stepOverH =
      realOption(argc, argv, "swift-step-over-h-tolerance", p.convergence.tolerance.stepOverH);

    p.quadrature.order = sizeOption(argc, argv, "quad-order", defaults.quadratureOrder);
    p.quadrature.surface = sizeOption(argc, argv, "surface-quadrature-order", 0);
    p.quadrature.volume = sizeOption(argc, argv, "volume-quadrature-order", 0);
    p.quadrature.quality = sizeOption(argc, argv, "quality-validation-order", 0);
    p.quadrature.validation =
      sizeOption(argc, argv, "geometric-validation-order", p.quadrature.validation);
    p.convergence.iterations.outer =
      sizeOption(argc, argv, "swift-outer-iterations", defaults.maxIterations);

    p.convergence.tolerance.linearRelative = realOption(
      argc, argv, "swift-linear-relative-tolerance", p.convergence.tolerance.linearRelative);
    p.convergence.iterations.linear =
      sizeOption(argc, argv, "swift-linear-iterations", p.convergence.iterations.linear);
    p.interfaceAttribute = interfaceAttribute;
    p.trace =
      boolOption(argc, argv, "trace", boolOption(argc, argv, "swift-trace", false));
    return p;
  }
}

#endif
