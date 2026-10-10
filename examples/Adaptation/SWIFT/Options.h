/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_SWIFT_OPTIONS_H
#define RODIN_SWIFT_OPTIONS_H

#include <charconv>
#include <cmath>
#include <iostream>
#include <locale>
#include <sstream>
#include <string>
#include <type_traits>
#include <variant>

#include <Rodin/Adaptation/SWIFT/Parameters.h>
#include <Rodin/Alert.h>

namespace Rodin::Examples
{
  /// Named reconstruction controls shared by every displacement degree.
  struct ReconstructionOptions
  {
      size_t n = 16, dimension = 2;
      Real lobes = 0, amplitude = Real(0.05), radius = Real(0.25), phase = 0;
      Real cx = Real(0.5), cy = Real(0.5), cz = Real(0.5);
      std::string output, solver;
      bool help = false;
      Adaptation::SWIFT::Parameters parameters;

      /**
       * @brief Parses the same named options for every reconstruction degree.
       * @param argc Number of command-line arguments.
       * @param argv Command-line arguments, including the executable name.
       */
      ReconstructionOptions(int argc, char** argv)
      {
        using Value = std::variant<Real*, size_t*, bool*, std::string*>;
        struct Option
        {
            Value value;
            const char* description;
        };
        auto& p = parameters;
        solver = p.linear.solver == Adaptation::SWIFT::Parameters::LinearSolver::MUMPS
          ? "mumps"
          : "sparse-lu";
        const FlatMap<std::string, Option> options = {
          {"n", {&n, "Grid points per axis (at least 4)"}},
          {"dimension", {&dimension, "Spatial dimension: 2 or 3"}},
          {"lobes", {&lobes, "Angular frequency; zero selects a circle/sphere"}},
          {"amp", {&amplitude, "Radial perturbation amplitude"}},
          {"R0", {&radius, "Base target radius"}},
          {"phase", {&phase, "Target angular phase in radians"}},
          {"cx", {&cx, "Target center, first coordinate"}},
          {"cy", {&cy, "Target center, second coordinate"}},
          {"cz", {&cz, "Target center, third coordinate"}},
          {"output", {&output, "XDMF output stem (no extension); automatic if empty"}},
          {"help", {&help, "Print options and exit"}},
          {"swift-fit", {&p.model.fit, "Fitting stiffness"}},
          {"swift-distribution-deviatoric",
            {&p.model.distribution.deviatoric, "Centered deviatoric-strain stiffness"}},
          {"swift-distribution-divergence",
            {&p.model.distribution.divergence, "Centered divergence stiffness"}},
          {"swift-hinge", {&p.model.hinge, "Hinge/model-decrease ratio"}},
          {"swift-jacobian", {&p.model.jacobian, "Relative Jacobian floor"}},
          {"swift-distortion", {&p.model.distortion, "Relative distortion budget"}},
          {"swift-quality-guard", {&p.model.qualityGuard, "Hinge activation guard"}},
          {"swift-jacobian-weight",
            {&p.model.jacobianWeight, "Jacobian hinge row weight"}},
          {"swift-distortion-weight",
            {&p.model.distortionWeight, "Distortion hinge row weight"}},
          {"swift-robust-scale",
            {&p.model.robustScale, "Welsch scale; zero selects automatic"}},
          {"swift-directional-newton",
            {&p.globalization.directionalNewton, "Enable directional Newton (0 or 1)"}},
          {"swift-max-step-over-h",
            {&p.globalization.maxStepOverH,
              "Predictor motion/h cap; zero is unrestricted"}},
          {"swift-armijo", {&p.globalization.armijo, "Armijo sufficient decrease"}},
          {"swift-linear-solver", {&solver, "Linear backend: cg, sparse-lu or mumps"}},
          {"swift-linear-threads",
            {&p.linear.threads, "Backend threads; zero keeps its policy"}},
          {"swift-geometric-tolerance",
            {&p.convergence.tolerance.geometric,
              "Geometric target; zero selects h^(p+1)"}},
          {"swift-inner-relative-tolerance",
            {&p.convergence.tolerance.innerRelative,
              "Relative inner stationarity tolerance"}},
          {"swift-inner-absolute-tolerance",
            {&p.convergence.tolerance.innerAbsolute,
              "Absolute inner stationarity allowance"}},
          {"swift-linear-relative-tolerance",
            {&p.convergence.tolerance.linearRelative, "Linear residual tolerance"}},
          {"swift-energy-tolerance",
            {&p.convergence.tolerance.energy, "Relative energy stagnation threshold"}},
          {"swift-step-tolerance",
            {&p.convergence.tolerance.step, "Absolute accepted-motion threshold"}},
          {"swift-step-over-h-tolerance",
            {&p.convergence.tolerance.stepOverH, "Accepted-motion/h threshold"}},
          {"swift-outer-iterations",
            {&p.convergence.iterations.outer, "Outer fitting iteration cap"}},
          {"swift-inner-iterations",
            {&p.convergence.iterations.inner,
              "Inner Newton correction cap per outer step"}},
          {"swift-linear-iterations",
            {&p.convergence.iterations.linear, "CG iteration cap per linear solve"}},
          {"swift-backtracks",
            {&p.convergence.iterations.backtracks, "Outer trial-halving cap"}},
          {"swift-stagnation-iterations",
            {&p.convergence.iterations.stagnation,
              "Consecutive small changes before stagnation"}},
          {"quad-order",
            {&p.quadrature.order, "Common integration override; zero is automatic"}},
          {"surface-quadrature-order",
            {&p.quadrature.surface, "Surface integration override"}},
          {"volume-quadrature-order",
            {&p.quadrature.volume, "Volume integration override"}},
          {"quality-validation-order",
            {&p.quadrature.quality, "Independent quality sampling order"}},
          {"geometric-validation-order",
            {&p.quadrature.validation, "Independent geometric sampling order"}},
          {"swift-trace", {&p.trace, "Print iteration diagnostics (0 or 1)"}},
          {"trace", {&p.trace, "Alias of swift-trace"}},
          {"swift-quality-witness",
            {&p.traceQualityWitness,
              "Record limiting-cell quality diagnostics (0 or 1)"}}};
        for (int i = 1; i < argc; ++i)
        {
          const std::string argument(argv[i]);
          if (!argument.starts_with("--"))
            Alert::Exception() << "Expected a named option: " << argument << Alert::Raise;
          const auto separator = argument.find('=');
          const auto name = argument.substr(
            2, separator == std::string::npos ? std::string::npos : separator - 2);
          const auto option = options.find(name);
          if (option == options.end())
            Alert::Exception() << "Unknown reconstruction option: " << argument
                               << Alert::Raise;
          std::string value;
          if (separator != std::string::npos)
            value = argument.substr(separator + 1);
          else if (i + 1 < argc && !std::string(argv[i + 1]).starts_with("--"))
            value = argv[++i];
          else if (std::holds_alternative<bool*>(option->second.value))
            value = "1";
          else
            Alert::Exception() << "Missing value for --" << name << Alert::Raise;
          std::visit(
            [&]<class T>(T* destination) {
              if constexpr (std::is_same_v<T, std::string>)
                *destination = value;
              else if constexpr (std::is_same_v<T, bool>)
              {
                if (value != "0" && value != "1")
                  Alert::Exception() << "Expected 0 or 1 for --" << name << Alert::Raise;
                *destination = value == "1";
              }
              else if constexpr (std::is_floating_point_v<T>)
              {
                std::istringstream input(value);
                input.imbue(std::locale::classic());
                T parsed{};
                input >> std::noskipws >> parsed;
                if (!input || input.peek() != std::char_traits<char>::eof() ||
                  !std::isfinite(parsed))
                  Alert::Exception() << "Invalid numeric value for --" << name << ": "
                                     << value << Alert::Raise;
                *destination = parsed;
              }
              else
              {
                T parsed{};
                const auto result =
                  std::from_chars(value.data(), value.data() + value.size(), parsed);
                if (result.ec != std::errc{} ||
                  result.ptr != value.data() + value.size() ||
                  !std::isfinite(static_cast<Real>(parsed)))
                  Alert::Exception() << "Invalid numeric value for --" << name << ": "
                                     << value << Alert::Raise;
                *destination = parsed;
              }
            },
            option->second.value);
        }
        if (help)
        {
          std::cout << "Usage: " << argv[0] << " [--name=value | --name value]...\n";
          for (const auto& [name, option] : options)
          {
            std::cout << "  --" << name << " : " << option.description << " (value: ";
            std::visit([](const auto* value) { std::cout << *value; }, option.value);
            std::cout << ")\n";
          }
          return;
        }
        if (n < 4 || (dimension != 2 && dimension != 3))
          Alert::Exception() << "Expected n >= 4 and dimension 2 or 3." << Alert::Raise;
        if (lobes < 0 || (dimension == 2 && lobes != std::floor(lobes)) ||
          amplitude < 0 || !(radius > amplitude))
          Alert::Exception()
            << "Expected nonnegative lobes/amp, R0 > amp, and integer lobes in 2D."
            << Alert::Raise;
        if (solver == "cg")
          p.linear.solver = Adaptation::SWIFT::Parameters::LinearSolver::CG;
        else if (solver == "sparse-lu")
          p.linear.solver = Adaptation::SWIFT::Parameters::LinearSolver::SparseLU;
        else if (solver == "mumps")
        {
#ifdef RODIN_USE_MUMPS
          p.linear.solver = Adaptation::SWIFT::Parameters::LinearSolver::MUMPS;
#else
          Alert::Exception() << "MUMPS is not enabled in this build." << Alert::Raise;
#endif
        }
        else
          Alert::Exception() << "Unknown SWIFT linear solver: " << solver << Alert::Raise;
        p.model.h = Real(1) / Real(n - 1);
      }
  };
}
#endif
