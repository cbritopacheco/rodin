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
#include <string>

#include <Rodin/Adaptation/WNGIRParameters.h>

namespace Rodin::Examples
{
  struct WNGIRExampleDefaults
  {
      std::size_t maxIterations = 60;
      std::size_t quadratureOrder = 0;
      Real kappaBulk = Real(0.00703125);
      Real rDiv = 1;
      Real kappaJ = 1;
      Real kappaQ = 1;
      Real tauRmsHFloor = Real(0.03);
      Real tauInfHFloor = Real(0.20);
      bool parseLegacyMaxIterations = false;
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

  inline Adaptation::WNGIRParameters makeWNGIRParameters(int argc, char** argv, Real h,
    Geometry::Attribute interfaceAttribute, const WNGIRExampleDefaults& defaults = {})
  {
    Adaptation::WNGIRParameters p;
    p.h = h;

    p.kappaBulk =
      realOption(argc, argv, "wngir-kappa-bulk", defaults.kappaBulk);
    p.rDiv = realOption(
      argc, argv, "wngir-r-div", "wngir-divergence-ratio", defaults.rDiv);

    p.kappaObs = realOption(argc, argv, "wngir-kappa-obs", Real(1));
    p.robustScale = realOption(argc, argv, "wngir-robust-scale", p.robustScale);

    p.kappaJ = realOption(argc, argv, "wngir-kappa-j", defaults.kappaJ);
    p.kappaQ = realOption(argc, argv, "wngir-kappa-q", defaults.kappaQ);
    p.jSafe = realOption(argc, argv, "wngir-jsafe", Real(1e-2));
    p.qMax = realOption(argc, argv, "wngir-qmax", Real(10));
    p.primalBarrierIterations = std::max<std::size_t>(1,
      sizeOption(
        argc, argv, "wngir-primal-barrier-iterations", p.primalBarrierIterations));
    p.primalBarrierRelativeTolerance = realOption(
      argc, argv, "wngir-primal-barrier-relative-tol", p.primalBarrierRelativeTolerance);
    p.muHat =
      realOption(argc, argv, "wngir-mu-hat", "wngir-primal-barrier-mu", p.muHat);
    p.thetaBoundary = realOption(argc, argv, "wngir-theta-boundary",
      "wngir-fraction-to-boundary", p.thetaBoundary);
    p.omegaMin = realOption(argc, argv, "wngir-omega-min", Real(0.1));
    p.alphaMin = realOption(argc, argv, "wngir-alpha-min", Real(1e-4));
    p.armijoCoefficient = realOption(argc, argv, "wngir-armijo", p.armijoCoefficient);
    p.descentFraction =
      realOption(argc, argv, "wngir-descent-fraction", p.descentFraction);
    p.directionNormFactor =
      realOption(argc, argv, "wngir-direction-norm-factor", p.directionNormFactor);

    p.jMinRatio = realOption(argc, argv, "j-min", Real(1e-8));
    const Real jSafeRatio = realOption(argc, argv, "j-safe", Real(1e-3));
    p.jLineSearchRatio =
      realOption(argc, argv, "j-ls", std::max(p.jMinRatio, Real(10) * jSafeRatio));
    p.tauRmsHFloor = realOption(argc, argv, "wngir-rms-floor", defaults.tauRmsHFloor);
    p.tauInfHFloor = realOption(argc, argv, "wngir-sup-floor", defaults.tauInfHFloor);
    p.tauJumpRms =
      realOption(argc, argv, "wngir-rms-normal-jump-factor", p.tauJumpRms);
    p.tauJumpInf =
      realOption(argc, argv, "wngir-sup-normal-jump-factor", p.tauJumpInf);
    // Zero delegates the physical tolerance to WNGIR, where the sampled
    // level-set gradient converts mesh length to field units.
    p.tauRms = realOption(argc, argv, "wngir-rms-tol", Real(0));
    p.tauInf = realOption(argc, argv, "wngir-sup-tol", Real(0));
    p.energyStagTol = realOption(argc, argv, "wngir-energy-stag-tol", Real(1e-4));
    p.stepTol = realOption(argc, argv, "wngir-step-tol", Real(1e-4) * h);
    p.acceptedStepOverHTol =
      realOption(argc, argv, "wngir-step-h-tol", p.acceptedStepOverHTol);

    p.quadratureOrder = sizeOption(argc, argv, "quad-order", defaults.quadratureOrder);
    p.maxIterations = sizeOption(argc, argv, "wngir-steps", defaults.maxIterations);
    if (defaults.parseLegacyMaxIterations)
      p.maxIterations = sizeOption(argc, argv, "wngir-max-iters", p.maxIterations);

    p.rigidStabilisationLevel = realOption(
      argc, argv, "wngir-rigid-stabilisation", p.rigidStabilisationLevel);
    p.cgRelativeTolerance =
      realOption(argc, argv, "wngir-cg-rtol", p.cgRelativeTolerance);
    p.cgMaxIterations = sizeOption(argc, argv, "wngir-cg-max-iters", p.cgMaxIterations);
    p.hasInterfaceAttribute = true;
    p.interfaceAttribute = interfaceAttribute;
    p.trace =
      boolOption(argc, argv, "trace", boolOption(argc, argv, "wngir-trace", false));
    return p;
  }
}

#endif
