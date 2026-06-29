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
    std::size_t quadratureOrder = 4;
    Real betaMax = 50;
    Real activeRMSOverHTol = 0;
    Real activeSupOverHTol = 0;
    Real gammaSize = 0;
    bool parseLegacyMaxIterations = false;
  };

  inline bool findOption(
      int argc,
      char** argv,
      const std::string& name,
      std::string* value)
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

  inline Real realOption(
      int argc,
      char** argv,
      const std::string& name,
      Real fallback)
  {
    std::string value;
    if (!findOption(argc, argv, name, &value) || value.empty())
      return fallback;
    return static_cast<Real>(std::atof(value.c_str()));
  }

  inline std::size_t sizeOption(
      int argc,
      char** argv,
      const std::string& name,
      std::size_t fallback)
  {
    std::string value;
    if (!findOption(argc, argv, name, &value) || value.empty())
      return fallback;
    return static_cast<std::size_t>(std::strtoull(value.c_str(), nullptr, 10));
  }

  inline bool boolOption(
      int argc,
      char** argv,
      const std::string& name,
      bool fallback)
  {
    std::string value;
    if (!findOption(argc, argv, name, &value))
      return fallback;
    if (value.empty())
      return true;
    return std::strtoull(value.c_str(), nullptr, 10) != 0;
  }

  inline Adaptation::WNGIRParameters makeWNGIRParameters(
      int argc,
      char** argv,
      Real h,
      Geometry::Attribute interfaceAttribute,
      const WNGIRExampleDefaults& defaults = {})
  {
    Adaptation::WNGIRParameters p;
    p.h = h;

    const Real gammaMFactor =
      realOption(argc, argv, "wngir-gamma-m", Real(1));
    const Real gammaHFactor =
      realOption(argc, argv, "wngir-gamma-h", Real(1));
    const Real gammaDivFactor =
      realOption(argc, argv, "wngir-gamma-div", gammaHFactor);
    p.gammaM = gammaMFactor / h;
    p.gammaH = gammaHFactor / h;
    p.gammaDiv = gammaDivFactor / h;
    p.ellM = realOption(argc, argv, "wngir-ell", Real(3)) * h;

    p.gammaObs = realOption(argc, argv, "wngir-gamma-obs", Real(1));
    p.residualStabilizedObservationMetric =
      boolOption(argc, argv, "wngir-residual-stabilized-obs", true);
    p.betaMax = realOption(argc, argv, "wngir-beta-max", defaults.betaMax);

    p.gammaJ = realOption(argc, argv, "wngir-gamma-j", Real(1));
    p.gammaQ = realOption(argc, argv, "wngir-gamma-q", Real(1));
    p.jSafe = realOption(argc, argv, "wngir-jsafe", Real(1e-2));
    p.qMax = realOption(argc, argv, "wngir-qmax", Real(10));
    p.s0J = realOption(argc, argv, "wngir-s0j", Real(0.25));
    p.s0Q = realOption(argc, argv, "wngir-s0q", Real(2));

    p.gammaQual = realOption(argc, argv, "wngir-gamma-qual", Real(1));
    p.qStar = realOption(argc, argv, "wngir-qstar", Real(1.75));
    p.gammaSize = realOption(argc, argv, "wngir-gamma-size", defaults.gammaSize);
    p.jStar = realOption(argc, argv, "wngir-jstar", Real(0.3));

    p.omegaMin = realOption(argc, argv, "wngir-omega-min", Real(0.1));
    p.alphaMin = realOption(argc, argv, "wngir-alpha-min", Real(1e-4));
    p.energyLineSearch = !boolOption(argc, argv, "wngir-no-energy-ls", false);

    p.jMinRatio = realOption(argc, argv, "j-min", Real(1e-8));
    const Real jSafeRatio = realOption(argc, argv, "j-safe", Real(1e-3));
    p.jLineSearchRatio =
      realOption(argc, argv, "j-ls", std::max(p.jMinRatio, Real(10) * jSafeRatio));

    p.activeRMSTol = realOption(argc, argv, "wngir-rms-tol", Real(4) * h * h);
    p.activeSupTol = realOption(argc, argv, "wngir-sup-tol", Real(10) * h * h);
    p.activeRMSOverHTol =
      realOption(argc, argv, "wngir-rms-h-tol", defaults.activeRMSOverHTol);
    p.activeSupOverHTol =
      realOption(argc, argv, "wngir-sup-h-tol", defaults.activeSupOverHTol);
    p.energyStagTol =
      realOption(argc, argv, "wngir-energy-stag-tol", Real(1e-4));
    p.stepTol = realOption(argc, argv, "wngir-step-tol", Real(1e-4) * h);

    p.quadratureOrder =
      sizeOption(argc, argv, "quad-order", defaults.quadratureOrder);
    p.maxIterations =
      sizeOption(argc, argv, "wngir-steps", defaults.maxIterations);
    if (defaults.parseLegacyMaxIterations)
      p.maxIterations =
        sizeOption(argc, argv, "wngir-max-iters", p.maxIterations);

    p.cgRelativeTolerance =
      realOption(argc, argv, "wngir-cg-rtol", p.cgRelativeTolerance);
    p.cgMaxIterations =
      sizeOption(argc, argv, "wngir-cg-max-iters", p.cgMaxIterations);
    p.andersonMemory =
      sizeOption(argc, argv, "wngir-aa-memory", p.andersonMemory);
    p.andersonStart =
      sizeOption(argc, argv, "wngir-aa-start", p.andersonStart);
    p.andersonDamping =
      realOption(argc, argv, "wngir-aa-damping", p.andersonDamping);
    p.andersonMinDamping =
      realOption(argc, argv, "wngir-aa-min-damping", p.andersonMinDamping);

    p.includeQualityGradient =
      boolOption(argc, argv, "wngir-quality-energy", false);
    p.includeQualityMetric =
      boolOption(argc, argv, "wngir-quality-metric", true);
    p.includeAdmissibilityMetric =
      boolOption(argc, argv, "wngir-admissibility-metric", true);
    p.includeAdmissibilityGradient =
      boolOption(argc, argv, "wngir-admissibility-gradient", false);
    p.splitQualityDirection =
      boolOption(argc, argv, "wngir-split-quality", false);
    p.qualityDirectionWeight =
      realOption(argc, argv, "wngir-quality-weight", Real(1));

    p.hasInterfaceAttribute = true;
    p.interfaceAttribute = interfaceAttribute;
    p.trace =
      boolOption(argc, argv, "trace",
        boolOption(argc, argv, "wngir-trace", false));
    return p;
  }
}

#endif
