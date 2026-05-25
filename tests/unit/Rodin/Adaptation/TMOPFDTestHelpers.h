/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_UNIT_ADAPTATION_TMOPFDTESTHELPERS_H
#define RODIN_TESTS_UNIT_ADAPTATION_TMOPFDTESTHELPERS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string_view>

#include <Rodin/Adaptation.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

namespace Rodin::Tests::Unit::TMOPFD
{
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;
  using namespace Rodin::Adaptation::TargetMatrixOptimization;

  static constexpr Attribute Negative = 1;
  static constexpr Attribute Positive = 2;
  static constexpr Attribute Boundary = 40;
  static constexpr Attribute Interface = 99;

  struct FDSweepResult
  {
    std::array<Real, 4> steps{{1e-5, 1e-6, 1e-7, 1e-8}};
    std::array<Real, 4> errors{{
      std::numeric_limits<Real>::infinity(),
      std::numeric_limits<Real>::infinity(),
      std::numeric_limits<Real>::infinity(),
      std::numeric_limits<Real>::infinity()
    }};

    Real bestError() const
    {
      return *std::min_element(errors.begin(), errors.end());
    }
  };

  inline bool hasEpsilonScalingTrend(const FDSweepResult& result)
  {
    const Real best = result.bestError();
    const Real e0 = result.errors[0];
    const Real e1 = result.errors[1];
    const Real e2 = result.errors[2];
    if (!std::isfinite(e0) || !std::isfinite(e1) || !std::isfinite(e2))
      return false;

    // Machine-floor regime: derivative consistency is already at numerical noise level.
    if (best <= Real(1e-10))
      return true;

    if (e0 <= Real(0) || e1 <= Real(0) || e2 <= Real(0))
      return false;

    const Real p01 = std::log(e0 / e1) / std::log(Real(10));
    const Real p12 = std::log(e1 / e2) / std::log(Real(10));
    return p01 >= Real(0.25) || p12 >= Real(0.25);
  }

  inline void printFDSweep(
      std::string_view name,
      std::string_view order,
      const FDSweepResult& result)
  {
    std::cout << "FD_TABLE " << name << " " << order << "\n";
    for (size_t i = 0; i < result.steps.size(); ++i)
      std::cout << "eps=" << result.steps[i] << " error=" << result.errors[i] << "\n";
  }

  inline LocalMesh makeTwoTriangleSquareWithVerticalInterface()
  {
    auto mesh = LocalMesh::Builder()
      .initialize(2)
      .nodes(6)
      .vertex({0, 0})
      .vertex({0.5, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({0.5, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 3})
      .polytope(Polytope::Type::Triangle, {1, 4, 3})
      .polytope(Polytope::Type::Triangle, {1, 2, 4})
      .polytope(Polytope::Type::Triangle, {2, 5, 4})
      .finalize();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(1, 2);
    mesh.setAttribute({2, 0}, Negative);
    mesh.setAttribute({2, 1}, Negative);
    mesh.setAttribute({2, 2}, Positive);
    mesh.setAttribute({2, 3}, Positive);
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto& edge = conn.getPolytope(1, e);
      const auto a = mesh.getVertexCoordinates(edge(0));
      const auto b = mesh.getVertexCoordinates(edge(1));
      if (std::abs(a[0] - Real(0.5)) <= Real(1e-12)
       && std::abs(b[0] - Real(0.5)) <= Real(1e-12))
      {
        mesh.setAttribute({1, e}, Interface);
      }
      else
      {
        const bool boundary =
          (std::abs(a[0]) <= Real(1e-12) && std::abs(b[0]) <= Real(1e-12))
       || (std::abs(a[0] - Real(1)) <= Real(1e-12)
        && std::abs(b[0] - Real(1)) <= Real(1e-12))
       || (std::abs(a[1]) <= Real(1e-12) && std::abs(b[1]) <= Real(1e-12))
       || (std::abs(a[1] - Real(1)) <= Real(1e-12)
        && std::abs(b[1] - Real(1)) <= Real(1e-12));
        if (boundary)
          mesh.setAttribute({1, e}, Boundary);
      }
    }
    return mesh;
  }

  inline void fillDisplacement(Eigen::VectorXd& data, Real amplitude)
  {
    for (Eigen::Index i = 0; i < data.size(); ++i)
    {
      const Real k = static_cast<Real>(i + 1);
      data(i) = amplitude * (std::sin(Real(0.37) * k)
               + Real(0.5) * std::cos(Real(0.19) * k));
    }
  }

  inline void fillDirection(Eigen::VectorXd& data)
  {
    for (Eigen::Index i = 0; i < data.size(); ++i)
    {
      const Real k = static_cast<Real>(i + 1);
      data(i) = std::cos(Real(0.23) * k)
              + Real(0.25) * std::sin(Real(0.71) * k);
    }
    data /= data.norm();
  }

  template <class AssembleResidual, class AssembleTangentAction>
  inline FDSweepResult finiteDifferenceConsistencySweep(
      Eigen::VectorXd& u,
      AssembleResidual assembleResidual,
      AssembleTangentAction assembleTangentAction)
  {
    FDSweepResult result;
    const auto u0 = u;
    auto direction = u;
    fillDirection(direction);
    for (size_t i = 0; i < result.steps.size(); ++i)
    {
      const Real eps = result.steps[i];
      const auto r0 = assembleResidual();
      u = u0 + eps * direction;
      const auto r1 = assembleResidual();
      u = u0;
      const auto fd = (r1 - r0) / eps;
      const auto jd = assembleTangentAction(direction);
      const Real denom = std::max<Real>({Real(1), fd.norm(), jd.norm()});
      result.errors[i] = (fd - jd).norm() / denom;
    }
    return result;
  }

  inline FDSweepResult tmopQualityFdSweep(bool useP2)
  {
    auto mesh = makeTwoTriangleSquareWithVerticalInterface();
    if (useP2)
    {
      VectorH1<2, LocalMesh> space(
          std::integral_constant<size_t, 2>{}, mesh, 2);
      GridFunction u(space);
      fillDisplacement(u.getData(), Real(0.002));
      TrialFunction du(space);
      TestFunction v(space);
      QualityTerm quality(
          ShapeSizeBlendMetric(Real(0.5)), IdentityTargetJacobian{}, Real(0.7));
      quality.setQuadratureOrder(4);

      auto assembleResidual = [&]()
      {
        LinearForm r(v);
        r = quality.residual(u, v);
        r.assemble();
        return r.getVector();
      };
      auto assembleTangentAction = [&](const Math::Vector<Real>& d)
      {
        BilinearForm j(du, v);
        j = quality.tangent(u, du, v);
        j.assemble();
        return (j.getOperator() * d).eval();
      };
        return finiteDifferenceConsistencySweep(
          u.getData(), assembleResidual, assembleTangentAction);
    }
    else
    {
      P1 space(mesh, 2);
      GridFunction u(space);
      fillDisplacement(u.getData(), Real(0.002));
      TrialFunction du(space);
      TestFunction v(space);
      QualityTerm quality(
          ShapeSizeBlendMetric(Real(0.5)), IdentityTargetJacobian{}, Real(0.7));
      quality.setQuadratureOrder(4);

      auto assembleResidual = [&]()
      {
        LinearForm r(v);
        r = quality.residual(u, v);
        r.assemble();
        return r.getVector();
      };
      auto assembleTangentAction = [&](const Math::Vector<Real>& d)
      {
        BilinearForm j(du, v);
        j = quality.tangent(u, du, v);
        j.assemble();
        return (j.getOperator() * d).eval();
      };
        return finiteDifferenceConsistencySweep(
          u.getData(), assembleResidual, assembleTangentAction);
    }
  }

      inline FDSweepResult deviationFdSweep(bool useP2)
  {
    auto mesh = makeTwoTriangleSquareWithVerticalInterface();
    if (useP2)
    {
      VectorH1<2, LocalMesh> space(
          std::integral_constant<size_t, 2>{}, mesh, 2);
      GridFunction u(space);
      fillDisplacement(u.getData(), Real(0.002));
      TrialFunction du(space);
      TestFunction v(space);
      DeviationTerm deviation(Real(0.25));

      auto assembleResidual = [&]()
      {
        LinearForm r(v);
        r = deviation.residual(u, v);
        r.assemble();
        return r.getVector();
      };
      auto assembleTangentAction = [&](const Math::Vector<Real>& d)
      {
        BilinearForm j(du, v);
        j = deviation.tangent(du, v);
        j.assemble();
        return (j.getOperator() * d).eval();
      };
        return finiteDifferenceConsistencySweep(
          u.getData(), assembleResidual, assembleTangentAction);
    }
    else
    {
      P1 space(mesh, 2);
      GridFunction u(space);
      fillDisplacement(u.getData(), Real(0.002));
      TrialFunction du(space);
      TestFunction v(space);
      DeviationTerm deviation(Real(0.25));

      auto assembleResidual = [&]()
      {
        LinearForm r(v);
        r = deviation.residual(u, v);
        r.assemble();
        return r.getVector();
      };
      auto assembleTangentAction = [&](const Math::Vector<Real>& d)
      {
        BilinearForm j(du, v);
        j = deviation.tangent(du, v);
        j.assemble();
        return (j.getOperator() * d).eval();
      };
        return finiteDifferenceConsistencySweep(
          u.getData(), assembleResidual, assembleTangentAction);
    }
  }

      inline FDSweepResult analyticFitFdSweep(bool useP2)
  {
    auto mesh = makeTwoTriangleSquareWithVerticalInterface();
    auto phi = [](const Math::SpatialPoint& x)
    {
      return x[0] * x[0] + Real(0.5) * x[1] * x[1] - Real(0.4);
    };
    auto gradPhi = [](const Math::SpatialPoint& x)
    {
      return Math::SpatialPoint{Real(2) * x[0], x[1]};
    };
    auto hessPhi = [](const Math::SpatialPoint&)
    {
      Math::SpatialMatrix<Real> h(2, 2);
      h(0, 0) = Real(2);
      h(0, 1) = Real(0);
      h(1, 0) = Real(0);
      h(1, 1) = Real(1);
      return h;
    };

    if (useP2)
    {
      VectorH1<2, LocalMesh> space(
          std::integral_constant<size_t, 2>{}, mesh, 2);
      GridFunction u(space);
      fillDisplacement(u.getData(), Real(0.002));
      TrialFunction du(space);
      TestFunction v(space);

      AnalyticLevelSetFitTerm fit(
        phi, gradPhi, hessPhi, Optional<Attribute>(Interface), Real(1.1));
      fit.setQuadratureOrder(4);
      fit.setNormalization(Real(1));

      auto assembleResidual = [&]()
      {
        LinearForm r(v);
        r = fit.residual(u, v);
        r.assemble();
        return r.getVector();
      };
      auto assembleTangentAction = [&](const Math::Vector<Real>& d)
      {
        BilinearForm j(du, v);
        j = fit.tangent(u, du, v);
        j.assemble();
        return (j.getOperator() * d).eval();
      };
        return finiteDifferenceConsistencySweep(
          u.getData(), assembleResidual, assembleTangentAction);
    }
    else
    {
      P1 space(mesh, 2);
      GridFunction u(space);
      fillDisplacement(u.getData(), Real(0.002));
      TrialFunction du(space);
      TestFunction v(space);

      AnalyticLevelSetFitTerm fit(
        phi, gradPhi, hessPhi, Optional<Attribute>(Interface), Real(1.1));
      fit.setQuadratureOrder(4);
      fit.setNormalization(Real(1));

      auto assembleResidual = [&]()
      {
        LinearForm r(v);
        r = fit.residual(u, v);
        r.assemble();
        return r.getVector();
      };
      auto assembleTangentAction = [&](const Math::Vector<Real>& d)
      {
        BilinearForm j(du, v);
        j = fit.tangent(u, du, v);
        j.assemble();
        return (j.getOperator() * d).eval();
      };
        return finiteDifferenceConsistencySweep(
          u.getData(), assembleResidual, assembleTangentAction);
    }
  }

      inline FDSweepResult phaseFdSweep(bool useP2)
  {
    auto mesh = makeTwoTriangleSquareWithVerticalInterface();
    auto phi = [](const Math::SpatialPoint& x)
    {
      return x[0] - Real(0.5);
    };
    auto gradPhi = [](const Math::SpatialPoint&)
    {
      return Math::SpatialPoint{1, 0};
    };

    if (useP2)
    {
      VectorH1<2, LocalMesh> space(
          std::integral_constant<size_t, 2>{}, mesh, 2);
      GridFunction u(space);
      fillDisplacement(u.getData(), Real(0.002));
      TrialFunction du(space);
      TestFunction v(space);

      VolumetricPhaseConsistencyTerm phase(
        phi, gradPhi, Negative, Positive, Real(0.9));
      phase
        .setQuadratureOrder(4)
        .setEpsilon(Real(0.5))
        .setMargin(Real(1))
        .setNormalization(Real(1));

      auto assembleResidual = [&]()
      {
        LinearForm r(v);
        r = phase.residual(u, v);
        r.assemble();
        return r.getVector();
      };
      auto assembleTangentAction = [&](const Math::Vector<Real>& d)
      {
        BilinearForm j(du, v);
        j = phase.tangent(u, du, v);
        j.assemble();
        return (j.getOperator() * d).eval();
      };
        return finiteDifferenceConsistencySweep(
          u.getData(), assembleResidual, assembleTangentAction);
    }
    else
    {
      P1 space(mesh, 2);
      GridFunction u(space);
      fillDisplacement(u.getData(), Real(0.002));
      TrialFunction du(space);
      TestFunction v(space);

      VolumetricPhaseConsistencyTerm phase(
        phi, gradPhi, Negative, Positive, Real(0.9));
      phase
        .setQuadratureOrder(4)
        .setEpsilon(Real(0.5))
        .setMargin(Real(1))
        .setNormalization(Real(1));

      auto assembleResidual = [&]()
      {
        LinearForm r(v);
        r = phase.residual(u, v);
        r.assemble();
        return r.getVector();
      };
      auto assembleTangentAction = [&](const Math::Vector<Real>& d)
      {
        BilinearForm j(du, v);
        j = phase.tangent(u, du, v);
        j.assemble();
        return (j.getOperator() * d).eval();
      };
      return finiteDifferenceConsistencySweep(
          u.getData(), assembleResidual, assembleTangentAction);
    }
  }

  inline Real tmopQualityFdError(bool useP2)
  {
    return tmopQualityFdSweep(useP2).bestError();
  }

  inline Real deviationFdError(bool useP2)
  {
    return deviationFdSweep(useP2).bestError();
  }

  inline Real analyticFitFdError(bool useP2)
  {
    return analyticFitFdSweep(useP2).bestError();
  }

  inline Real phaseFdError(bool useP2)
  {
    return phaseFdSweep(useP2).bestError();
  }
}

#endif
