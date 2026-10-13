/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Reference-to-physical gradient preparation through the H1 expression.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include "convergence/Convergence.h"

namespace Rodin::Tests::Benchmarks
{
  bool physicalGradientFailed = false;

  /**
   * @brief Measures the public gradient expression used by generic integrators.
   *
   * ## Architecture
   *
   * One retained affine cell and its owned mapped quadrature are used. The
   * inverse Jacobians and reference table are warmed before timing. Changing
   * the integration-point index requires physical-gradient preparation at
   * every sample, not a repeated-point cache-hit measurement. Basis extraction
   * and mesh/element queries inside the expression remain included; mapped
   * quadrature construction, local matrix contraction and assembly are not.
   * Specialized diffusion rules perform their own transformation loops, so
   * this expression scope is not a timing decomposition of those rules.
   *
   * An affine physical field @f$f(x)=\eta(1+\sum_j(j+1)x_j)@f$ is represented
   * through element DOF functionals. Reconstructing its gradient must recover
   * @f$\eta(1,\ldots,d)@f$; perturbing one contributing basis must fail.
   * Every physical basis component is also retained before timing and must
   * remain exactly unchanged after the sweep.
   *
   * @tparam K Scalar H1 degree.
   * @tparam Scalar Real or complex coefficient field.
   */
  template <size_t K, class Scalar>
  class PhysicalGradient
  {
    public:
    /// @brief Dimensionless affine-reproduction roundoff budget.
    static constexpr Real Tolerance = 1e-10;

    /**
     * @brief Sweeps prepared quadrature samples through the gradient expression.
     * @param state Requested polynomial quadrature order and repetitions.
     * @param geometry Unit-box cell family.
     */
    static void run(benchmark::State& state, Geometry::Polytope::Type geometry)
    {
      using namespace Variational;
      auto mesh = Convergence::UniformGrid(geometry).makeMesh(2);
      const auto cell = *mesh.getCell();
      const size_t d = cell.getDimension();
      H1<K, Scalar> space(std::integral_constant<size_t, K>{}, mesh);
      TrialFunction u(space);
      auto gradient = Grad(u);
      const auto& fe = space.getFiniteElement(d, cell.getIndex());
      const auto& formula = QF::PolytopeQuadratureFormula::get(state.range(0), geometry);
      const Geometry::PolytopeQuadrature quadrature(cell, formula);
      const Scalar multiplier = []() {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(1, 0.5);
        else
          return Scalar(1);
      }();
      Math::Vector<Scalar> coefficients(fe.getCount());
      for (size_t a = 0; a < fe.getCount(); ++a)
      {
        coefficients(a) = fe.getLinearForm(a)([&](const auto& reference) {
          const Geometry::Point point(cell, reference);
          Scalar value = multiplier;
          for (size_t j = 0; j < d; ++j)
            value += multiplier * Real(j + 1) * point.getPhysicalCoordinates()(j);
          return value;
        });
      }
      std::vector<Math::SpatialVector<Scalar>> initial;
      const auto valid = [&](const auto& values) {
        for (size_t qp = 0; qp < formula.getSize(); ++qp)
        {
          for (size_t j = 0; j < d; ++j)
          {
            Scalar value = 0;
            for (size_t a = 0; a < fe.getCount(); ++a)
              value += coefficients(a) * values[qp * fe.getCount() + a](j);
            if (!std::isfinite(std::abs(value)) ||
              std::abs(value - multiplier * Real(j + 1)) > Tolerance)
              return false;
          }
        }
        return true;
      };
      for (size_t qp = 0; qp < formula.getSize(); ++qp)
      {
        const IntegrationPoint ip(quadrature.getPoint(qp), &formula, qp);
        gradient.setIntegrationPoint(ip);
        for (size_t a = 0; a < fe.getCount(); ++a)
          initial.push_back(gradient.getBasis(a));
      }
      auto wrong = initial;
      size_t perturbed = 0;
      while (perturbed < fe.getCount() && std::abs(coefficients(perturbed)) == 0)
        ++perturbed;
      if (perturbed == fe.getCount())
      {
        physicalGradientFailed = true;
        state.SkipWithError("Affine coefficient representation is zero");
        return;
      }
      wrong[perturbed](0) += Scalar(1) / coefficients(perturbed);
      if (!valid(initial) || valid(wrong))
      {
        physicalGradientFailed = true;
        state.SkipWithError("Affine gradient or perturbation-rejection guard failed");
        return;
      }
      for (auto _ : state)
      {
        for (size_t qp = 0; qp < formula.getSize(); ++qp)
        {
          const IntegrationPoint ip(quadrature.getPoint(qp), &formula, qp);
          gradient.setIntegrationPoint(ip);
          for (size_t a = 0; a < fe.getCount(); ++a)
            benchmark::DoNotOptimize(gradient.getBasis(a));
        }
        benchmark::ClobberMemory();
      }
      for (size_t qp = 0; qp < formula.getSize(); ++qp)
      {
        const IntegrationPoint ip(quadrature.getPoint(qp), &formula, qp);
        gradient.setIntegrationPoint(ip);
        for (size_t a = 0; a < fe.getCount(); ++a)
        {
          const auto value = gradient.getBasis(a);
          for (size_t j = 0; j < d; ++j)
          {
            if (value(j) != initial[qp * fe.getCount() + a](j))
            {
              physicalGradientFailed = true;
              state.SkipWithError("Repeated gradient preparation changed a component");
              return;
            }
          }
        }
      }
      state.counters["dimension"] = d;
      state.counters["cells"] = 1;
      state.counters["degree"] = K;
      state.counters["basis_functions"] = fe.getCount();
      state.counters["quadrature_points"] = formula.getSize();
      state.counters["gradient_transform_products"] =
        formula.getSize() * fe.getCount() * d * d;
    }

    /// @brief Registers seven volume families and three quadrature orders.
    static void registerCases()
    {
      using G = Geometry::Polytope::Type;
      for (const auto geometry : {G::Segment, G::Triangle, G::Quadrilateral,
             G::Tetrahedron, G::Pyramid, G::Hexahedron, G::Wedge})
      {
        const std::string name = std::string("PhysicalGradient/H1/P") +
          std::to_string(K) + (std::is_same_v<Scalar, Complex> ? "/Complex/" : "/Real/") +
          std::string(Convergence::UniformGrid::getGeometryName(geometry));
        benchmark::RegisterBenchmark(
          name.c_str(), [geometry](auto& state) { run(state, geometry); })
          ->Arg(2)
          ->Arg(8)
          ->Arg(16)
          ->UseRealTime();
      }
    }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  PhysicalGradient<1, Real>::registerCases();
  PhysicalGradient<2, Real>::registerCases();
  PhysicalGradient<3, Real>::registerCases();
  PhysicalGradient<1, Complex>::registerCases();
  PhysicalGradient<2, Complex>::registerCases();
  PhysicalGradient<3, Complex>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "One-cell physical basis-gradient sweep; warmed inverses");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || physicalGradientFailed ? 1 : 0;
}
