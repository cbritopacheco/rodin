/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Reference-to-physical derivative preparation through H1 expressions.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include "convergence/Convergence.h"

namespace Rodin::Tests::Benchmarks
{
  bool physicalGradientFailed = false;

  enum class PhysicalDerivative
  {
    Gradient,
    Jacobian,
    Divergence
  };

  /**
   * @brief Measures public physical derivative expressions used by generic integrators.
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
   * through element DOF functionals. The scalar gradient must recover
   * @f$\eta(1,\ldots,d)@f$. Vector scopes use @f$u_c=(c+1)f@f$, with
   * @f$(\nabla u)_{cj}=\eta(c+1)(j+1)@f$ and
   * @f$\operatorname{div}u=\eta\sum_j(j+1)^2@f$.
   * Perturbing one contributing basis must fail.
   * Every physical basis component is also retained before timing and must
   * remain exactly unchanged after the sweep.
   * Native P1 scopes additionally compare every indexed basis derivative with
   * fresh direct-point evaluation. Identity-map scopes distinguish geometric
   * map order from reference-basis sample dependence. On affine simplices the
   * warmed native P1 cache is intentionally reused throughout timing; the
   * transform-product counter is a logical work proxy, not an instruction count.
   * Native P1 basis access is materialized into an owned value before the
   * mutable benchmark barrier; no deprecated const-reference barrier is used.
   * The retained H1 value-return protocol is unchanged.
   *
   * @tparam K Scalar H1 degree.
   * @tparam Scalar Real or complex coefficient field.
   * @tparam Operation Physical differential expression to prepare.
   * @tparam NativeP1 Whether to use native P1 rather than H1.
   * @tparam IdentityMap Whether the retained cell uses an identity transformation.
   */
  template <size_t K, class Scalar,
    PhysicalDerivative Operation = PhysicalDerivative::Gradient, bool NativeP1 = false,
    bool IdentityMap = false>
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
        if constexpr (IdentityMap)
        {
          mesh.setPolytopeTransformation(
            {d, cell.getIndex()}, new Geometry::IdentityTransformation(d));
        }
        const auto space = [&]() {
          if constexpr (NativeP1 && Operation == PhysicalDerivative::Gradient)
            return P1<Scalar>(mesh);
          else if constexpr (NativeP1)
            return P1<Math::SpatialVector<Scalar>>(mesh, d);
          else if constexpr (Operation == PhysicalDerivative::Gradient)
            return H1<K, Scalar>(std::integral_constant<size_t, K>{}, mesh);
          else
          {
            return H1<K, Math::SpatialVector<Scalar>>(
              std::integral_constant<size_t, K>{}, mesh, d);
          }
        }();
        TrialFunction u(space);
        auto gradient = [&]() {
          if constexpr (Operation == PhysicalDerivative::Gradient)
            return Grad(u);
          else if constexpr (Operation == PhysicalDerivative::Jacobian)
            return Jacobian(u);
          else
            return Div(u);
        }();
        gradient.setPolytope(cell);
        using Value = std::remove_cvref_t<decltype(gradient.getBasis(size_t(0)))>;
        const auto& fe = space.getFiniteElement(d, cell.getIndex());
        const auto& formula = QF::QuadratureFormula::get(state.range(0), geometry);
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
            if constexpr (Operation == PhysicalDerivative::Gradient)
              return value;
            else
            {
              Math::SpatialVector<Scalar> vector(d);
              for (size_t c = 0; c < d; ++c)
                vector(c) = Scalar(c + 1) * value;
              return vector;
            }
          });
        }
        std::vector<Value> initial;
        const auto valid = [&](const auto& values) {
          for (size_t qp = 0; qp < formula.getSize(); ++qp)
          {
            if constexpr (Operation == PhysicalDerivative::Divergence)
            {
              Scalar value = 0;
              for (size_t a = 0; a < fe.getCount(); ++a)
                value += coefficients(a) * values[qp * fe.getCount() + a];
              Scalar expected = 0;
              for (size_t j = 0; j < d; ++j)
                expected += multiplier * Real((j + 1) * (j + 1));
              if (!std::isfinite(std::abs(value)) ||
                std::abs(value - expected) > Tolerance)
                return false;
            }
            else
            {
              for (size_t c = 0; c < (Operation == PhysicalDerivative::Gradient ? 1 : d);
                   ++c)
              {
                for (size_t j = 0; j < d; ++j)
                {
                  Scalar value = 0;
                  for (size_t a = 0; a < fe.getCount(); ++a)
                  {
                    const auto& basis = values[qp * fe.getCount() + a];
                    if constexpr (Operation == PhysicalDerivative::Gradient)
                      value += coefficients(a) * basis(j);
                    else
                      value += coefficients(a) * basis(c, j);
                  }
                  if (!std::isfinite(std::abs(value)) ||
                    std::abs(value - multiplier * Real((c + 1) * (j + 1))) > Tolerance)
                    return false;
                }
              }
            }
          }
          return true;
        };
        for (size_t qp = 0; qp < formula.getSize(); ++qp)
        {
          const IntegrationPoint ip(quadrature, qp);
          gradient.setIntegrationPoint(ip);
          for (size_t a = 0; a < fe.getCount(); ++a)
            initial.push_back(gradient.getBasis(a));
          if constexpr (NativeP1)
          {
            auto direct = [&]() {
              if constexpr (Operation == PhysicalDerivative::Gradient)
                return Grad(u);
              else if constexpr (Operation == PhysicalDerivative::Jacobian)
                return Jacobian(u);
              else
                return Div(u);
            }();
            const IntegrationPoint pointwise(quadrature.getPoint(qp));
            direct.setIntegrationPoint(pointwise);
            for (size_t a = 0; a < fe.getCount(); ++a)
            {
              const auto expected = direct.getBasis(a);
              const auto& actual = initial[qp * fe.getCount() + a];
              const bool equal = [&]() {
                if constexpr (Operation == PhysicalDerivative::Divergence)
                  return actual == expected;
                else
                {
                  for (size_t c = 0;
                       c < (Operation == PhysicalDerivative::Gradient ? 1 : d); ++c)
                  {
                    for (size_t j = 0; j < d; ++j)
                    {
                      if constexpr (Operation == PhysicalDerivative::Gradient)
                      {
                        if (actual(j) != expected(j))
                          return false;
                      }
                      else if (actual(c, j) != expected(c, j))
                        return false;
                    }
                  }
                  return true;
                }
              }();
              if (!equal)
              {
                physicalGradientFailed = true;
                state.SkipWithError(
                  "Indexed basis differs from fresh direct-point derivative");
                return;
              }
            }
          }
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
        if constexpr (Operation == PhysicalDerivative::Gradient)
          wrong[perturbed](0) += Scalar(1) / coefficients(perturbed);
        else if constexpr (Operation == PhysicalDerivative::Jacobian)
          wrong[perturbed](0, 0) += Scalar(1) / coefficients(perturbed);
        else
          wrong[perturbed] += Scalar(1) / coefficients(perturbed);
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
            const IntegrationPoint ip(quadrature, qp);
            gradient.setIntegrationPoint(ip);
            for (size_t a = 0; a < fe.getCount(); ++a)
            {
              if constexpr (NativeP1)
              {
                auto value = gradient.getBasis(a);
                benchmark::DoNotOptimize(value);
              }
              else
                benchmark::DoNotOptimize(gradient.getBasis(a));
            }
          }
          benchmark::ClobberMemory();
        }
        for (size_t qp = 0; qp < formula.getSize(); ++qp)
        {
          const IntegrationPoint ip(quadrature, qp);
          gradient.setIntegrationPoint(ip);
          for (size_t a = 0; a < fe.getCount(); ++a)
          {
            const auto value = gradient.getBasis(a);
            const auto& expected = initial[qp * fe.getCount() + a];
            const bool unchanged = [&]() {
              if constexpr (Operation == PhysicalDerivative::Divergence)
                return value == expected;
              else
              {
                for (size_t c = 0;
                     c < (Operation == PhysicalDerivative::Gradient ? 1 : d); ++c)
                {
                  for (size_t j = 0; j < d; ++j)
                  {
                    if constexpr (Operation == PhysicalDerivative::Gradient)
                    {
                      if (value(j) != expected(j))
                        return false;
                    }
                    else if (value(c, j) != expected(c, j))
                      return false;
                  }
                }
                return true;
              }
            }();
            if (!unchanged)
            {
              physicalGradientFailed = true;
              state.SkipWithError("Repeated derivative preparation changed a component");
              return;
            }
          }
        }
        state.counters["dimension"] = d;
        state.counters["geometry_order"] = cell.getTransformation().getOrder();
        state.counters["identity_map"] = IdentityMap;
        const bool affineReference = geometry == Geometry::Polytope::Type::Segment ||
          geometry == Geometry::Polytope::Type::Triangle ||
          geometry == Geometry::Polytope::Type::Tetrahedron;
        state.counters["warm_affine_reuse"] =
          NativeP1 && affineReference && cell.getTransformation().getOrder() == 1;
        state.counters["cells"] = 1;
        state.counters["degree"] = K;
        state.counters["basis_functions"] = fe.getCount();
        state.counters["quadrature_points"] = formula.getSize();
        state.counters["field_components"] =
          Operation == PhysicalDerivative::Gradient ? 1 : d;
        state.counters["local_dofs"] = fe.getCount();
        const size_t scalarBasis =
          Operation == PhysicalDerivative::Gradient ? fe.getCount() : fe.getCount() / d;
        state.counters["scalar_basis_functions"] = scalarBasis;
        state.counters["gradient_transform_products"] =
          formula.getSize() * scalarBasis * d * d;
      }

      /// @brief Registers seven volume families, including elasticity-matched orders.
      static void registerCases()
      {
        using G = Geometry::Polytope::Type;
        for (const auto geometry : {G::Segment, G::Triangle, G::Quadrilateral,
               G::Tetrahedron, G::Pyramid, G::Hexahedron, G::Wedge})
        {
          const std::string operation = Operation == PhysicalDerivative::Gradient
            ? "PhysicalGradient"
            : Operation == PhysicalDerivative::Jacobian ? "PhysicalJacobian"
                                                        : "PhysicalDivergence";
          const std::string family = NativeP1
            ? (IdentityMap ? "/NativeP1Identity" : "/NativeP1")
            : "/H1/P" + std::to_string(K);
          const std::string name = operation + family +
            (std::is_same_v<Scalar, Complex> ? "/Complex/" : "/Real/") +
            std::string(Convergence::UniformGrid::getGeometryName(geometry));
          benchmark::RegisterBenchmark(
            name.c_str(), [geometry](auto& state) { run(state, geometry); })
            ->Arg(2)
            ->Arg(4)
            ->Arg(6)
            ->Arg(8)
            ->Arg(12)
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
  PhysicalGradient<1, Real, PhysicalDerivative::Jacobian>::registerCases();
  PhysicalGradient<2, Real, PhysicalDerivative::Jacobian>::registerCases();
  PhysicalGradient<3, Real, PhysicalDerivative::Jacobian>::registerCases();
  PhysicalGradient<1, Complex, PhysicalDerivative::Jacobian>::registerCases();
  PhysicalGradient<2, Complex, PhysicalDerivative::Jacobian>::registerCases();
  PhysicalGradient<3, Complex, PhysicalDerivative::Jacobian>::registerCases();
  PhysicalGradient<1, Real, PhysicalDerivative::Divergence>::registerCases();
  PhysicalGradient<2, Real, PhysicalDerivative::Divergence>::registerCases();
  PhysicalGradient<3, Real, PhysicalDerivative::Divergence>::registerCases();
  PhysicalGradient<1, Complex, PhysicalDerivative::Divergence>::registerCases();
  PhysicalGradient<2, Complex, PhysicalDerivative::Divergence>::registerCases();
  PhysicalGradient<3, Complex, PhysicalDerivative::Divergence>::registerCases();
  const auto registerNative = []<class Scalar, bool Identity>() {
    PhysicalGradient<1, Scalar, PhysicalDerivative::Gradient, true,
      Identity>::registerCases();
    PhysicalGradient<1, Scalar, PhysicalDerivative::Jacobian, true,
      Identity>::registerCases();
    PhysicalGradient<1, Scalar, PhysicalDerivative::Divergence, true,
      Identity>::registerCases();
  };
  registerNative.template operator()<Real, false>();
  registerNative.template operator()<Complex, false>();
  registerNative.template operator()<Real, true>();
  registerNative.template operator()<Complex, true>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "One-cell physical basis-derivative sweep; warmed inverses");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || physicalGradientFailed ? 1 : 0;
}
