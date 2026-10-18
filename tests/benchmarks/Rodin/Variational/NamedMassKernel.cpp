/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Named mass-kernel lifecycle and within-assembly reference reuse.
 */
#include <benchmark/benchmark.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <type_traits>

#include "Rodin/Variational.h"
#include "../../../QuadratureReference.h"
#include "../../../convergence/Convergence.h"

namespace Rodin::Tests::Benchmarks
{
  /// @brief Nonempty H1 value shapes measured independently of spatial dimension.
  enum class NamedMassShape
  {
    Scalar,
    Vector,
    Matrix
  };

  /// @brief Propagates numerical qualification failures to the executable result.
  bool namedMassKernelFailed = false;

  /**
   * @brief Measures the production named mass kernel without global assembly.
   *
   * ## Mathematical foundation
   *
   * For the component basis @f$\Phi_{a,c}=\phi_a E_c@f$, the local matrix is
   * @f$A_{bm+e,am+c}=\delta_{ec}S_{ba}@f$. The scalar block is evaluated
   * independently by the original expression-entry quadrature loop. Weighted
   * scopes use @f$\rho=s\alpha@f$, with @f$\alpha=2@f$ over the reals and
   * @f$\alpha=2+i/2@f$ over the complex numbers. The borrowed scale is one
   * during timing and two during the final coefficient-replacement check.
   * Unweighted scopes use @f$\rho=1@f$. A complex coefficient does not imply
   * Hermitian symmetry. The requested order is the production basis-order sum;
   * comparison certifies that discrete operator, not rational-pyramid exactness.
   *
   * ## Architecture and measurement boundary
   *
   * FreshKernel includes construction, coefficient cloning where applicable,
   * reference-cache creation, computation and kernel destruction. Output storage
   * is retained. CachedKernel retains the same kernel and its prepared reference
   * data across computations. This models reuse across cells within one named
   * assembly, not retention of that cache across MassForm::assemble calls.
   * Geometry and canonical quadrature caches are primed in both scopes.
   * Mesh/space construction, the scalar oracle and numerical checks are untimed;
   * scatter, constraints, solving and process-first geometry are not measured.
   *
   * Reference-cache payload is a source-derived model, not measured resident
   * memory. Degrees one through three bound this initial profiling matrix;
   * the high-degree component-expanded cache is investigated separately before
   * launching its multi-gigabyte configurations. Point cases are structural
   * counting-measure operations, not refinement studies. Empty shapes, mixed
   * degrees, other space families and curved/embedded cells are separate scopes.
   *
   * @tparam K Nominal H1 degree.
   * @tparam Scalar Real or complex scalar field.
   * @tparam Shape Scalar, vector or matrix value shape.
   * @tparam Weighted Whether the kernel owns a scalar coefficient.
   * @tparam Warm Whether the same prepared kernel is retained during timing.
   */
  template <size_t K, class Scalar, NamedMassShape Shape, bool Weighted, bool Warm>
  class NamedMassKernel
  {
    public:
      /// @brief Existing named-form normalized Frobenius comparison budget.
      static constexpr Real OperatorBudget = 1e-11;

      /**
       * @brief Measures one production local-kernel lifecycle or reuse scope.
       * @param state Field rows and columns, repetitions and timing state.
       * @param geometry Reference-cell family, including Point.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry)
      {
        using namespace Geometry;
        using namespace Variational;
        const size_t rows = state.range(0), columns = state.range(1);
        const size_t components = rows * columns;
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        auto mesh = [&]() {
          if (dimension == 0)
            return LocalMesh::Builder()
              .initialize(3)
              .nodes(1)
              .vertex({0, 0, 0})
              .finalize();
          return LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
        }();
        for (size_t d = 1; d <= dimension; ++d)
        {
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
        const auto cell = *mesh.getCell();
        auto space = [&]() {
          if constexpr (Shape == NamedMassShape::Scalar)
            return H1<K, Scalar>(std::integral_constant<size_t, K>{}, mesh);
          else if constexpr (Shape == NamedMassShape::Vector)
            return H1<K, Math::SpatialVector<Scalar>>(
              std::integral_constant<size_t, K>{}, mesh, rows);
          else
            return H1<K, Math::SpatialMatrix<Scalar>>(
              std::integral_constant<size_t, K>{}, mesh, rows, columns);
        }();
        H1<K, Scalar> scalarSpace(std::integral_constant<size_t, K>{}, mesh);
        TrialFunction u(space);
        TestFunction v(space);
        TrialFunction scalarU(scalarSpace);
        TestFunction scalarV(scalarSpace);
        Real scale = 1;
        const Scalar alpha = [] {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Scalar(2, 0.5);
          else
            return Scalar(2);
        }();
        [[maybe_unused]] const auto coefficient = [&]() {
          const auto evaluate = [alpha, &scale]([[maybe_unused]] const Point& point) {
            return alpha * scale;
          };
          if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction(evaluate);
          else
            return RealFunction(evaluate);
        }();
        const auto makeKernel = [&]() {
          if constexpr (Weighted)
          {
            using Form = decltype(MassForm(coefficient, u, v));
            return typename Form::KernelType(space, space, coefficient);
          }
          else
          {
            using Form = decltype(MassForm(u, v));
            return typename Form::KernelType(space, space);
          }
        };
        const auto scalarIntegral = [&]() {
          if constexpr (Weighted)
            return Integral(coefficient * scalarU, scalarV);
          else
            return Integral(scalarU, scalarV);
        }();
        const auto& element = space.getFiniteElement(dimension, cell.getIndex());
        const size_t count = element.getCount();
        const size_t scalarCount =
          scalarSpace.getFiniteElement(dimension, cell.getIndex()).getCount();
        if (count != components * scalarCount || count == 0)
        {
          namedMassKernelFailed = true;
          state.SkipWithError("Named mass fixture has an invalid component layout");
          return;
        }
        const size_t order = 2 * element.getOrder();
        QuadratureReference reference(scalarIntegral);
        reference.assemble(cell, order);
        Math::Matrix<Scalar> expected = Math::Matrix<Scalar>::Zero(count, count);
        for (size_t b = 0; b < scalarCount; ++b)
        {
          for (size_t a = 0; a < scalarCount; ++a)
          {
            for (size_t c = 0; c < components; ++c)
            {
              expected(b * components + c, a * components + c) =
                reference.getOperator()(b, a);
            }
          }
        }
        const auto error = [&](const Math::Matrix<Scalar>& value,
                             Real multiplier) -> Optional<Real> {
          if (value.rows() != expected.rows() || value.cols() != expected.cols())
            return {};
          return (value - multiplier * expected).norm() /
            std::max(Real(1), multiplier * expected.norm());
        };
        Math::Matrix<Scalar> output;
        Real maxError = 0;
        {
          auto kernel = makeKernel();
          kernel.compute(output, cell);
          const auto initialError = error(output, 1);
          if (!initialError || !std::isfinite(*initialError) ||
            *initialError > OperatorBudget || expected.norm() == 0)
          {
            namedMassKernelFailed = true;
            state.SkipWithError("Named mass kernel failed its independent block oracle");
            return;
          }
          maxError = *initialError;
          Math::Matrix<Scalar> perturbed = output;
          perturbed(0, 0) += Scalar(1);
          const auto perturbedError = error(perturbed, 1);
          if (!perturbedError || !std::isfinite(*perturbedError) ||
            *perturbedError <= OperatorBudget)
          {
            namedMassKernelFailed = true;
            state.SkipWithError("Named mass oracle accepted a changed matrix entry");
            return;
          }
        }
        const auto check = [&](Real multiplier) {
          const auto observed = error(output, multiplier);
          if (!observed || !std::isfinite(*observed) || *observed > OperatorBudget)
          {
            namedMassKernelFailed = true;
            state.SkipWithError(
              "Named mass kernel failed repetition or coefficient replacement");
            return false;
          }
          maxError = std::max(maxError, *observed);
          return true;
        };
        if constexpr (Warm)
        {
          auto kernel = makeKernel();
          kernel.compute(output, cell);
          for (auto _ : state)
          {
            kernel.compute(output, cell);
            benchmark::DoNotOptimize(output.data());
            benchmark::ClobberMemory();
          }
          if (!check(1))
            return;
          if constexpr (Weighted)
          {
            scale = 2;
            kernel.compute(output, cell);
            if (!check(2))
              return;
          }
        }
        else
        {
          for (auto _ : state)
          {
            auto kernel = makeKernel();
            kernel.compute(output, cell);
            benchmark::DoNotOptimize(output.data());
            benchmark::ClobberMemory();
          }
          if (!check(1))
            return;
          if constexpr (Weighted)
          {
            scale = 2;
            auto kernel = makeKernel();
            kernel.compute(output, cell);
            if (!check(2))
              return;
          }
        }
        const auto& qf = QF::PolytopeQuadratureFormula::get(order, geometry);
        state.counters["degree"] = K;
        state.counters["dimension"] = dimension;
        state.counters["ambient_dimension"] = mesh.getSpaceDimension();
        state.counters["components"] = components;
        state.counters["field_rows"] = rows;
        state.counters["field_columns"] = columns;
        state.counters["cells"] = 1;
        state.counters["fixture_cells"] = mesh.getPolytopeCount(dimension);
        state.counters["global_dofs"] = space.getSize();
        state.counters["geometry_factor_order"] =
          cell.getTransformation().getFactorOrder();
        state.counters["scalar_basis"] = scalarCount;
        state.counters["local_basis"] = count;
        state.counters["local_entries"] = count * count;
        state.counters["quadrature_order"] = order;
        state.counters["quadrature_points"] = qf.getSize();
        state.counters["scalar_pair_sample_proxy"] =
          qf.getSize() * scalarCount * scalarCount;
        state.counters["component_pair_sample_proxy"] =
          qf.getSize() * scalarCount * scalarCount * components;
        state.counters["dense_entry_sample_proxy"] = qf.getSize() * count * count;
        state.counters["expanded_reference_payload_bytes_model"] =
          qf.getSize() * count * count * sizeof(Scalar);
        state.counters["scalar_reference_payload_bytes_model"] =
          qf.getSize() * scalarCount * scalarCount * sizeof(Scalar);
        state.counters["output_bytes"] = count * count * sizeof(Scalar);
        state.counters["scalar_bytes"] = sizeof(Scalar);
        state.counters["reference_cold"] = !Warm;
        state.counters["output_storage_retained"] = 1;
        state.counters["geometry_cache_primed"] = 1;
        state.counters["weighted"] = Weighted;
        state.counters["coefficient_replacement_checked"] = Weighted;
        state.counters["reference_error"] = maxError;
      }
  };

  /**
   * @brief Registers independent named-kernel reuse and nonempty shape scopes.
   * @tparam K H1 degree, bounded by three in this initial matrix.
   * @tparam Scalar Real or complex field scalar.
   * @tparam Shape Scalar, vector or matrix range.
   * @tparam Weighted Whether a scalar coefficient is owned.
   */
  template <size_t K, class Scalar, NamedMassShape Shape, bool Weighted>
  void registerNamedMassKernel()
  {
    using Geometry::Polytope;
    for (const auto geometry :
      {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
        Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      for (const bool warm : {false, true})
      {
        const std::string name = "NamedMassKernel/H1/P" + std::to_string(K) + "/" +
          (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
          (Shape == NamedMassShape::Scalar      ? "Scalar/"
              : Shape == NamedMassShape::Vector ? "Vector/"
                                                : "Matrix/") +
          (geometry == Polytope::Type::Point
              ? std::string("Point")
              : std::string(Convergence::UniformGrid::getGeometryName(geometry))) +
          (Weighted ? "/Weighted/" : "/Unit/") + (warm ? "CachedKernel" : "FreshKernel");
        auto* registration = benchmark::RegisterBenchmark(
          name.c_str(), [geometry, warm](benchmark::State& state) {
            if (warm)
              NamedMassKernel<K, Scalar, Shape, Weighted, true>::run(state, geometry);
            else
              NamedMassKernel<K, Scalar, Shape, Weighted, false>::run(state, geometry);
          });
        if constexpr (Shape == NamedMassShape::Scalar)
          registration->Args({1, 1});
        else if constexpr (Shape == NamedMassShape::Vector)
        {
          for (Integer components : {1, 2, 3})
            registration->Args({components, 1});
        }
        else
        {
          registration->Args({1, 1});
          registration->Args({1, 2});
          registration->Args({2, 3});
        }
        registration->UseRealTime();
      }
    }
  }

  /**
   * @brief Registers one degree across scalar fields, shapes and coefficients.
   * @tparam K H1 degree.
   */
  template <size_t K>
  void registerNamedMassDegree()
  {
    registerNamedMassKernel<K, Real, NamedMassShape::Scalar, false>();
    registerNamedMassKernel<K, Real, NamedMassShape::Scalar, true>();
    registerNamedMassKernel<K, Real, NamedMassShape::Vector, false>();
    registerNamedMassKernel<K, Real, NamedMassShape::Vector, true>();
    registerNamedMassKernel<K, Real, NamedMassShape::Matrix, false>();
    registerNamedMassKernel<K, Real, NamedMassShape::Matrix, true>();
    registerNamedMassKernel<K, Complex, NamedMassShape::Scalar, false>();
    registerNamedMassKernel<K, Complex, NamedMassShape::Scalar, true>();
    registerNamedMassKernel<K, Complex, NamedMassShape::Vector, false>();
    registerNamedMassKernel<K, Complex, NamedMassShape::Vector, true>();
    registerNamedMassKernel<K, Complex, NamedMassShape::Matrix, false>();
    registerNamedMassKernel<K, Complex, NamedMassShape::Matrix, true>();
  }
}

/**
 * @brief Runs opt-in named mass-kernel measurements and numerical guards.
 * @param argc Argument count supplied to Google Benchmark.
 * @param argv Argument values supplied to Google Benchmark.
 * @returns Nonzero when arguments or a numerical guard are rejected.
 */
int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerNamedMassDegree<1>();
  registerNamedMassDegree<2>();
  registerNamedMassDegree<3>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return namedMassKernelFailed ? 1 : 0;
}
