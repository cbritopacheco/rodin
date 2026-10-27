/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Cell-local named elasticity kernels with explicit quadrature accounting.
 *
 * The benchmark belongs to the Solid module. Sparse assembly and solving are
 * deliberately excluded; their owning modules provide separate consumers.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include <string>
#include <type_traits>

#include "Rodin/Solid.h"
#include "convergence/CurvedGeometry.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  bool elasticityKernelFailed = false;

  /**
   * @brief Measures one local stiffness matrix without scatter or global setup.
   *
   * ## Measurement contract
   *
   * Warm retains the kernel, output capacity, reference derivatives and mapped
   * geometry. FirstKernel constructs, evaluates and destroys a kernel and its
   * output on every iteration, with mapped geometry already warm. Neither scope
   * represents process-first quadrature lookup or cold geometry. Exactly one
   * selected cell is evaluated per iteration; mesh and space construction,
   * oracle evaluation and guards remain outside timing.
   *
   * Native and prescribed quadratic geometry use the same unit-box mesh and
   * full-dimensional vector space. Constant real Lamé parameters are fixed at
   * @f$\lambda=2@f$ and @f$\mu=3@f$. The public integral supplies an independent
   * implementation of the same requested quadrature rule, not a claim of exact
   * integration on curved or rational cells. The dimensionless comparison is
   * @f$\|A_K-A_K^{\mathrm{integral}}\|_F/
   * \max(1,\|A_K^{\mathrm{integral}}\|_F)\le10^{-11}@f$.
   * A perturbed matrix must be rejected before timing starts.
   *
   * Counters distinguish family order, polynomial degree when known, requested
   * quadrature order and actual point count. The pair-component work proxy is
   * @f$n_q n_{\mathrm{local}}^2 s^2@f$; it is not an instruction count and
   * excludes reference preparation and geometry costs. No source optimization,
   * quadrature reduction, new cache, collective or CTest timing gate is added.
   *
   * @tparam Degree Native P1 or H1 degree two.
   * @tparam Curved Whether the prescribed quadratic map is installed.
   * @tparam First Whether kernel construction and destruction are timed.
   * @param state Google Benchmark iteration state.
   * @param geometry Full-dimensional volume-cell family.
   */
  template <size_t Degree, bool Curved, bool First>
  void elasticityKernel(benchmark::State& state, Polytope::Type geometry)
  {
    using Range = Math::SpatialVector<Real>;
    using Space = std::conditional_t<Degree == 1, P1<Range>, H1<Degree, Range>>;
    auto mesh = Convergence::UniformGrid(geometry).makeMesh(2);
    Convergence::CurvedGeometry<LocalMesh> mapping(mesh);
    if constexpr (Curved)
      mapping.template install<2>();
    const size_t dimension = mesh.getSpaceDimension();
    auto space = [&] {
      if constexpr (Degree == 1)
        return Space(mesh, dimension);
      else
        return Space(std::integral_constant<size_t, Degree>{}, mesh, dimension);
    }();
    TrialFunction trial(space);
    TestFunction test(space);
    const RealFunction lambda(Real(2)), mu(Real(3));
    using Form = decltype(LinearElasticityForm(lambda, mu, trial, test));
    using Kernel = typename Form::Kernel;
    using Matrix = typename Kernel::MatrixType;
    const auto cell = *mesh.getCell();
    const auto& element = space.getFiniteElement(cell.getDimension(), cell.getIndex());
    const size_t count = element.getCount();
    const size_t order = 2 * element.getDegree().value_or(2 * element.getOrder());
    const auto& formula = QF::QuadratureFormula::get(order, geometry);
    auto integral = LinearElasticityIntegral(trial, test)(lambda, mu);
    integral.setPolytope(cell);
    Matrix expected(count, count);
    for (size_t row = 0; row < count; ++row)
    {
      for (size_t col = 0; col < count; ++col)
        expected(row, col) = integral.integrate(col, row);
    }
    constexpr Real OperatorTolerance = 1e-11;
    const Real scale = std::max(Real(1), expected.norm());
    const auto valid = [&](const Matrix& matrix) {
      return matrix.rows() == expected.rows() && matrix.cols() == expected.cols() &&
        matrix.allFinite() && (matrix - expected).norm() <= OperatorTolerance * scale;
    };
    Kernel kernel(space, space, lambda, mu);
    Matrix result;
    kernel.compute(result, cell);
    Matrix wrong = result;
    wrong(0, 0) += scale;
    if (!expected.allFinite() || expected.norm() == 0 || !valid(result) || valid(wrong))
    {
      elasticityKernelFailed = true;
      state.SkipWithError("Local elasticity operator or sensitivity guard failed");
      return;
    }
    for (auto _ : state)
    {
      if constexpr (First)
      {
        Kernel first(space, space, lambda, mu);
        Matrix local;
        first.compute(local, cell);
        benchmark::DoNotOptimize(local.data());
        benchmark::ClobberMemory();
        // Keep the latest result for an untimed postcondition without copying.
        result.swap(local);
      }
      else
      {
        kernel.compute(result, cell);
        benchmark::DoNotOptimize(result.data());
        benchmark::ClobberMemory();
      }
    }
    if (!valid(result))
    {
      elasticityKernelFailed = true;
      state.SkipWithError("Repeated local elasticity evaluation changed the operator");
    }
    state.counters["dimension"] = dimension;
    state.counters["family_order"] = element.getOrder();
    if (const auto degree = element.getDegree())
      state.counters["reference_total_degree"] = *degree;
    state.counters["requested_quadrature_order"] = order;
    state.counters["quadrature_points"] = formula.getSize();
    state.counters["local_dofs"] = count;
    state.counters["matrix_entries"] = count * count;
    state.counters["pair_component_proxy"] =
      formula.getSize() * count * count * dimension * dimension;
    state.counters["global_dofs"] = space.getSize();
    state.counters["mesh_cells"] = mesh.getCellCount();
    state.counters["cells_per_iteration"] = 1;
    state.counters["geometry_degree"] = Curved ? 2 : 1;
    state.counters["relative_error"] = (result - expected).norm() / scale;
  }

  template <size_t Degree, bool Curved, bool First>
  void registerElasticityKernel()
  {
    const std::string prefix = "LinearElasticity/CellKernel/P" + std::to_string(Degree) +
      (Curved ? "/CurvedQ2/" : "/Native/") + (First ? "FirstKernel/" : "Warm/");
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Hexahedron, Polytope::Type::Wedge, Polytope::Type::Pyramid})
    {
      const std::string name =
        prefix + std::string(Convergence::UniformGrid::getGeometryName(geometry));
      benchmark::RegisterBenchmark(name.c_str(), [geometry](benchmark::State& state) {
        elasticityKernel<Degree, Curved, First>(state, geometry);
      })->Unit(benchmark::kMicrosecond);
    }
  }

  const bool registeredElasticityKernel = [] {
    registerElasticityKernel<1, false, false>();
    registerElasticityKernel<1, false, true>();
    registerElasticityKernel<1, true, false>();
    registerElasticityKernel<1, true, true>();
    registerElasticityKernel<2, false, false>();
    registerElasticityKernel<2, false, true>();
    registerElasticityKernel<2, true, false>();
    registerElasticityKernel<2, true, true>();
    return true;
  }();
}

int main(int argc, char** argv)
{
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return Rodin::Tests::Benchmarks::elasticityKernelFailed ? 1 : 0;
}
