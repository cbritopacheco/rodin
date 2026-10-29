/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief CPU-time benchmarks for complete cell-local hyperelastic tangents.
 *
 * Mesh, fields, quadrature and reference tabulation are prepared before timing.
 * Each iteration binds one cell and computes all local matrix entries, including
 * constitutive binding, tangent actions and test-gradient contractions. Sparse
 * assembly, solving and matrix reads are excluded. Each case runs on one thread;
 * use the CPU column and run benchmark executables serially. The postcondition
 * checks repeatability; the unit tests supply the component-contraction oracle.
 */
#include <benchmark/benchmark.h>
#include <string>
#include "Rodin/Solid.h"
#include "Rodin/Variational.h"

namespace Rodin::Tests::Benchmarks
{
  bool hyperElasticityKernelFailed = false;

  /**
   * @brief Measures one complete warm tangent matrix in process CPU time.
   * @tparam Degree Degree of the vector displacement, trial and test spaces.
   * @param state Google Benchmark iteration state.
   * @param geometry Full-dimensional cell family.
   * @param dynamic Whether independent active fields and midpoint strain are supplied.
   */
  template <size_t Degree>
  void activeTangent(
    benchmark::State& state, Geometry::Polytope::Type geometry, bool dynamic)
  {
    using namespace Geometry;
    using namespace Variational;
    const size_t dimension = Polytope::Traits(geometry).getDimension();
    auto mesh = LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
    for (size_t i = 1; i <= dimension; ++i)
    {
      for (size_t j = 0; j <= dimension; ++j)
      {
        if (i != j)
        {
          mesh.getConnectivity().compute(i, j);
        }
      }
    }
    using Space = std::conditional_t<Degree == 1, P1<Math::SpatialVector<Real>>,
      H1<Degree, Math::SpatialVector<Real>>>;
    auto space = [&] {
      if constexpr (Degree == 1)
      {
        return Space(mesh, dimension);
      }
      else
      {
        return Space(std::integral_constant<size_t, Degree>{}, mesh, dimension);
      }
    }();
    TrialFunction trial(space);
    TestFunction test(space);
    GridFunction displacement(space);
    for (size_t i = 0; i < space.getSize(); ++i)
    {
      displacement.getData()(i) = 0.001 * (i % 11);
    }
    Solid::ActiveFiberLaw::Parameters parameters;
    parameters.stiffness = 2;
    parameters.damping = 0.3;
    parameters.contractility = 0.8;
    parameters.crossBridgeStiffness = 0.5;
    parameters.initial.stiffness = 1.2;
    parameters.initial.stress = 0.4;
    Solid::ActiveContraction law(
      Solid::NeoHookean(2, 3), Solid::ActiveFiberLaw(parameters));
    Solid::InternalVirtualWork work(law, displacement);
    work.setInput([dynamic](Solid::ConstitutivePoint& point) {
      point.set<Solid::Tags::ActiveExtension>(-0.01);
      if (dynamic)
      {
        point.set<Solid::Tags::TimeStep>(0.02);
        point.set<Solid::Tags::PreviousActiveExtension>(-0.015);
        point.set<Solid::Tags::PreviousActiveGamma>(1.2);
        point.set<Solid::Tags::PreviousActiveBeta>(0.8);
        point.set<Solid::Tags::ElectricalActivation>(0.4);
        point.set<Solid::Tags::ActiveGamma>(1.1);
        point.set<Solid::Tags::ActiveBeta>(0.7);
        point.set<Solid::Tags::PreviousFiberStrain>(0.015);
      }
    });
    auto tangent = work.Tangent(trial, test);
    const auto cell = *mesh.getCell();
    const auto& element = space.getFiniteElement(cell.getDimension(), cell.getIndex());
    const size_t count = element.getCount();
    const size_t order = 2 * element.getDegree().value_or(2 * element.getOrder());
    const auto& formula = QF::QuadratureFormula::get(order, geometry);
    tangent.setPolytope(cell);
    Math::Matrix<Real> expected(count, count);
    for (size_t row = 0; row < count; ++row)
    {
      for (size_t col = 0; col < count; ++col)
      {
        expected(row, col) = tangent.integrate(col, row);
      }
    }
    if (!expected.allFinite() || expected.norm() == 0)
    {
      hyperElasticityKernelFailed = true;
      state.SkipWithError("Initial tangent matrix must be finite and nonzero");
      return;
    }
    for (auto _ : state)
    {
      tangent.setPolytope(cell);
      benchmark::DoNotOptimize(tangent.integrate(0, 0));
      benchmark::ClobberMemory();
    }
    Math::Matrix<Real> result(count, count);
    for (size_t row = 0; row < count; ++row)
    {
      for (size_t col = 0; col < count; ++col)
      {
        result(row, col) = tangent.integrate(col, row);
      }
    }
    if (!result.allFinite() || (result - expected).norm() > 1e-12 * expected.norm())
    {
      hyperElasticityKernelFailed = true;
      state.SkipWithError("Repeated cell binding changed the tangent matrix");
    }
    state.counters["dimension"] = dimension;
    state.counters["factor_degree"] = Degree;
    state.counters["reported_element_order"] = element.getOrder();
    state.counters["requested_quadrature_order"] = order;
    state.counters["quadrature_points"] = formula.getSize();
    state.counters["local_dofs"] = count;
    state.counters["matrix_entries"] = count * count;
    state.counters["cells_per_iteration"] = 1;
    state.counters["relative_repeat_error"] =
      (result - expected).norm() / expected.norm();
  }

  /**
   * @brief Registers serial static and dynamic tangent benchmarks across cell families.
   * @tparam Degree Degree of the displacement and shape-function spaces.
   */
  template <size_t Degree>
  void registerActiveTangent()
  {
    using Geometry::Polytope;
    for (const auto& [geometry, name] : {std::pair{Polytope::Type::Triangle, "Triangle"},
           std::pair{Polytope::Type::Quadrilateral, "Quadrilateral"},
           std::pair{Polytope::Type::Tetrahedron, "Tetrahedron"},
           std::pair{Polytope::Type::Hexahedron, "Hexahedron"}})
    {
      for (bool dynamic : {false, true})
      {
        const std::string label = "InternalVirtualWork/Tangent/Active/P" +
          std::to_string(Degree) + "/" + name + (dynamic ? "/Dynamic" : "/Static");
        benchmark::RegisterBenchmark(label.c_str(),
          [geometry, dynamic](
            benchmark::State& state) { activeTangent<Degree>(state, geometry, dynamic); })
          ->Unit(benchmark::kMicrosecond);
      }
    }
  }

  const bool registeredActiveTangent = [] {
    registerActiveTangent<1>();
    registerActiveTangent<2>();
    registerActiveTangent<3>();
    return true;
  }();
}

int main(int argc, char** argv)
{
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
  {
    return 1;
  }
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return Rodin::Tests::Benchmarks::hyperElasticityKernelFailed ? 1 : 0;
}
