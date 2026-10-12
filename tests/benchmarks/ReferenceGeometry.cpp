/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Paired reference-table and pointwise geometry Jacobian measurements.
 */
#include <benchmark/benchmark.h>
#include "Rodin/Geometry/ParametricTransformation.h"
#include "Rodin/Variational/P1/P1Element.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

namespace Rodin::Tests::Benchmarks
{
  bool referenceGeometryFailed = false;

  /**
   * @brief Measures bulk Jacobians against the original pointwise fallback.
   *
   * Both paths evaluate @f$J_{ji}(\widehat x_q)=
   * \sum_a X_{ja}\partial_i\phi_a(\widehat x_q)@f$ with the same vertex
   * basis and accumulation order. The pointwise measurement explicitly calls
   * the base transformation's bulk fallback; no second production path is added.
   * Setup, warmup and exact entry comparisons are outside timing. Repeated
   * calls retain output-vector capacity. Physical Jacobians are timed, but no
   * physical-coordinate evaluation, FE assembly or solve is timed. Results are
   * local diagnostics, not CTest or CI gates.
   *
   * @param state Quadrature order and benchmark iterations.
   * @param geometry Reference-cell family, including a point.
   * @param pointwise Whether to measure the base-class pointwise fallback.
   */
  void referenceGeometry(
    benchmark::State& state, Geometry::Polytope::Type geometry, bool pointwise)
  {
    using namespace Geometry;
    Variational::P1Element<Real> element(geometry);
    PointCloud controls(3, element.getCount());
    for (size_t a = 0; a < element.getCount(); ++a)
    {
      for (size_t j = 0; j < 3; ++j)
        controls(j, a) = Real((a + 1) * (j + 1)) / element.getCount();
    }
    ParametricTransformation transformation(controls, element);
    const auto& formula = QF::PolytopeQuadratureFormula::get(state.range(0), geometry);
    std::vector<Math::SpatialMatrix<Real>> expected, measured;
    transformation.PolytopeTransformation::jacobian(expected, formula);
    const auto evaluate = [&]() {
      if (pointwise)
        transformation.PolytopeTransformation::jacobian(measured, formula);
      else
        transformation.jacobian(measured, formula);
    };
    const auto valid = [&]() {
      if (measured.size() != expected.size())
        return false;
      for (size_t qp = 0; qp < expected.size(); ++qp)
      {
        for (size_t j = 0; j < 3; ++j)
        {
          for (size_t i = 0; i < transformation.getReferenceDimension(); ++i)
          {
            if (measured[qp](j, i) != expected[qp](j, i))
              return false;
          }
        }
      }
      return true;
    };
    evaluate();
    if (!valid())
    {
      referenceGeometryFailed = true;
      state.SkipWithError("Bulk geometry differs from pointwise Jacobians");
      return;
    }
    for (auto _ : state)
    {
      evaluate();
      benchmark::DoNotOptimize(measured.data());
      benchmark::ClobberMemory();
    }
    if (!valid())
    {
      referenceGeometryFailed = true;
      state.SkipWithError("Timed geometry changed Jacobian entries");
    }
    state.counters["dimension"] = transformation.getReferenceDimension();
    state.counters["physical_dimension"] = 3;
    state.counters["basis_functions"] = element.getCount();
    state.counters["quadrature_points"] = formula.getSize();
    state.counters["quadrature_order"] = state.range(0);
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Tests::Benchmarks;
  using G = Polytope::Type;
  for (const auto& [geometry, label] :
    {std::pair{G::Point, "Point"}, {G::Segment, "Segment"}, {G::Triangle, "Triangle"},
      {G::Quadrilateral, "Quadrilateral"}, {G::Tetrahedron, "Tetrahedron"},
      {G::Hexahedron, "Hexahedron"}, {G::Pyramid, "Pyramid"}, {G::Wedge, "Wedge"}})
  {
    for (bool pointwise : {false, true})
    {
      const std::string name =
        std::string("ReferenceGeometry/") + (pointwise ? "Pointwise/" : "Table/") + label;
      benchmark::RegisterBenchmark(name.c_str(),
        [geometry, pointwise](
          auto& state) { referenceGeometry(state, geometry, pointwise); })
        ->Arg(8)
        ->Arg(12)
        ->Arg(20)
        ->UseRealTime();
    }
  }
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || referenceGeometryFailed ? 1 : 0;
}
