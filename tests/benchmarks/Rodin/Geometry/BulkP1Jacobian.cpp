/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Bulk P1 geometry contractions and their mapped-quadrature consumer.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include "Rodin/Geometry.h"
#include "Rodin/Variational/P1.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

namespace Rodin::Tests::Benchmarks
{
  bool bulkP1Failed = false;

  /**
   * @brief Isolates the bulk geometry operation used by public assembly.
   *
   * ## Architecture
   *
   * The actual native P1 element and its existing reference table are used.
   * Bulk contractions overwrite retained output matrices; OwnedQuadrature
   * constructs/replaces the mapped points used by variational integrators.
   * Neither scope includes mesh construction, reference-table initialization,
   * analytic guards, operator assembly or solving. Full assembly is measured
   * separately by its owning Assembly targets on matched quadrature orders.
   *
   * Affine controls represent @f$x=b+A\widehat x@f$, with nonzero mixed
   * components, so the independent Jacobian is @f$A@f$. Perturbed controls
   * exercise the same native basis without assuming an affine Jacobian on
   * tensor-product, wedge or rational pyramid cells. Simplex P1 maps remain
   * affine. Every result must exactly reproduce the original pointwise
   * contraction and its dimensions; affine cases additionally satisfy an
   * independent analytic roundoff budget. Wrong entries or dimensions must
   * be rejected. No new cache or alternative production geometry path is added.
   *
   * @tparam OwnedQuadrature Whether mapped-point construction is measured.
   */
  template <bool OwnedQuadrature>
  class BulkP1Jacobian
  {
    public:
    /// @brief Roundoff allowance for affine reproduction, not a mapping predicate.
    static constexpr Real AffineTolerance = 1e-12;

    /**
     * @brief Measures one qualified native P1 geometry operation.
     * @param state Embedding dimension, quadrature order and repetition protocol.
     * @param geometry Reference-cell family.
     * @param perturbed Whether control coefficients depart from the affine map.
     */
    static void run(
      benchmark::State& state, Geometry::Polytope::Type geometry, bool perturbed)
    {
      using namespace Geometry;
      const Polytope::Traits traits(geometry);
      const size_t d = traits.getDimension(), s = state.range(0);
      const Variational::RealP1Element element(geometry);
      Math::SpatialMatrix<Real> affine(s, d);
      for (size_t j = 0; j < s; ++j)
      {
        for (size_t i = 0; i < d; ++i)
        {
          affine(j, i) = j == i ? Real(1) + Real(i + 1) / 8
            : (j + i) % 2       ? -Real(1) / 32
                                : Real(1) / 16;
        }
      }
      PointCloud controls(s, element.getCount());
      for (size_t a = 0; a < element.getCount(); ++a)
      {
        for (size_t j = 0; j < s; ++j)
        {
          controls(j, a) = element.getLinearForm(a)([&](const auto& reference) {
            Real value = Real(j + 1) / 8;
            for (size_t i = 0; i < d; ++i)
              value += affine(j, i) * reference(i);
            return value;
          });
          if (perturbed)
            controls(j, a) += Real(Integer((a + j) % 3) - 1) / 64;
        }
      }
      Mesh<Context::Local>::Builder builder;
      builder.initialize(s).nodes(element.getCount());
      IndexArray vertices(element.getCount());
      for (size_t a = 0; a < element.getCount(); ++a)
      {
        Math::SpatialPoint vertex(s);
        for (size_t j = 0; j < s; ++j)
          vertex(j) = controls(j, a);
        builder.vertex(vertex);
        vertices(a) = a;
      }
      if (d > 0)
        builder.polytope(geometry, vertices);
      auto mesh = builder.finalize();
      mesh.setPolytopeTransformation(
        {d, 0}, new ParametricTransformation(controls, element));
      const auto cell = *mesh.getPolytope(d, 0);
      const auto& transformation = cell.getTransformation();
      const QF::PolytopeQuadratureFormula formula(state.range(1), geometry);
      std::vector<Math::SpatialMatrix<Real>> reference(formula.getSize()), result;
      for (size_t qp = 0; qp < formula.getSize(); ++qp)
        transformation.jacobian(reference[qp], formula.getPoint(qp));
      Optional<PolytopeQuadrature> mapped;
      const auto evaluate = [&]() {
        if constexpr (OwnedQuadrature)
        {
          mapped.emplace(cell, formula);
          benchmark::DoNotOptimize(mapped->getSize());
        }
        else
        {
          transformation.jacobian(result, formula);
          benchmark::DoNotOptimize(result.data());
        }
      };
      const auto extract = [&]() {
        if constexpr (OwnedQuadrature)
        {
          result.resize(mapped->getSize());
          for (size_t qp = 0; qp < result.size(); ++qp)
            result[qp] = mapped->getPoint(qp).getJacobian();
        }
      };
      const auto valid = [&](const auto& matrices) {
        if (matrices.size() != reference.size())
          return false;
        for (size_t qp = 0; qp < matrices.size(); ++qp)
        {
          const auto& matrix = matrices[qp];
          if (matrix.rows() != s || matrix.cols() != d)
            return false;
          for (size_t j = 0; j < s; ++j)
          {
            for (size_t i = 0; i < d; ++i)
            {
              if (!std::isfinite(matrix(j, i)) || matrix(j, i) != reference[qp](j, i) ||
                (!perturbed && std::abs(matrix(j, i) - affine(j, i)) > AffineTolerance))
                return false;
            }
          }
        }
        return true;
      };
      evaluate();
      extract();
      auto wrong = result;
      if (d > 0)
        wrong.front()(0, 0) += 1;
      else
        wrong.front().resize(s, 1);
      if (!valid(result) || valid(wrong))
      {
        bulkP1Failed = true;
        state.SkipWithError("Initial bulk P1 geometry or sensitivity guard failed");
        return;
      }
      for (auto _ : state)
      {
        evaluate();
        benchmark::ClobberMemory();
      }
      extract();
      if (!valid(result))
      {
        bulkP1Failed = true;
        state.SkipWithError("Repeated bulk P1 geometry changed the qualified values");
      }
      state.counters["reference_dimension"] = d;
      state.counters["embedding_dimension"] = s;
      state.counters["basis_functions"] = element.getCount();
      state.counters["quadrature_order"] = state.range(1);
      state.counters["quadrature_points"] = formula.getSize();
      state.counters["contraction_products"] =
        formula.getSize() * element.getCount() * d * s;
      state.counters["output_entries"] = formula.getSize() * d * s;
    }

    /// @brief Registers all geometries, admissible embeddings and three orders.
    static void registerCases()
    {
      using G = Geometry::Polytope::Type;
      for (const auto& [geometry, label] :
        {std::pair{G::Point, "Point"}, {G::Segment, "Segment"}, {G::Triangle, "Triangle"},
          {G::Quadrilateral, "Quadrilateral"}, {G::Tetrahedron, "Tetrahedron"},
          {G::Pyramid, "Pyramid"}, {G::Hexahedron, "Hexahedron"}, {G::Wedge, "Wedge"}})
      {
        const size_t d = Geometry::Polytope::Traits(geometry).getDimension();
        for (const bool perturbed : {false, true})
        {
          const std::string name = std::string("BulkP1Geometry/") +
            (OwnedQuadrature ? "OwnedQuadrature/" : "BulkJacobian/") +
            (perturbed ? "PerturbedControls/" : "Affine/") + label;
          auto* registration = benchmark::RegisterBenchmark(name.c_str(),
            [geometry, perturbed](auto& state) { run(state, geometry, perturbed); });
          for (size_t s = std::max(size_t(1), d); s <= 3; ++s)
          {
            for (const int order : {2, 8, 16})
              registration->Args({static_cast<int>(s), order});
          }
          registration->UseRealTime();
        }
      }
    }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  BulkP1Jacobian<false>::registerCases();
  BulkP1Jacobian<true>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("geometry_element", "Native P1, not H1<1>");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || bulkP1Failed ? 1 : 0;
}
