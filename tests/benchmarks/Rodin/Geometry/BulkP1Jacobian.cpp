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
#include <type_traits>
#include "Rodin/Geometry.h"
#include "Rodin/Variational/P1.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Tests::Benchmarks
{
  bool bulkP1Failed = false;

  /// @brief Source stages within mapped-quadrature construction.
  enum class P1GeometryOperation
  {
    BulkJacobian, ///< Refill retained Jacobian output.
    OwnedQuadrature, ///< Construct and replace complete mapped quadrature.
    PreparedPoint, ///< Construct one stack point from a prepared Jacobian.
    RetainedPoints ///< Refill retained point storage from prepared Jacobians.
  };

  /**
   * @brief Isolates the bulk geometry operation used by public assembly.
   *
   * ## Architecture
   *
   * The actual native P1 element and its existing reference table are used.
   * Bulk contractions overwrite retained output matrices; OwnedQuadrature
   * constructs/replaces the mapped points used by variational integrators.
   * PreparedPoint constructs and destroys stack points from untimed prepared
   * Jacobians. RetainedPoints clears and refills a reserved point vector from
   * the same inputs, without timed vector allocation or Jacobian evaluation.
   * These scopes expose public source operations, not additive exclusive-time
   * fractions: their storage lifetimes and cache traffic differ.
   * No scope includes mesh construction, reference-table initialization,
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
   * @tparam Operation Source stage to measure.
   */
  template <P1GeometryOperation Operation>
  class BulkP1Jacobian
  {
    public:
      /// @brief Roundoff allowance for affine reproduction, not a mapping predicate.
      static constexpr Real AffineTolerance = 1e-12;

      /**
       * @brief Measures one qualified native P1 geometry operation.
       * @tparam Formula Default dispatcher or directly selected Gaussian family.
       * @param state Embedding dimension, polynomial degree (dispatcher) or
       * points per factor (Gaussian family), and repetition protocol.
       * @param geometry Reference-cell family.
       * @param perturbed Whether control coefficients depart from the affine map.
       */
      template <class Formula = QF::QuadratureFormula>
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
        const QF::QuadratureFormula& formula = [&]() -> const QF::QuadratureFormula& {
          if constexpr (std::is_same_v<Formula, QF::GaussLegendre>)
            return QF::QuadratureFormula::get<Formula>(geometry, size_t(state.range(1)));
          else
            return QF::QuadratureFormula::get(size_t(state.range(1)), geometry);
        }();
        std::vector<Math::SpatialMatrix<Real>> reference(formula.getSize()), result;
        for (size_t qp = 0; qp < formula.getSize(); ++qp)
          transformation.jacobian(reference[qp], formula.getPoint(qp));
        constexpr bool preparedPoints = Operation == P1GeometryOperation::PreparedPoint ||
          Operation == P1GeometryOperation::RetainedPoints;
        std::vector<Math::SpatialMatrix<Real>> prepared;
        std::vector<Point> points;
        if constexpr (preparedPoints)
          transformation.jacobian(prepared, formula);
        if constexpr (Operation == P1GeometryOperation::RetainedPoints)
          points.reserve(prepared.size());
        Optional<PolytopeQuadrature> mapped;
        const auto evaluate = [&]() {
          if constexpr (Operation == P1GeometryOperation::OwnedQuadrature)
          {
            mapped.emplace(cell, formula);
            benchmark::DoNotOptimize(mapped->getSize());
          }
          else if constexpr (Operation == P1GeometryOperation::PreparedPoint)
          {
            for (size_t qp = 0; qp < prepared.size(); ++qp)
            {
              Point point(prepared[qp], cell, formula, qp);
              benchmark::DoNotOptimize(point);
            }
          }
          else if constexpr (Operation == P1GeometryOperation::RetainedPoints)
          {
            points.clear();
            for (size_t qp = 0; qp < prepared.size(); ++qp)
              points.emplace_back(prepared[qp], cell, formula, qp);
            benchmark::DoNotOptimize(points.data());
          }
          else
          {
            transformation.jacobian(result, formula);
            benchmark::DoNotOptimize(result.data());
          }
        };
        bool samplesValid = true;
        const auto extract = [&]() {
          samplesValid = true;
          if constexpr (Operation == P1GeometryOperation::OwnedQuadrature)
          {
            result.resize(mapped->getSize());
            for (size_t qp = 0; qp < result.size(); ++qp)
              result[qp] = mapped->getPoint(qp).getJacobian();
          }
          else if constexpr (preparedPoints)
          {
            result.resize(prepared.size());
            const auto consume = [&](const Point& point, size_t qp) {
              result[qp] = point.getJacobian();
              const auto& actual = point.getReferenceCoordinates();
              const auto& expected = formula.getPoint(qp);
              if (actual.size() != expected.size())
                samplesValid = false;
              else
              {
                for (Eigen::Index i = 0; i < actual.size(); ++i)
                  samplesValid = samplesValid && actual(i) == expected(i);
              }
            };
            for (size_t qp = 0; qp < prepared.size(); ++qp)
            {
              if constexpr (Operation == P1GeometryOperation::PreparedPoint)
              {
                const Point point(prepared[qp], cell, formula, qp);
                consume(point, qp);
              }
              else
                consume(points[qp], qp);
            }
          }
        };
        const auto valid = [&](const auto& matrices) {
          if constexpr (std::is_same_v<Formula, QF::GaussLegendre>)
          {
            const size_t factorPoints = state.range(1);
            if (formula.getSize() != factorPoints * factorPoints * factorPoints)
              return false;
          }
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
        if (!samplesValid || !valid(result) || valid(wrong))
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
        if (!samplesValid || !valid(result))
        {
          bulkP1Failed = true;
          state.SkipWithError("Repeated bulk P1 geometry changed the qualified values");
        }
        state.counters["reference_dimension"] = d;
        state.counters["embedding_dimension"] = s;
        state.counters["basis_functions"] = element.getCount();
        if constexpr (std::is_same_v<Formula, QF::GaussLegendre>)
        {
          state.counters["quadrature_order"] = 2 * state.range(1) - 1;
          state.counters["quadrature_factor_points"] = state.range(1);
        }
        else
          state.counters["quadrature_order"] = state.range(1);
        state.counters["quadrature_points"] = formula.getSize();
        state.counters["contraction_products"] =
          preparedPoints ? 0 : formula.getSize() * element.getCount() * d * s;
        state.counters["point_constructions"] =
          Operation == P1GeometryOperation::BulkJacobian ? 0 : formula.getSize();
        state.counters["prepared_jacobian_count"] = prepared.size();
        state.counters["retained_point_capacity"] = points.capacity();
        state.counters["point_size_bytes"] = sizeof(Point);
        state.counters["output_entries"] = formula.getSize() * d * s;
        state.counters["reference_table_payload_bytes"] =
          formula.getSize() * element.getCount() * (1 + d) * sizeof(Real);
        state.counters["jacobian_output_payload_bytes"] =
          formula.getSize() * sizeof(Math::SpatialMatrix<Real>);
      }

      /**
       * @brief Registers the all-geometry baseline and dispatcher controls.
       *
       * Orders two, eight, fifteen and sixteen cover every admissible embedding.
       * Pyramids additionally retain orders ten and eleven across their table
       * boundary. Tetrahedral Gaussian-factor sweeps isolate sample-count
       * scaling without a change of quadrature family.
       */
      static void registerCases()
      {
        using G = Geometry::Polytope::Type;
        const std::string scope = Operation == P1GeometryOperation::BulkJacobian
          ? "BulkJacobian/"
          : Operation == P1GeometryOperation::OwnedQuadrature ? "OwnedQuadrature/"
          : Operation == P1GeometryOperation::PreparedPoint   ? "PreparedPoint/"
                                                              : "RetainedPoints/";
        for (const auto& [geometry, label] :
          {std::pair{G::Point, "Point"}, {G::Segment, "Segment"},
            {G::Triangle, "Triangle"}, {G::Quadrilateral, "Quadrilateral"},
            {G::Tetrahedron, "Tetrahedron"}, {G::Pyramid, "Pyramid"},
            {G::Hexahedron, "Hexahedron"}, {G::Wedge, "Wedge"}})
        {
          const size_t d = Geometry::Polytope::Traits(geometry).getDimension();
          for (const bool perturbed : {false, true})
          {
            const std::string name = std::string("BulkP1Geometry/") + scope +
              (perturbed ? "PerturbedControls/" : "Affine/") + label;
            auto* registration = benchmark::RegisterBenchmark(name.c_str(),
              [geometry, perturbed](auto& state) { run(state, geometry, perturbed); });
            for (size_t s = std::max(size_t(1), d); s <= 3; ++s)
            {
              for (const int order : {2, 8, 15, 16})
                registration->Args({static_cast<int>(s), order});
              if (geometry == G::Pyramid)
              {
                for (const int order : {10, 11})
                  registration->Args({static_cast<int>(s), order});
              }
            }
            registration->UseRealTime();
          }
        }
        for (const bool perturbed : {false, true})
        {
          const std::string name = std::string("BulkP1Geometry/") + scope +
            "GaussianFactor/" + (perturbed ? "PerturbedControls/" : "Affine/") +
            "Tetrahedron";
          auto* registration =
            benchmark::RegisterBenchmark(name.c_str(), [perturbed](auto& state) {
              run<QF::GaussLegendre>(state, G::Tetrahedron, perturbed);
            });
          for (const int factorPoints : {3, 4, 5, 6, 7, 8, 9, 10, 12, 16})
            registration->Args({3, factorPoints});
          registration->UseRealTime();
        }
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  BulkP1Jacobian<P1GeometryOperation::BulkJacobian>::registerCases();
  BulkP1Jacobian<P1GeometryOperation::OwnedQuadrature>::registerCases();
  BulkP1Jacobian<P1GeometryOperation::PreparedPoint>::registerCases();
  BulkP1Jacobian<P1GeometryOperation::RetainedPoints>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("geometry_element", "Native P1, not H1<1>");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || bulkP1Failed ? 1 : 0;
}
