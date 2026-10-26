/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Separates mapped-point caches, geometry kernels and quadrature construction.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include <type_traits>
#include <vector>
#include "Rodin/Geometry.h"
#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/QF/QuadratureFormula.h"
#include "../QF/RepeatedFormula.h"

namespace Rodin::Tests::Benchmarks
{
  bool mappedPointFailed = false;

  enum class PointOperation
  {
    Coordinates,
    Jacobian,
    Inverse,
    Distortion,
    Construction,
    ConstructionAndCoordinates,
    ConstructionAndMetrics,
    ConstructionAndMetricsReferenceMiss,
    BulkJacobian
  };

  /**
   * @brief Measures geometric point operations with an independent regular map.
   *
   * ## Mathematical foundation
   *
   * For reference dimension @f$d@f$ and embedding dimension @f$s\ge d@f$,
   * the map is @f$x_j=r_j+\alpha r_j^2@f$ for @f$j<d@f$ and
   * @f$x_j=(j+1)/4@f$ otherwise. Affine maps have @f$\alpha=0@f$;
   * quadratic maps have @f$\alpha=1/8@f$ and require geometry degree at least
   * two. The Jacobian has diagonal entries @f$1+2\alpha r_j>0@f$ and
   * zero off-diagonal entries. Its left pseudoinverse has reciprocal diagonal
   * entries, and the measure factor is their original diagonal product.
   * A point has an empty Jacobian and measure one; inverse sweeps are omitted.
   * Control coefficients are obtained through the element's DOF functionals,
   * rather than assuming a nodal layout. Analytic map and metric values are
   * compared entrywise with a dimensionless roundoff budget. A perturbed
   * output must be rejected before timing.
   *
   * ## Measurement scopes
   *
   * Rebound scopes clear a preconstructed point's caches through the public
   * setPolytope operation before every evaluation. They include invalidation,
   * the selected lazy computation and output writes, but no point allocation.
   * Warm scopes retain the selected cached quantity. These ordinary points
   * do not carry quadrature provenance: the rebound map/Jacobian is pointwise.
   * Construction measures fresh owned mapped quadrature, including allocation,
   * bulk Jacobians and destruction/replacement, with reference tables warmed.
   * ConstructionAndCoordinates additionally consumes every lazy physical map;
   * ConstructionAndMetrics consumes distortion followed by inverse, matching
   * diffusion binding, while retaining all inverse entries and measure values.
   * A point consumes only its unit measure. These combined scopes are not
   * isolated access latencies. No mesh-owned quadrature
   * lookup, assembly, solver or MPI operation is timed. The separately named
   * ConstructionAndMetricsReferenceMiss scope rotates nine logically distinct
   * copies of the same rule through the eight-entry reference store. It includes
   * reference-table rebuilding as well as mapped construction and metrics.
   * Rule creation, capacity warmup and exact miss-protocol checks are untimed.
   * BulkJacobian calls the transformation's bulk overload with retained matrix
   * output and warmed reference tables. No mapped points or metrics are created;
   * output extraction and analytic qualification remain outside timing.
   * BulkJacobian and ConstructionAndMetrics also sample adjacent dispatcher
   * transition orders. These studies hold the map and geometry degree fixed;
   * they do not replace a required integration degree with a cheaper rule.
   * The reported Jacobian contraction terms count
   * @f$N_qN_\phi ds@f$ scalar control-coordinate/derivative products, not
   * executed floating-point instructions or metric work.
   * Mesh/space creation, analytic guards and sample preparation are untimed.
   *
   * @tparam K Geometry polynomial degree.
   */
  template <size_t K>
  class MappedPoint
  {
    public:
      /// @brief Positive curvature chosen to keep every sampled map regular.
      static constexpr Real Curvature = Real(1) / 8;

      /// @brief Roundoff allowance for independently represented dimensionless maps.
      static constexpr Real OracleTolerance = 1e-9;

      /// @brief One more live rule than the reference store's eight-entry capacity.
      static constexpr size_t ReferenceVersions = 9;

      /**
       * @brief Times a selected geometry scope and checks its complete output.
       * @tparam Formula Default dispatcher or directly selected Gaussian family.
       * @param state Embedding dimension, polynomial order (dispatcher) or points
       * per factor (Gaussian family), and timing iterations.
       * @param geometry Reference polytope family.
       * @param operation Selected cache or construction operation.
       * @param warm Whether point caches remain valid between sweeps.
       * @param curved Whether the analytic map has quadratic curvature.
       */
      template <class Formula = QF::QuadratureFormula>
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        PointOperation operation, bool warm, bool curved)
      {
        using namespace Geometry;
        const Polytope::Traits traits(geometry);
        const size_t d = traits.getDimension(), s = state.range(0);
        const Real alpha = curved ? Curvature : 0;
        const auto coordinate = [&](const Math::SpatialPoint& r, size_t j) {
          return j < d ? r(j) + alpha * r(j) * r(j) : Real(j + 1) / 4;
        };
        Mesh<Context::Local>::Builder builder;
        builder.initialize(s).nodes(traits.getVertexCount());
        IndexArray vertices(traits.getVertexCount());
        for (size_t a = 0; a < traits.getVertexCount(); ++a)
        {
          Math::SpatialPoint vertex(static_cast<std::uint8_t>(s));
          for (size_t j = 0; j < s; ++j)
            vertex(j) = coordinate(traits.getVertex(a), j);
          builder.vertex(vertex);
          vertices(a) = a;
        }
        if (d > 0)
          builder.polytope(geometry, vertices);
        auto mesh = builder.finalize();
        const auto cell = *mesh.getPolytope(d, 0);
        Variational::RealH1Element<K> element(geometry);
        PointCloud controls(s, element.getCount());
        for (size_t a = 0; a < element.getCount(); ++a)
        {
          for (size_t j = 0; j < s; ++j)
          {
            controls(j, a) =
              element.getLinearForm(a)([&](const auto& r) { return coordinate(r, j); });
          }
        }
        mesh.setPolytopeTransformation(
          {d, 0}, new ParametricTransformation(controls, element));
        const QF::QuadratureFormula& formula = [&]() -> const QF::QuadratureFormula& {
          if constexpr (std::is_same_v<Formula, QF::GaussLegendre>)
            return QF::QuadratureFormula::get<Formula>(geometry, size_t(state.range(1)));
          else
            return QF::QuadratureFormula::get(size_t(state.range(1)), geometry);
        }();
        const bool referenceMiss =
          operation == PointOperation::ConstructionAndMetricsReferenceMiss;
        const bool bulk = operation == PointOperation::BulkJacobian;
        const bool metrics =
          operation == PointOperation::ConstructionAndMetrics || referenceMiss;
        std::vector<std::reference_wrapper<const QF::QuadratureFormula>> versions;
        if (referenceMiss)
        {
          versions.reserve(ReferenceVersions);
          for (size_t i = 0; i < ReferenceVersions; ++i)
          {
            versions.emplace_back(
              QF::QuadratureFormula::get<RepeatedFormula>(formula, i));
          }
        }
        size_t nextVersion = 0;
        const QF::QuadratureFormula* activeFormula = &formula;
        std::vector<Point> points;
        for (size_t qp = 0; qp < formula.getSize(); ++qp)
          points.emplace_back(cell, formula.getPoint(qp));
        const bool construction = operation == PointOperation::Construction ||
          operation == PointOperation::ConstructionAndCoordinates || metrics;
        Optional<PolytopeQuadrature> mapped;
        std::vector<Math::SpatialMatrix<Real>> jacobians;
        std::vector<Real> expected;
        for (size_t qp = 0; qp < formula.getSize(); ++qp)
        {
          const auto& r = formula.getPoint(qp);
          if (operation == PointOperation::Coordinates ||
            operation == PointOperation::ConstructionAndCoordinates)
          {
            for (size_t j = 0; j < s; ++j)
              expected.push_back(coordinate(r, j));
          }
          else if (operation == PointOperation::Jacobian ||
            operation == PointOperation::Construction || bulk)
          {
            for (size_t j = 0; j < s; ++j)
            {
              for (size_t i = 0; i < d; ++i)
                expected.push_back(j == i ? 1 + 2 * alpha * r(i) : 0);
            }
          }
          else if (operation == PointOperation::Inverse || metrics)
          {
            for (size_t i = 0; i < d; ++i)
            {
              for (size_t j = 0; j < s; ++j)
                expected.push_back(j == i ? 1 / (1 + 2 * alpha * r(i)) : 0);
            }
            if (metrics)
            {
              Real measure = 1;
              for (size_t i = 0; i < d; ++i)
                measure *= 1 + 2 * alpha * r(i);
              expected.push_back(measure);
            }
          }
          else
          {
            Real measure = 1;
            for (size_t i = 0; i < d; ++i)
              measure *= 1 + 2 * alpha * r(i);
            expected.push_back(measure);
          }
        }
        std::vector<Real> output(expected.size());
        const auto consume = [&](const Point& point, size_t& pos) {
          if (operation == PointOperation::Coordinates ||
            operation == PointOperation::ConstructionAndCoordinates)
          {
            const auto& x = point.getPhysicalCoordinates();
            for (size_t j = 0; j < s; ++j)
              output[pos++] = x(j);
          }
          else if (operation == PointOperation::Jacobian ||
            operation == PointOperation::Construction)
          {
            const auto& jacobian = point.getJacobian();
            for (size_t j = 0; j < s; ++j)
            {
              for (size_t i = 0; i < d; ++i)
                output[pos++] = jacobian(j, i);
            }
          }
          else if (metrics)
          {
            const Real measure = point.getDistortion();
            if (d > 0)
            {
              const auto& inverse = point.getJacobianInverse();
              for (size_t i = 0; i < d; ++i)
              {
                for (size_t j = 0; j < s; ++j)
                  output[pos++] = inverse(i, j);
              }
            }
            output[pos++] = measure;
          }
          else if (operation == PointOperation::Inverse)
          {
            const auto& inverse = point.getJacobianInverse();
            for (size_t i = 0; i < d; ++i)
            {
              for (size_t j = 0; j < s; ++j)
                output[pos++] = inverse(i, j);
            }
          }
          else
            output[pos++] = point.getDistortion();
        };
        const auto evaluate = [&]() {
          size_t pos = 0;
          if (bulk)
          {
            cell.getTransformation().jacobian(jacobians, formula);
            benchmark::DoNotOptimize(jacobians.data());
          }
          else if (construction)
          {
            activeFormula = referenceMiss ? &versions[nextVersion].get() : &formula;
            mapped.emplace(cell, *activeFormula);
            if (referenceMiss)
              nextVersion = (nextVersion + 1) % versions.size();
            if (operation == PointOperation::ConstructionAndCoordinates || metrics)
            {
              for (size_t qp = 0; qp < mapped->getSize(); ++qp)
                consume(mapped->getPoint(qp), pos);
            }
            benchmark::DoNotOptimize(mapped->getSize());
          }
          else
          {
            for (auto& point : points)
            {
              if (!warm)
                point.setPolytope(cell);
              consume(point, pos);
            }
          }
        };
        const auto outputValid = [&]() {
          for (size_t i = 0; i < expected.size(); ++i)
          {
            if (!std::isfinite(output[i]) ||
              std::abs(output[i] - expected[i]) >
                OracleTolerance * std::max(Real(1), std::abs(expected[i])))
              return false;
          }
          return true;
        };
        const auto valid = [&]() {
          if constexpr (std::is_same_v<Formula, QF::GaussLegendre>)
          {
            const size_t factorPoints = state.range(1);
            if (formula.getSize() != factorPoints * factorPoints * factorPoints)
              return false;
          }
          for (const auto& point : points)
          {
            if (point.getDimension(PointBase::Coordinates::Reference) != d ||
              point.getDimension(PointBase::Coordinates::Physical) != s)
              return false;
          }
          if (bulk)
          {
            if (jacobians.size() != formula.getSize())
              return false;
            size_t pos = 0;
            for (const auto& jacobian : jacobians)
            {
              if (jacobian.rows() != s || jacobian.cols() != d)
                return false;
              for (size_t j = 0; j < s; ++j)
              {
                for (size_t i = 0; i < d; ++i)
                  output[pos++] = jacobian(j, i);
              }
            }
          }
          else if (construction)
          {
            if (!mapped || mapped->getSize() != formula.getSize())
              return false;
            if (&mapped->getQuadratureFormula() != activeFormula)
              return false;
            for (size_t qp = 0; qp < mapped->getSize(); ++qp)
            {
              const auto& point = mapped->getPoint(qp);
              if (point.getReferenceCoordinates().size() != d ||
                point.getJacobian().rows() != s || point.getJacobian().cols() != d)
                return false;
              if (metrics && d > 0 &&
                (point.getJacobianInverse().rows() != d ||
                  point.getJacobianInverse().cols() != s))
                return false;
            }
            if (operation == PointOperation::Construction)
            {
              size_t pos = 0;
              for (size_t qp = 0; qp < mapped->getSize(); ++qp)
                consume(mapped->getPoint(qp), pos);
            }
          }
          return outputValid();
        };
        const auto fail = [&](const char* message) {
          mappedPointFailed = true;
          state.SkipWithError(message);
        };
        evaluate();
        if (!valid())
        {
          fail("Mapped geometry disagrees with independent regular-map oracle");
          return;
        }
        // Validate sensitivity after extracting construction outputs, so the
        // extraction itself cannot overwrite the deliberate perturbation.
        if (!output.empty())
        {
          output.front() += 1;
          if (outputValid())
          {
            fail("Geometry oracle accepted a perturbed output");
            return;
          }
        }
        evaluate();
        if (referenceMiss)
        {
          for (const auto& version : versions)
            element.getTabulation(version.get());
          for (const auto& version : versions)
          {
            if (element.findTabulation(version.get()))
            {
              fail("Reference-miss protocol found a resident table before rebuilding");
              return;
            }
            element.getTabulation(version.get());
          }
          nextVersion = 0;
        }
        for (auto _ : state)
        {
          evaluate();
          benchmark::DoNotOptimize(output.data());
          benchmark::ClobberMemory();
        }
        if (!valid())
          fail("Timed geometric operation changed analytic values");
        if (referenceMiss && element.findTabulation(versions[nextVersion].get()))
          fail("Reference-miss cycle unexpectedly retained the next table");
        state.counters["reference_dimension"] = d;
        state.counters["cells"] = 1;
        state.counters["embedding_dimension"] = s;
        state.counters["geometry_degree"] = K;
        state.counters["curvature"] = alpha;
        if constexpr (std::is_same_v<Formula, QF::GaussLegendre>)
        {
          state.counters["quadrature_order"] = 2 * state.range(1) - 1;
          state.counters["quadrature_factor_points"] = state.range(1);
        }
        else
          state.counters["quadrature_order"] = state.range(1);
        state.counters["quadrature_points"] = formula.getSize();
        state.counters["basis_functions"] = element.getCount();
        state.counters["output_entries"] = output.size();
        state.counters["formula_identities"] = referenceMiss ? versions.size() : 1;
        state.counters["reference_table_miss"] = referenceMiss;
        if (bulk || construction)
        {
          state.counters["jacobian_contraction_terms"] =
            formula.getSize() * element.getCount() * d * s;
          state.counters["reference_table_payload_bytes"] =
            formula.getSize() * element.getCount() * (1 + d) * sizeof(Real);
          state.counters["jacobian_output_payload_bytes"] =
            formula.getSize() * sizeof(Math::SpatialMatrix<Real>);
        }
      }

      /**
       * @brief Registers all families, embeddings and quadrature-order studies.
       *
       * Every scope retains orders two, four and eight. Bulk Jacobians and
       * mapped construction with metrics additionally cover orders fifteen
       * through seventeen, with adjacent compact/product transitions on
       * tetrahedra, pyramids, triangles and wedges. Other geometries provide
       * Gaussian-product and zero-dimensional controls at the same orders.
       * Thirty additional tetrahedral bulk scopes select only the Gaussian
       * family. Affine degree-one/two geometry and quadratic degree-two geometry
       * hold each map and basis fixed while varying sample count. Their second
       * argument is the one-dimensional factor count.
       */
      static void registerCases()
      {
        using G = Geometry::Polytope::Type;
        for (const auto& [geometry, label] :
          {std::pair{G::Point, "Point"}, {G::Segment, "Segment"},
            {G::Triangle, "Triangle"}, {G::Quadrilateral, "Quadrilateral"},
            {G::Tetrahedron, "Tetrahedron"}, {G::Pyramid, "Pyramid"},
            {G::Hexahedron, "Hexahedron"}, {G::Wedge, "Wedge"}})
        {
          const size_t d = Geometry::Polytope::Traits(geometry).getDimension();
          for (const bool curved : {false, true})
          {
            if (curved && (K < 2 || d == 0))
              continue;
            for (const auto& [operation, operationLabel] :
              {std::pair{PointOperation::Coordinates, "Coordinates"},
                {PointOperation::Jacobian, "Jacobian"},
                {PointOperation::Inverse, "Inverse"},
                {PointOperation::Distortion, "Distortion"},
                {PointOperation::Construction, "Construction"},
                {PointOperation::ConstructionAndCoordinates,
                  "ConstructionAndCoordinates"},
                {PointOperation::ConstructionAndMetrics, "ConstructionAndMetrics"},
                {PointOperation::ConstructionAndMetricsReferenceMiss,
                  "ConstructionAndMetricsReferenceMiss"},
                {PointOperation::BulkJacobian, "BulkJacobian"}})
            {
              if (d == 0 && operation == PointOperation::Inverse)
                continue;
              for (const bool warm : {false, true})
              {
                const bool construction = operation == PointOperation::Construction ||
                  operation == PointOperation::ConstructionAndCoordinates ||
                  operation == PointOperation::ConstructionAndMetrics ||
                  operation == PointOperation::ConstructionAndMetricsReferenceMiss ||
                  operation == PointOperation::BulkJacobian;
                if (construction && warm)
                  continue;
                const std::string name = "MappedPoint/" + std::string(operationLabel) +
                  (construction ? "/"
                      : warm    ? "/Warm/"
                                : "/Rebound/") +
                  (curved ? "Quadratic/" : "Affine/") + "P" + std::to_string(K) + "/" +
                  label;
                auto* registration = benchmark::RegisterBenchmark(
                  name.c_str(), [geometry, operation, warm, curved](auto& state) {
                    run(state, geometry, operation, warm, curved);
                  });
                std::vector<int> orders{2, 4, 8};
                if (operation == PointOperation::BulkJacobian ||
                  operation == PointOperation::ConstructionAndMetrics)
                {
                  orders.insert(orders.end(), {15, 16, 17});
                  if (geometry == G::Tetrahedron)
                    orders.push_back(14);
                  else if (geometry == G::Pyramid)
                    orders.insert(orders.end(), {9, 10, 11});
                  else if (geometry == G::Triangle || geometry == G::Wedge)
                    orders.insert(orders.end(), {49, 50, 51});
                }
                for (size_t s = std::max(size_t(1), d); s <= 3; ++s)
                {
                  for (const int order : orders)
                    registration->Args({static_cast<int>(s), order});
                }
                registration->UseRealTime();
              }
            }
          }
        }
        if constexpr (K <= 2)
        {
          for (const bool curved : {false, true})
          {
            if (curved && K < 2)
              continue;
            const std::string name = "MappedPoint/BulkJacobian/GaussianFactor/" +
              std::string(curved ? "Quadratic/" : "Affine/") + "P" + std::to_string(K) +
              "/Tetrahedron";
            auto* registration =
              benchmark::RegisterBenchmark(name.c_str(), [curved](auto& state) {
                run<QF::GaussLegendre>(
                  state, G::Tetrahedron, PointOperation::BulkJacobian, false, curved);
              });
            for (const int factorPoints : {3, 4, 5, 6, 7, 8, 9, 10, 12, 16})
              registration->Args({3, factorPoints});
            registration->UseRealTime();
          }
        }
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  MappedPoint<1>::registerCases();
  MappedPoint<2>::registerCases();
  MappedPoint<3>::registerCases();
  MappedPoint<4>::registerCases();
  MappedPoint<5>::registerCases();
  MappedPoint<6>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || mappedPointFailed ? 1 : 0;
}
