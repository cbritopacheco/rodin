/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_BENCHMARKS_H1_BASISLIFECYCLE_H
#define RODIN_TESTS_BENCHMARKS_H1_BASISLIFECYCLE_H

#include <benchmark/benchmark.h>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include "Rodin/Variational/H1.h"
#include "../../../../convergence/Convergence.h"

namespace Rodin::Tests::Benchmarks
{
  inline bool h1BasisLifecycleFailed = false;

  /// @brief Reference authority and formula-identity working-set protocol.
  enum class H1BasisIdentityProtocol
  {
    AssociatedFresh,
    AssociatedOne,
    AssociatedNine,
    UnprovenNine
  };

  /**
   * @brief Measures reference association, identity churn and sample-sweep costs.
   *
   * The H1 component basis is @f$\Phi_{a,c}=\phi_a E_c@f$. Scalar and vector
   * controls retain the same scalar basis and component convention. Associated
   * Single scopes bind sample zero; Sweep scopes bind every mapped sample.
   * AssociatedFresh measures one first binding before any numerical warmup;
   * its direct checks run afterwards. AssociatedOne repeats one owned formula.
   * AssociatedNine cycles nine
   * distinct formula identities, targeting the current eight-entry reference
   * store. UnprovenNine uses the same owned points but only origin metadata,
   * which must never authorize reference-table substitution. Unproven Single
   * binds the last sample while origin metadata remains index zero; Unproven
   * Sweep likewise retains origin index zero while actual coordinates vary.
   *
   * ## Measurement boundary
   *
   * Mesh, space, formulas and mapped samples are constructed outside timing.
   * Except for AssociatedFresh, two untimed priming cycles establish the identity
   * protocol. Fresh scopes force one iteration and one repetition; a controlled
   * process-first study must additionally select one case per new process.
   * The timed
   * operation includes identity selection, binding and all active-value stores.
   * Identity churn can include allocation, value and derivative tabulation for
   * implementations using the reference store. It is not process-first timing;
   * shared element initialization is already warm. No query-count instrumentation
   * is inserted into timed formula access. Public table-presence queries occur
   * only between priming cycles and report participation rather than requiring
   * a particular implementation to use the cache.
   *
   * ## Numerical qualification
   *
   * Every output is checked exactly against direct scalar basis evaluation.
   * Logical dimensions and reference authority are checked outside timing;
   * Fresh scopes perform all numerical qualification after the first timed bind.
   * Changed output must fail. Only one active-value output buffer is retained.
   *
   * @tparam K H1 degree.
   * @tparam Range Real/complex scalar, vector or matrix field range.
   */
  template <size_t K, class Range>
  class H1BasisLifecycle
  {
    public:
      /// @brief Scalar field used by basis values and the direct oracle.
      using Scalar = typename FormLanguage::Traits<Range>::ScalarType;

      /**
       * @brief Measures a qualified binding protocol on one cell.
       * @param state Requested quadrature order and iteration state.
       * @param geometry Reference-cell family, including Point.
       * @param protocol Reference authority and formula-identity working set.
       * @param whole Whether the full quadrature sample set is traversed.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        H1BasisIdentityProtocol protocol, bool whole)
      {
        using namespace Geometry;
        using namespace Variational;
        // Protocol size, not a public cache-capacity guarantee.
        constexpr size_t ChurnIdentities = 9;
        constexpr size_t Rows = FormLanguage::IsMatrixRange<Range>::Value ? 2
          : FormLanguage::IsVectorRange<Range>::Value                     ? 3
                                                                          : 1;
        constexpr size_t Columns = FormLanguage::IsMatrixRange<Range>::Value ? 3 : 1;
        constexpr size_t Components = Rows * Columns;
        const bool fresh = protocol == H1BasisIdentityProtocol::AssociatedFresh;
        const bool associated = protocol != H1BasisIdentityProtocol::UnprovenNine;
        const bool churn = protocol == H1BasisIdentityProtocol::AssociatedNine ||
          protocol == H1BasisIdentityProtocol::UnprovenNine;
        const size_t order = state.range(0), identities = churn ? ChurnIdentities : 1;
        auto mesh = [&]() {
          if (geometry == Polytope::Type::Point)
            return LocalMesh::Builder()
              .initialize(3)
              .nodes(1)
              .vertex({0, 0, 0})
              .finalize();
          return Convergence::UniformGrid(geometry).makeMesh(2);
        }();
        for (size_t d = 1; d <= mesh.getDimension(); ++d)
        {
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
        using Space = H1<K, Range, decltype(mesh)>;
        auto space = [&]() {
          if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
            return Space(std::integral_constant<size_t, K>{}, mesh, Rows, Columns);
          else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
            return Space(std::integral_constant<size_t, K>{}, mesh, Rows);
          else
            return Space(std::integral_constant<size_t, K>{}, mesh);
        }();
        const auto cell = *mesh.getCell();
        const H1Element<K, Scalar> scalar(geometry);
        const auto& prototype = QF::PolytopeQuadratureFormula::get(order, geometry);
        std::vector<std::unique_ptr<const QF::QuadratureFormulaBase>> formulas;
        std::vector<std::unique_ptr<PolytopeQuadrature>> quadratures;
        formulas.reserve(identities);
        quadratures.reserve(identities);
        for (size_t i = 0; i < identities; ++i)
        {
          formulas.emplace_back(prototype.copy());
          quadratures.emplace_back(
            std::make_unique<PolytopeQuadrature>(cell, *formulas.back()));
        }
        IntegrationPoint cursor(*quadratures.front(), 0);
        TrialFunction u(space);
        auto& leaf = u.getDerived();
        const size_t count = leaf.getDOFs(cell);
        const size_t samples = whole ? prototype.getSize() : 1;
        const size_t singleSample = associated ? 0 : prototype.getSize() - 1;
        std::vector<Scalar> output(samples * count * Components);
        const auto fail = [&](const char* message) {
          h1BasisLifecycleFailed = true;
          state.SkipWithError(message);
        };
        if (count != scalar.getCount() * Components || prototype.getSize() == 0)
        {
          fail("Basis dimensions or quadrature sample count are invalid");
          return;
        }
        for (size_t i = 0; i < identities; ++i)
        {
          for (size_t j = 0; j < i; ++j)
          {
            if (formulas[i]->getCacheIdentity() == formulas[j]->getCacheIdentity())
            {
              fail("Formula copies did not receive distinct logical identities");
              return;
            }
          }
        }
        const auto component = [](const auto& value, size_t c) -> Scalar {
          if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
            return value(c / Columns, c % Columns);
          else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
            return value(c);
          else
            return value;
        };
        const auto sweep = [&]<bool Check>(size_t identity) {
          const auto& quadrature = *quadratures[identity];
          for (size_t qp = 0; qp < samples; ++qp)
          {
            const size_t sample = whole ? qp : singleSample;
            cursor = associated ? IntegrationPoint(quadrature, sample)
                                : IntegrationPoint(quadrature.getPoint(sample),
                                    formulas[identity].get(), 0);
            if constexpr (Check)
            {
              if (cursor.getReferenceQuadratureFormula() !=
                (associated ? formulas[identity].get() : nullptr))
                return false;
            }
            leaf.setIntegrationPoint(cursor);
            for (size_t a = 0; a < count; ++a)
            {
              decltype(auto) value = leaf.getBasis(a);
              if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
              {
                if constexpr (Check)
                {
                  if (value.rows() != Rows || value.cols() != Columns)
                    return false;
                }
              }
              else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
              {
                if constexpr (Check)
                {
                  if (value.size() != Rows)
                    return false;
                }
              }
              for (size_t c = 0; c < Components; ++c)
                output[(qp * count + a) * Components + c] = component(value, c);
            }
          }
          return true;
        };
        const auto valid = [&](size_t identity) {
          for (size_t qp = 0; qp < samples; ++qp)
          {
            const size_t sample = whole ? qp : singleSample;
            const auto& reference =
              quadratures[identity]->getPoint(sample).getReferenceCoordinates();
            for (size_t a = 0; a < scalar.getCount(); ++a)
            {
              const Scalar phi = scalar.getBasis(a)(reference);
              for (size_t c = 0; c < Components; ++c)
              {
                for (size_t entry = 0; entry < Components; ++entry)
                {
                  if (output[(qp * count + a * Components + c) * Components + entry] !=
                    (entry == c ? phi : Scalar(0)))
                    return false;
                }
              }
            }
          }
          return true;
        };
        size_t retained = 0;
        if (!fresh)
        {
          for (size_t i = 0; i < identities; ++i)
          {
            if (!sweep.template operator()<true>(i) || !valid(i))
            {
              fail("Binding protocol differs from direct component bases");
              return;
            }
          }
          output.front() += Scalar(1);
          if (valid(identities - 1))
          {
            fail("Direct basis oracle accepted changed output");
            return;
          }
          for (const auto& formula : formulas)
            retained += bool(scalar.findTabulation(formula->getCacheIdentity()));
          if (churn && retained == identities)
          {
            fail("All churn identities are retained; review the cache-size protocol");
            return;
          }
          // Re-prime after inspection so public presence queries cannot determine
          // the measured initial cache order or the last bound shape sample.
          for (size_t i = 0; i < identities; ++i)
            sweep.template operator()<false>(i);
        }
        size_t next = 0, last = identities - 1;
        for (auto _ : state)
        {
          last = next;
          sweep.template operator()<false>(last);
          if (churn)
            next = (next + 1) % identities;
          benchmark::DoNotOptimize(output.data());
          benchmark::ClobberMemory();
        }
        if (!valid(last))
          fail("Timed binding protocol changed the qualified basis values");
        if (fresh)
        {
          if (!sweep.template operator()<true>(last) || !valid(last))
            fail("Fresh binding dimensions or reference authority are invalid");
          output.front() += Scalar(1);
          if (valid(last))
            fail("Fresh binding oracle accepted changed output");
        }
        size_t afterTables = 0;
        for (const auto& formula : formulas)
          afterTables += bool(scalar.findTabulation(formula->getCacheIdentity()));
        state.counters["dimension"] = mesh.getDimension();
        state.counters["ambient_dimension"] = mesh.getSpaceDimension();
        state.counters["degree"] = K;
        state.counters["field_rows"] = Rows;
        state.counters["field_columns"] = Columns;
        state.counters["components"] = Components;
        state.counters["local_basis"] = count;
        state.counters["scalar_basis"] = scalar.getCount();
        state.counters["quadrature_order"] = order;
        state.counters["quadrature_points"] = prototype.getSize();
        state.counters["bound_samples"] = samples;
        state.counters["formula_identities"] = identities;
        state.counters["primed_reference_tables"] = retained;
        state.counters["post_timed_reference_tables"] = afterTables;
        state.counters["first_binding"] = fresh;
        state.counters["reference_authority"] = associated;
        state.counters["stored_scalar_values"] = output.size();
        state.counters["retained_output_bytes"] = output.size() * sizeof(Scalar);
      }

      /// @brief Registers association/identity protocols at three sample orders.
      static void registerCases()
      {
        using Geometry::Polytope;
        const std::string range = FormLanguage::IsMatrixRange<Range>::Value ? "Matrix"
          : FormLanguage::IsVectorRange<Range>::Value                       ? "Vector"
                                                                            : "Scalar";
        for (auto geometry :
          {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          for (auto protocol : {H1BasisIdentityProtocol::AssociatedFresh,
                 H1BasisIdentityProtocol::AssociatedOne,
                 H1BasisIdentityProtocol::AssociatedNine,
                 H1BasisIdentityProtocol::UnprovenNine})
          {
            for (bool whole : {false, true})
            {
              const std::string name = std::string("H1/BasisLifecycle/") +
                (protocol == H1BasisIdentityProtocol::AssociatedFresh ? "AssociatedFresh/"
                    : protocol == H1BasisIdentityProtocol::AssociatedOne
                    ? "AssociatedOne/"
                    : protocol == H1BasisIdentityProtocol::AssociatedNine
                    ? "AssociatedNine/"
                    : "UnprovenNine/") +
                (whole ? "Sweep/" : "Single/") + "P" + std::to_string(K) + "/" +
                (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") + range + "/" +
                (geometry == Polytope::Type::Point
                    ? "Point"
                    : std::string(Convergence::UniformGrid::getGeometryName(geometry)));
              auto* entry = benchmark::RegisterBenchmark(name.c_str(),
                [geometry, protocol, whole](
                  auto& state) { run(state, geometry, protocol, whole); })
                              ->Arg(2)
                              ->Arg(4)
                              ->Arg(8)
                              ->UseRealTime();
              if (protocol == H1BasisIdentityProtocol::AssociatedFresh)
                entry->Iterations(1)->Repetitions(1);
            }
          }
        }
      }
  };

  /// @brief Registers one degree across scalar fields and supported range shapes.
  template <size_t K>
  void registerH1BasisLifecycleDegree()
  {
    H1BasisLifecycle<K, Real>::registerCases();
    H1BasisLifecycle<K, Complex>::registerCases();
    H1BasisLifecycle<K, Math::SpatialVector<Real>>::registerCases();
    H1BasisLifecycle<K, Math::SpatialVector<Complex>>::registerCases();
    H1BasisLifecycle<K, Math::SpatialMatrix<Real>>::registerCases();
    H1BasisLifecycle<K, Math::SpatialMatrix<Complex>>::registerCases();
  }
}

#endif
