/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Isolates scalar reference-element dispatch, evaluation and tabulation.
 *
 * No physical map, mesh, global assembly or solver is included. These local
 * kernels are shared by local, threaded and distributed consumers; this
 * executable does not measure backend scheduling or communication.
 */
#include <benchmark/benchmark.h>
#include <algorithm>
#include <array>
#include <functional>
#include <type_traits>
#include <vector>
#include "Rodin/Variational/P0/P0Element.h"
#include "Rodin/Variational/P0g/P0gElement.h"
#include "Rodin/Variational/P1/P1Element.h"
#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

namespace Rodin::Tests::Benchmarks
{
  bool referenceElementFailed = false;

  /// @brief Documented store capacity, checked by the resident-table guards.
  constexpr int ReferenceCacheCapacity = 8;

  enum class ReferenceOperation
  {
    Values,
    Derivatives,
    BulkDerivatives,
    PrefetchedValues,
    PrefetchedDerivatives,
    WarmTable,
    CapacityWarmMiss,
    WarmTablePosition
  };

  /**
   * @brief Measures one reference operation with explicit dispatch/cache scope.
   *
   * Direct and prefetched sweeps write the same quadrature-major output. The
   * former obtains each descriptor in the sample loop; the latter obtains
   * descriptors before timing. Neither path changes the library implementation.
   * BulkDerivatives uses the element's typed gradient output, including its
   * coordinate validation and owned point-local modal work. Its output has the
   * same quadrature-major layout, without physical mapping or field contraction.
   * WarmTable includes lookup only, not a sweep over table entries. Misses rotate
   * nine independently identified formulas through the eight-entry reference
   * cache. Initial allocation and basis initialization are outside timing;
   * CapacityWarmMiss therefore measures rebuilding, not first-use allocation.
   * WarmTable records its one-based search position outside timing: resident
   * table addresses order the cache's array slots, not formula identities.
   * WarmTablePosition selects each of the eight resident search positions
   * explicitly before timing. Formula contents are identical across versions;
   * lookup continues to use logical identity and geometry, never an address key.
   *
   * ## Numerical qualification
   *
   * Every value and derivative is compared exactly against direct evaluation
   * before and after timing. Separately, the public DOF functionals reconstruct
   * @f$f(\widehat x)=1+\sum_i(i+1)\widehat x_i@f$ and its gradient in
   * positive-degree spaces; constant spaces reconstruct @f$f=1@f$. This avoids
   * assuming nodal coefficients for modal elements. The reconstruction tolerance
   * is an absolute roundoff allowance on a dimensionless reference cell, not
   * a performance acceptance threshold. A perturbed output must be rejected.
   *
   * @tparam Element Scalar reference-element family.
   * @tparam Constant Whether only constants are represented.
   * @param state Quadrature order, optional cache position and timing iterations.
   * @param geometry Reference-cell family.
   * @param operation Measured kernel or cache protocol.
   */
  template <class Element, bool Constant>
  void referenceElement(benchmark::State& state, Geometry::Polytope::Type geometry,
    ReferenceOperation operation)
  {
    using Scalar = typename Element::ScalarType;
    using Store = Variational::ReferenceTabulationCache<Element, Scalar>;
    constexpr size_t FormulaCount = ReferenceCacheCapacity + 1;
    constexpr Real ReconstructionTolerance = 1e-8;
    Element element(geometry);
    const size_t ndof = element.getCount();
    const size_t dim = Geometry::Polytope::Traits(geometry).getDimension();
    std::vector<QF::PolytopeQuadratureFormula> formulas;
    formulas.reserve(FormulaCount);
    for (size_t i = 0; i < FormulaCount; ++i)
      formulas.emplace_back(state.range(0), geometry);
    const auto& formula = formulas.front();
    const size_t nqp = formula.getSize();
    using BasisResult = decltype(element.getBasis(size_t(0)));
    using Basis = std::remove_cvref_t<BasisResult>;
    using StoredBasis = std::conditional_t<std::is_lvalue_reference_v<BasisResult>,
      std::reference_wrapper<const Basis>, Basis>;
    std::vector<StoredBasis> bases;
    bases.reserve(ndof);
    for (size_t a = 0; a < ndof; ++a)
      bases.emplace_back(element.getBasis(a));
    const bool derivatives = operation == ReferenceOperation::Derivatives ||
      operation == ReferenceOperation::BulkDerivatives ||
      operation == ReferenceOperation::PrefetchedDerivatives;
    const bool prefetched = operation == ReferenceOperation::PrefetchedValues ||
      operation == ReferenceOperation::PrefetchedDerivatives;
    const bool warmLookup = operation == ReferenceOperation::WarmTable ||
      operation == ReferenceOperation::WarmTablePosition;
    const auto expected = Store::get(element, formula);
    std::vector<Scalar> output(derivatives ? nqp * ndof * dim : nqp * ndof);
    const auto sweep = [&]() {
      for (size_t qp = 0; qp < nqp; ++qp)
      {
        const auto& point = formula.getPoint(qp);
        if constexpr (requires {
                        element.tabulate(typename Element::Gradients{{}, 0}, point);
                      })
        {
          if (operation == ReferenceOperation::BulkDerivatives)
          {
            element.tabulate(
              typename Element::Gradients{
                std::span(output).subspan(qp * ndof * dim, ndof * dim), dim},
              point);
            continue;
          }
        }
        for (size_t a = 0; a < ndof; ++a)
        {
          const auto evaluate = [&](const auto& basis) {
            if (derivatives)
            {
              for (size_t i = 0; i < dim; ++i)
              {
                output[(qp * ndof + a) * dim + i] =
                  basis.template getDerivative<1>(i)(point);
              }
            }
            else
              output[qp * ndof + a] = basis(point);
          };
          if (!prefetched)
            evaluate(element.getBasis(a));
          else if constexpr (std::is_lvalue_reference_v<BasisResult>)
            evaluate(bases[a].get());
          else
            evaluate(bases[a]);
        }
      }
    };
    const auto fail = [&](const char* message) {
      referenceElementFailed = true;
      state.SkipWithError(message);
    };
    const Scalar phase = []() {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, Real(1) / 2);
      else
        return Scalar(1);
    }();
    const auto affine = [&](const Math::SpatialPoint& point) -> Scalar {
      Scalar value = 1;
      if constexpr (!Constant)
        for (size_t i = 0; i < dim; ++i)
          value += Scalar(i + 1) * point(i);
      return phase * value;
    };
    std::vector<Scalar> coefficients(ndof);
    for (size_t a = 0; a < ndof; ++a)
      coefficients[a] = element.getLinearForm(a)(affine);
    for (size_t qp = 0; qp < nqp; ++qp)
    {
      Scalar value = 0;
      std::array<Scalar, 3> gradient{};
      for (size_t a = 0; a < ndof; ++a)
      {
        value += coefficients[a] * expected.getBasis(qp, a);
        for (size_t i = 0; i < dim; ++i)
          gradient[i] += coefficients[a] * expected.template getDerivative<1>(qp, a, i);
      }
      const Real error = std::abs(value - affine(formula.getPoint(qp)));
      if (!std::isfinite(error) || error > ReconstructionTolerance)
      {
        fail("Reference-field reconstruction failed");
        return;
      }
      for (size_t i = 0; i < dim; ++i)
      {
        const Scalar exact = Constant ? Scalar(0) : phase * Scalar(i + 1);
        const Real derivativeError = std::abs(gradient[i] - exact);
        if (!std::isfinite(derivativeError) || derivativeError > ReconstructionTolerance)
        {
          fail("Reference-gradient reconstruction failed");
          return;
        }
      }
    }
    const auto& reference =
      derivatives ? expected.getDerivatives() : expected.getValues();
    sweep();
    if (output != reference)
    {
      fail("Direct/prefetched entries differ from reference tabulation");
      return;
    }
    if (!output.empty())
    {
      output.front() += Scalar(1);
      if (output == reference)
      {
        fail("Perturbed-output sensitivity check failed");
        return;
      }
      sweep();
    }
    for (const auto& qf : formulas)
      Store::get(element, qf);
    const QF::QuadratureFormulaBase* lookupFormula = &formula;
    if (operation == ReferenceOperation::WarmTablePosition)
    {
      using Resident =
        std::pair<const typename Store::Table*, const QF::QuadratureFormulaBase*>;
      std::vector<Resident> residents;
      residents.reserve(FormulaCount - 1);
      for (const auto& qf : formulas)
      {
        if (const auto table = Store::find(element, qf.getCacheIdentity()))
          residents.emplace_back(&table->get(), &qf);
      }
      const size_t requested = state.range(1);
      if (residents.size() != FormulaCount - 1 || requested == 0 ||
        requested > residents.size())
      {
        fail("Pinned lookup requires eight resident tables and a valid position");
        return;
      }
      std::sort(residents.begin(), residents.end(), [](const auto& lhs, const auto& rhs) {
        return std::less<const typename Store::Table*>{}(lhs.first, rhs.first);
      });
      lookupFormula = residents[requested - 1].second;
    }
    if (warmLookup)
    {
      const auto* selected = &Store::get(element, *lookupFormula);
      size_t resident = 0, position = 1;
      for (const auto& qf : formulas)
      {
        if (const auto table = Store::find(element, qf.getCacheIdentity()))
        {
          ++resident;
          position += std::less<const typename Store::Table*>{}(&table->get(), selected);
        }
      }
      if (resident != FormulaCount - 1)
      {
        fail("Warm lookup instrumentation did not observe all eight resident tables");
        return;
      }
      // Entries are stored in one array and searched in increasing slot order.
      // Address order diagnoses that search position; it is never a cache key.
      state.counters["warm_lookup_position"] = position;
      if (operation == ReferenceOperation::WarmTablePosition &&
        position != static_cast<size_t>(state.range(1)))
      {
        fail("Pinned lookup did not select the requested search position");
        return;
      }
    }
    if (operation == ReferenceOperation::CapacityWarmMiss)
    {
      for (const auto& qf : formulas)
      {
        if (Store::find(element, qf.getCacheIdentity()))
        {
          fail("Cache-miss protocol unexpectedly found a resident table");
          return;
        }
        const auto& table = Store::get(element, qf);
        if (table.getValues() != expected.getValues() ||
          table.getDerivatives() != expected.getDerivatives())
        {
          fail("Cache-miss warmup changed reference entries");
          return;
        }
      }
    }
    size_t next = 0;
    for (auto _ : state)
    {
      if (warmLookup)
      {
        const auto& table = Store::get(element, *lookupFormula);
        benchmark::DoNotOptimize(table.getValues().data());
      }
      else if (operation == ReferenceOperation::CapacityWarmMiss)
      {
        const auto& table = Store::get(element, formulas[next]);
        benchmark::DoNotOptimize(table.getValues().data());
        next = (next + 1) % FormulaCount;
      }
      else
      {
        sweep();
        benchmark::DoNotOptimize(output.data());
      }
      benchmark::ClobberMemory();
    }
    if (output != reference)
      fail("Timed sweep changed reference entries");
    const auto& table = Store::get(element, *lookupFormula);
    if (table.getValues() != expected.getValues() ||
      table.getDerivatives() != expected.getDerivatives())
      fail("Timed tabulation changed reference entries");
    state.counters["dimension"] = dim;
    state.counters["basis_functions"] = ndof;
    state.counters["quadrature_points"] = nqp;
    state.counters["table_bytes"] = sizeof(Scalar) * nqp * ndof * (1 + dim);
    state.counters["output_entries"] = output.size();
    state.counters["written_entries"] = warmLookup
      ? 0
      : (operation == ReferenceOperation::CapacityWarmMiss ? nqp * ndof * (1 + dim)
                                                           : output.size());
    state.counters["formula_identities"] =
      operation == ReferenceOperation::CapacityWarmMiss ? FormulaCount : 1;
  }

  /// @brief Registers reference operations and explicit cache search positions.
  template <class Element, bool Constant = false>
  void registerReferenceElement(const char* family, const char* scalar)
  {
    using G = Geometry::Polytope::Type;
    for (const auto& [geometry, label] :
      {std::pair{G::Point, "Point"}, {G::Segment, "Segment"}, {G::Triangle, "Triangle"},
        {G::Quadrilateral, "Quadrilateral"}, {G::Tetrahedron, "Tetrahedron"},
        {G::Hexahedron, "Hexahedron"}, {G::Pyramid, "Pyramid"}, {G::Wedge, "Wedge"}})
    {
      for (const auto& [operation, operationLabel] :
        {std::pair{ReferenceOperation::Values, "Values"},
          {ReferenceOperation::Derivatives, "Derivatives"},
          {ReferenceOperation::BulkDerivatives, "BulkDerivatives"},
          {ReferenceOperation::PrefetchedValues, "PrefetchedValues"},
          {ReferenceOperation::PrefetchedDerivatives, "PrefetchedDerivatives"},
          {ReferenceOperation::WarmTable, "WarmTable"},
          {ReferenceOperation::WarmTablePosition, "WarmTablePosition"},
          {ReferenceOperation::CapacityWarmMiss, "CapacityWarmMiss"}})
      {
        if (geometry == G::Point &&
          (operation == ReferenceOperation::Derivatives ||
            operation == ReferenceOperation::BulkDerivatives ||
            operation == ReferenceOperation::PrefetchedDerivatives))
          continue;
        if constexpr (!requires(const Element& element, const Math::SpatialPoint& point) {
                        element.tabulate(typename Element::Gradients{{}, 0}, point);
                      })
        {
          if (operation == ReferenceOperation::BulkDerivatives)
            continue;
        }
        const std::string name = std::string("ReferenceElement/") + family + "/" +
          scalar + "/" + operationLabel + "/" + label;
        auto* registration =
          benchmark::RegisterBenchmark(name.c_str(), [geometry, operation](auto& state) {
            referenceElement<Element, Constant>(state, geometry, operation);
          });
        for (const int order : {2, 4, 6, 8, 12})
        {
          if (operation == ReferenceOperation::WarmTablePosition)
          {
            for (int position = 1; position <= ReferenceCacheCapacity; ++position)
              registration->Args({order, position});
          }
          else
            registration->Arg(order);
        }
        registration->UseRealTime();
      }
    }
  }

  /// @brief Registers both scalar fields for one supported H1 degree.
  template <size_t K>
  void registerReferenceH1()
  {
    const std::string family = "H1/P" + std::to_string(K);
    registerReferenceElement<Variational::H1Element<K, Real>>(family.c_str(), "Real");
    registerReferenceElement<Variational::H1Element<K, Complex>>(
      family.c_str(), "Complex");
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Variational;
  using namespace Rodin::Tests::Benchmarks;
  registerReferenceElement<P0Element<Real>, true>("P0", "Real");
  registerReferenceElement<P0Element<Complex>, true>("P0", "Complex");
  registerReferenceElement<P0gElement<Real>, true>("P0g", "Real");
  registerReferenceElement<P0gElement<Complex>, true>("P0g", "Complex");
  registerReferenceElement<P1Element<Real>>("P1", "Real");
  registerReferenceElement<P1Element<Complex>>("P1", "Complex");
  registerReferenceH1<1>();
  registerReferenceH1<2>();
  registerReferenceH1<3>();
  registerReferenceH1<4>();
  registerReferenceH1<5>();
  registerReferenceH1<6>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || referenceElementFailed ? 1 : 0;
}
