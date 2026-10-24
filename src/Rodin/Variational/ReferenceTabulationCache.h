/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_REFERENCETABULATIONCACHE_H
#define RODIN_VARIATIONAL_REFERENCETABULATIONCACHE_H

/**
 * @file
 * @brief Shared reference-basis tabulation for scalar finite elements.
 */
#include <array>
#include <cassert>
#include <functional>
#include <span>
#include <type_traits>
#include <vector>
#include "Rodin/Types.h"
#include "ReferenceTabulation.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Variational
{
  /**
   * @brief Reuses scalar basis values and first derivatives at reference samples.
   *
   * For the element's own basis @f$\phi_a@f$ and reference samples
   * @f$\widehat x_q@f$, the table stores @f$\phi_a(\widehat x_q)@f$ and
   * @f$\partial_i\phi_a(\widehat x_q)@f$ in quadrature-major order.
   * No basis conversion, physical coordinates or mapped derivatives are stored.
   *
   * ## Architecture
   *
   * The bounded thread-local store is distinct for each element/scalar type.
   * Persistent formula address and reference geometry identify an entry.
   * Published formulas own immutable samples until process termination.
   * Borrowed tables are valid only until a subsequent insertion evicts their
   * slot; callers must consume them immediately.
   * Arbitrary coefficient callbacks can insert tables even when their returned
   * mathematical value is constant. Evaluate such callbacks before borrowing
   * each table, and finish consuming it before invoking another callback.
   * Eviction and cross-thread misses permit direct pointwise evaluation.
   * Cells share reference data, not geometry-dependent matrices.
   * Both public lookup operations use the same formula search. Construction
   * is a separate stage: it retires an evicted key before changing storage and
   * publishes the new key only after the element has filled the complete table.
   * A request obtains its thread-local cache once and shares that reference
   * between its lookup and construction stages.
   * An element's basis must be determined solely by its type and reference
   * geometry, as for P1 and fixed-degree H1; mutable basis parameters are not
   * represented by this key.
   * An element may supply a public bulk tabulation operation. It sizes and fills the
   * caller-owned table; cache ownership, lookup, eviction and publication stay
   * here. Elements without that operation use direct basis evaluation.
   *
   * @tparam Element Scalar finite element defining basis order and evaluation.
   * @tparam Scalar Scalar stored in the reference table.
   */
  template <class Element, class Scalar>
  class ReferenceTabulationCache final
  {
    public:
      /// @brief Owned numerical table, independent of cache policy.
      using Table = ReferenceTabulation<Scalar>;

    private:
      struct Entry
      {
        Optional<const QF::QuadratureFormula*> formula;
        Geometry::Polytope::Type geometry;
        Table table;

        Entry()
          : formula(),
            geometry(Geometry::Polytope::Type::Point),
            table()
        {}
      };

      struct Cache
      {
        /// @brief Retains the existing H1 cache's eight-entry capacity.
        static constexpr size_t Capacity = 8;
        std::array<Entry, Capacity> entries;
        size_t next;

        Cache()
          : entries(),
            next(0)
        {}
      };

      /**
       * @brief Gets this element/scalar type's thread-local store.
       * @returns The bounded cache owned by this thread.
       */
      static Cache& getCache()
      {
        static thread_local Cache cache;
        return cache;
      }

    public:
      /**
       * @brief Borrows an existing table without evaluating a formula.
       * @param element Element specifying the reference geometry.
       * @param formula Persistent formula defining the samples.
       * @returns A matching table in this thread, or no value.
       */
      static Optional<std::reference_wrapper<const Table>> find(
        const Element& element, const QF::QuadratureFormula& formula)
      {
        return find(getCache(), element, formula);
      }

      /**
       * @brief Gets or constructs the reference table for a formula's contents.
       *
       * A matching formula address and reference geometry identify immutable
       * samples. Hits do not query the formula's size or coordinates. Samples
       * are read only while constructing a missing table.
       *
       * @param element Element providing scalar basis values and derivatives.
       * @param formula Formula defining the reference samples.
       * @returns A borrowed table, consumed before another insertion can evict it.
       */
      static const Table& get(
        const Element& element, const QF::QuadratureFormula& formula)
      {
        auto& cache = getCache();
        if (const auto existing = find(cache, element, formula))
          return existing->get();
        return insert(cache, element, formula);
      }

    private:
      /**
       * @brief Searches one already obtained cache without evaluating samples.
       * @param cache This thread's existing reference-table store.
       * @param element Element specifying the reference geometry.
       * @param formula Persistent formula defining the samples.
       * @returns A matching table or no value, without changing the store.
       */
      static Optional<std::reference_wrapper<const Table>> find(
        const Cache& cache, const Element& element, const QF::QuadratureFormula& formula)
      {
        for (const auto& entry : cache.entries)
        {
          if (entry.formula == std::addressof(formula) &&
            entry.geometry == element.getGeometry())
            return std::cref(entry.table);
        }
        return {};
      }

      /**
       * @brief Constructs and publishes a table after a formula lookup misses.
       *
       * Eviction retires the
       * previous formula association before resizing or filling, so an interrupted fill
       * cannot publish partial data under either key. The round-robin policy
       * and element-owned arithmetic are unchanged.
       *
       * @param cache This thread's store, already obtained during lookup.
       * @param element Element providing the reference geometry and basis.
       * @param formula Immutable samples whose table is absent from this store.
       * @returns A newly filled table borrowed until its slot is evicted.
       */
      static const Table& insert(
        Cache& cache, const Element& element, const QF::QuadratureFormula& formula)
      {
        const auto* pointer = std::addressof(formula);
        const auto geometry = element.getGeometry();
        auto& entry = cache.entries[cache.next];
        cache.next = (cache.next + 1) % cache.entries.size();
        entry.formula.reset();
        auto& table = entry.table;
        if constexpr (requires { element.tabulate(table, formula); })
          element.tabulate(table, formula);
        else
        {
          table.resize(formula.getSize(), element.getCount(),
            Geometry::Polytope::Traits(geometry).getDimension());
          fill(element, table, formula);
        }
        entry.geometry = geometry;
        entry.formula = pointer;
        return table;
      }

      /**
       * @brief Fills an already-sized table by direct basis evaluation.
       *
       * Elements without a public bulk operation use this cache-local fallback.
       * It consumes only public basis operations and changes no cache key.
       *
       * @param element Element providing values and directional derivatives.
       * @param table Owned table with dimensions and storage already established.
       * @param formula Reference samples defining the table rows.
       */
      static void fill(
        const Element& element, Table& table, const QF::QuadratureFormulaBase& formula)
      {
        // Retain immutable descriptors by reference and own returned values.
        // Acquire either representation once per basis.
        using BasisResult = decltype(element.getBasis(size_t(0)));
        using Basis = std::remove_cvref_t<BasisResult>;
        using StoredBasis = std::conditional_t<std::is_lvalue_reference_v<BasisResult>,
          std::reference_wrapper<const Basis>, Basis>;
        std::vector<StoredBasis> bases;
        bases.reserve(table.getCount());
        for (size_t a = 0; a < table.getCount(); ++a)
          bases.emplace_back(element.getBasis(a));
        for (size_t qp = 0; qp < table.getSize(); ++qp)
        {
          const auto& reference = formula.getPoint(qp);
          for (size_t a = 0; a < table.getCount(); ++a)
          {
            const Basis& basis = [&]() -> const Basis& {
              if constexpr (std::is_lvalue_reference_v<BasisResult>)
                return bases[a].get();
              else
                return bases[a];
            }();
            table.getValues(qp).data[a] = basis(reference);
            for (size_t i = 0; i < table.getDimension(); ++i)
            {
              table.getGradients(qp).data[a * table.getDimension() + i] =
                basis.template getDerivative<1>(i)(reference);
            }
          }
        }
      }
  };
}
#endif
