/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_REFERENCETABULATION_H
#define RODIN_VARIATIONAL_REFERENCETABULATION_H

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
   * Logical formula identity and reference geometry identify an entry, never
   * approximate coordinate equality. Assignment or reconstruction of a formula
   * changes its identity. Lookup by identity neither dereferences the formula
   * nor creates an entry. Borrowed tables are valid only until a subsequent
   * insertion evicts their slot; callers must consume them immediately.
   * Eviction and cross-thread misses permit direct pointwise evaluation.
   * Cells share reference data, not geometry-dependent matrices.
   * An element's basis must be determined solely by its type and reference
   * geometry, as for P1 and fixed-degree H1; mutable basis parameters are not
   * represented by this key.
   *
   * @tparam Element Scalar finite element defining basis order and evaluation.
   * @tparam Scalar Scalar stored in the reference table.
   */
  template <class Element, class Scalar>
  class ReferenceTabulation final
  {
    public:
    /// @brief Reference samples and derivatives in the element's basis order.
    struct Table
    {
      size_t nqp; ///< Number of quadrature points.
      size_t ndof; ///< Number of scalar basis functions.
      size_t dim; ///< Reference dimension, including zero.
      std::vector<Scalar> phi; ///< Values indexed by qp*ndof+a.
      std::vector<Scalar> dphi; ///< Derivatives indexed by (qp*ndof+a)*dim+i.

      /// @brief Constructs an empty reference table.
      Table()
        : nqp(0),
          ndof(0),
          dim(0),
          phi(),
          dphi()
      {}

      /**
       * @brief Gets a reference basis value.
       * @param qp Quadrature-point index.
       * @param a Local basis-function index.
       * @returns The tabulated basis value.
       */
      const Scalar& getBasis(size_t qp, size_t a) const noexcept
      {
        return phi[qp * ndof + a];
      }

      /**
       * @brief Gets a tabulated reference derivative.
       * @tparam Order Zero selects values; one selects first derivatives.
       * Higher orders retain the H1 table's historical zero result.
       * @param qp Quadrature-point index.
       * @param a Local basis-function index.
       * @param i Reference coordinate, zero when selecting a value.
       * @returns The selected value, derivative, or historical higher-order zero.
       */
      template <size_t Order>
      const Scalar& getDerivative(size_t qp, size_t a, size_t i) const noexcept
      {
        if constexpr (Order == 0)
        {
          assert(i == 0);
          return getBasis(qp, a);
        }
        else if constexpr (Order == 1)
        {
          assert(i < dim);
          return dphi[(qp * ndof + a) * dim + i];
        }
        else
        {
          static const Scalar zero = Scalar(0);
          return zero;
        }
      }

      /**
       * @brief Borrows the reference gradient of one basis function.
       * @param qp Quadrature-point index.
       * @param a Local basis-function index.
       * @returns A span of reference dimension, empty on a point.
       */
      std::span<const Scalar> getGradient(size_t qp, size_t a) const noexcept
      {
        if (dim == 0)
          return {};
        return {dphi.data() + (qp * ndof + a) * dim, dim};
      }
    };

    private:
    struct Entry
    {
      Optional<size_t> identity;
      Geometry::Polytope::Type geometry;
      Table table;

      Entry()
        : identity(),
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
     * @param identity Logical formula lifetime/assignment identity.
     * @returns A matching table in this thread, or no value.
     */
    static Optional<std::reference_wrapper<const Table>> find(
      const Element& element, size_t identity)
    {
      for (const auto& entry : getCache().entries)
      {
        if (entry.identity == identity && entry.geometry == element.getGeometry())
          return std::cref(entry.table);
      }
      return {};
    }

    /**
     * @brief Gets or constructs the reference table for a formula's contents.
     * @param element Element providing scalar basis values and derivatives.
     * @param formula Formula defining the reference samples.
     * @returns A borrowed table, consumed before another insertion can evict it.
     */
    static const Table& get(
      const Element& element, const QF::QuadratureFormulaBase& formula)
    {
      auto& cache = getCache();
      const auto identity = formula.getCacheIdentity();
      const auto geometry = element.getGeometry();
      const size_t nqp = formula.getSize();
      for (const auto& entry : cache.entries)
      {
        if (entry.identity == identity && entry.geometry == geometry &&
          entry.table.nqp == nqp)
          return entry.table;
      }
      auto& entry = cache.entries[cache.next];
      cache.next = (cache.next + 1) % cache.entries.size();
      entry.identity.reset();
      auto& table = entry.table;
      table.nqp = nqp;
      table.ndof = element.getCount();
      table.dim = Geometry::Polytope::Traits(geometry).getDimension();
      table.phi.resize(nqp * table.ndof);
      table.dphi.resize(nqp * table.ndof * table.dim);
      // H1 returns immutable descriptors by reference; P1 returns values.
      // Acquire either representation once per basis, as in the H1 store.
      using BasisResult = decltype(element.getBasis(size_t(0)));
      using Basis = std::remove_cvref_t<BasisResult>;
      using StoredBasis = std::conditional_t<std::is_lvalue_reference_v<BasisResult>,
        std::reference_wrapper<const Basis>, Basis>;
      std::vector<StoredBasis> bases;
      bases.reserve(table.ndof);
      for (size_t a = 0; a < table.ndof; ++a)
        bases.emplace_back(element.getBasis(a));
      for (size_t qp = 0; qp < nqp; ++qp)
      {
        const auto& reference = formula.getPoint(qp);
        for (size_t a = 0; a < table.ndof; ++a)
        {
          const Basis& basis = [&]() -> const Basis& {
            if constexpr (std::is_lvalue_reference_v<BasisResult>)
              return bases[a].get();
            else
              return bases[a];
          }();
          table.phi[qp * table.ndof + a] = basis(reference);
          for (size_t i = 0; i < table.dim; ++i)
          {
            table.dphi[(qp * table.ndof + a) * table.dim + i] =
              basis.template getDerivative<1>(i)(reference);
          }
        }
      }
      entry.geometry = geometry;
      entry.identity = identity;
      return table;
    }
  };
}
#endif
