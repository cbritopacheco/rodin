/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_QUADRATUREFORMULA_H
#define RODIN_QF_QUADRATUREFORMULA_H

/**
 * @file
 * @brief Reference quadrature access and persistent immutable formulas.
 */

#include <array>
#include <cassert>
#include <functional>
#include <memory>
#include <mutex>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "Rodin/Math.h"
#include "Rodin/Types.h"
#include "Rodin/Geometry/Polytope.h"

#ifndef RODIN_QF_QUADRATUREFORMULA_LOOKUP_CAPACITY
/// @brief Positive per-thread lookup budget, configured consistently across translation units.
#define RODIN_QF_QUADRATUREFORMULA_LOOKUP_CAPACITY 8
#endif

namespace Rodin::QF
{
  /**
   * @defgroup RodinQuadrature Quadrature formulae
   * @brief Quadrature formulas on reference polytopes.
   */

  /**
   * @brief Read-only access to reference quadrature points and weights.
   * @ingroup RodinQuadrature
   *
   * Concrete families own their samples. Immediate numerical operations may
   * borrow this interface while its owner remains alive and unchanged.
   * Persistent consumers require @ref QuadratureFormula instead.
   */
  class QuadratureFormulaBase
  {
    public:
      /// @brief Destroys the formula through its value-access interface.
      virtual ~QuadratureFormulaBase() = default;

      /// @returns Number of reference samples.
      virtual size_t getSize() const = 0;

      /**
       * @param index Reference sample index.
       * @returns Weight of the reference sample.
       */
      virtual Real getWeight(size_t index) const = 0;

      /**
       * @param index Reference sample index.
       * @returns Coordinates borrowed until the formula is changed or destroyed.
       */
      virtual const Math::SpatialVector<Real>& getPoint(size_t index) const = 0;
  };

  /**
   * @brief Persistent immutable reference quadrature points and weights.
   * @ingroup RodinQuadrature
   *
   * ## Architecture
   *
   * Each formula family and constructor-key type has a shared owning pool.
   * Integer and enum arguments are copied into a tuple; published formula
   * arguments are borrowed through their permanent addresses. Tuple comparison
   * uses total pointer ordering for those references. No caller-owned formula
   * can enter a persistent cache through this interface.
   *
   * A bounded thread-local lookup borrows published formulas. Shared lookups and
   * insertion hold the pool mutex. Construction occurs outside the mutex so
   * generators may request other formulas. Concurrent misses may construct
   * duplicate candidates, but publish only one. Pools retain their entries until
   * process termination, including throughout static destruction. All published
   * samples are owned and immutable; reading them requires no lock.
   *
   * Constructor arguments must determine all generated samples. Generators must
   * not depend on unsynchronized mutable state or cyclic factory requests.
   */
  class QuadratureFormula final : public QuadratureFormulaBase
  {
    /// @brief Borrowed persistent formula in a constructor-argument tuple.
    struct Reference
    {
      /// @brief Address of the borrowed persistent formula.
      const QuadratureFormula* pointer;

      /**
       * @param other Formula reference compared with this reference.
       * @returns Whether both references designate the same formula.
       */
      bool operator==(const Reference& other) const = default;

      /**
       * @param other Formula reference compared with this reference.
       * @returns Whether this address precedes the other in the total ordering.
       */
      bool operator<(const Reference& other) const noexcept
      {
        return std::less<const QuadratureFormula*>{}(pointer, other.pointer);
      }
    };

    /**
     * @brief Shared owning storage protected by its mutex.
     * @tparam Key Owned constructor-argument tuple.
     */
    template <class Key>
    struct Pool
    {
      /// @brief Protects all accesses to the owning map.
      std::mutex mutex;

      /// @brief Published formulas, never erased or replaced.
      FlatMap<Key, std::unique_ptr<const QuadratureFormula>> formulas;
    };

    /**
     * @brief Per-thread borrowed lookup with round-robin replacement.
     * @tparam Key Owned constructor-argument tuple.
     */
    template <class Key>
    struct Lookup
    {
      /// @brief One constructor key and its borrowed published formula.
      struct Entry
      {
        /// @brief Owned constructor-argument values.
        Key key;

        /// @brief Persistent formula selected by the key.
        const QuadratureFormula* pointer;
      };

      static_assert(RODIN_QF_QUADRATUREFORMULA_LOOKUP_CAPACITY > 0);

      /// @brief Bounded slots owned exclusively by this thread.
      std::array<Optional<Entry>, RODIN_QF_QUADRATUREFORMULA_LOOKUP_CAPACITY>
        entries{};

      /// @brief Slot selected for the next round-robin replacement.
      size_t next = 0;
    };

    public:
      /**
       * @brief Published formulas cannot be copied.
       * @param other Formula whose copying is disabled.
       */
      QuadratureFormula(const QuadratureFormula& other) = delete;

      /**
       * @brief Published formulas cannot be moved.
       * @param other Formula whose moving is disabled.
       */
      QuadratureFormula(QuadratureFormula&& other) = delete;

      /**
       * @brief Published formulas cannot be reassigned.
       * @param other Formula whose assignment is disabled.
       * @returns Unavailable because assignment is disabled.
       */
      QuadratureFormula& operator=(const QuadratureFormula& other) = delete;

      /**
       * @brief Published formulas cannot be move-assigned.
       * @param other Formula whose assignment is disabled.
       * @returns Unavailable because assignment is disabled.
       */
      QuadratureFormula& operator=(QuadratureFormula&& other) = delete;

      /// @brief Destroys the owned reference samples.
      ~QuadratureFormula() override = default;

      /// @returns Number of immutable reference samples.
      size_t getSize() const override
      {
        return m_points.size();
      }

      /**
       * @param index Reference sample index.
       * @returns Immutable reference weight.
       */
      Real getWeight(size_t index) const override
      {
        assert(index < m_weights.size());
        return m_weights[index];
      }

      /**
       * @param index Reference sample index.
       * @returns Reference coordinates retained until process termination.
       */
      const Math::SpatialVector<Real>& getPoint(size_t index) const override
      {
        assert(index < m_points.size());
        return m_points[index];
      }

      /**
       * @brief Constructs or retrieves a persistent formula of a concrete family.
       * @tparam Formula Concrete generator implementing QuadratureFormulaBase.
       * @tparam Args Integer, enum, or persistent-formula constructor arguments.
       * @param args Arguments forwarded to the generator on a cache miss.
       * @returns Immutable formula retained until process termination.
       *
       * Argument value types select distinct pools; family factories normalize
       * counts and orders to size_t. Generator exceptions publish no entry.
       * Returned references are independent of constructor-argument lifetimes.
       */
      template <class Formula, class... Args>
#if __has_cpp_attribute(gnu::no_dangling)
      [[gnu::no_dangling]]
#endif
      static const QuadratureFormula& get(Args&&... args)
      {
        static_assert(std::is_base_of_v<QuadratureFormulaBase, Formula>);
        static_assert(!std::is_same_v<Formula, QuadratureFormula>);
        static_assert(std::is_constructible_v<Formula, Args&&...>);
        auto key = std::tuple{[]<class Argument>(const Argument& argument) {
          using Value = std::remove_cvref_t<Argument>;
          if constexpr (std::is_same_v<Value, QuadratureFormula>)
          {
            return Reference{std::addressof(argument)};
          }
          else
          {
            static_assert(std::is_integral_v<Value> || std::is_enum_v<Value>,
              "Unsupported quadrature factory argument type.");
            return Value(argument);
          }
        }(args)...};
        using Key = decltype(key);
        auto& lookup = getLookup<Formula, Key>();
        for (const auto& entry : lookup.entries)
        {
          if (entry && entry->key == key)
            return *entry->pointer;
        }
        auto& pool = getPool<Formula, Key>();
        const QuadratureFormula* pointer = nullptr;
        {
          const std::lock_guard lock(pool.mutex);
          const auto found = pool.formulas.find(key);
          if (found != pool.formulas.end())
            pointer = found->second.get();
        }
        if (!pointer)
        {
          Formula generated(std::forward<Args>(args)...);
          auto candidate =
            std::unique_ptr<const QuadratureFormula>(new QuadratureFormula(generated));
          {
            const std::lock_guard lock(pool.mutex);
            const auto found = pool.formulas.find(key);
            if (found != pool.formulas.end())
              pointer = found->second.get();
            else
            {
              const auto inserted = pool.formulas.emplace(key, std::move(candidate));
              pointer = inserted.first->second.get();
            }
          }
        }
        using Entry = typename Lookup<Key>::Entry;
        lookup.entries[lookup.next] = Entry{key, pointer};
        lookup.next = (lookup.next + 1) % lookup.entries.size();
        return *pointer;
      }

    private:
      /**
       * @brief Owns all generated samples before publishing the formula.
       * @param formula Generator borrowed only during construction.
       */
      explicit QuadratureFormula(const QuadratureFormulaBase& formula)
      {
        const size_t count = formula.getSize();
        m_points.reserve(count);
        m_weights.reserve(count);
        for (size_t index = 0; index < count; ++index)
        {
          m_points.push_back(formula.getPoint(index));
          m_weights.push_back(formula.getWeight(index));
        }
      }

      /**
       * @tparam Formula Concrete family selecting its own pool.
       * @tparam Key Owned constructor-argument tuple.
       * @returns Shared pool retained throughout static destruction.
       */
      template <class Formula, class Key>
      static Pool<Key>& getPool()
      {
        static auto* const pool = new Pool<Key>;
        return *pool;
      }

      /**
       * @tparam Formula Concrete family selecting its own lookup.
       * @tparam Key Owned constructor-argument tuple.
       * @returns This thread's bounded borrowed lookup.
       */
      template <class Formula, class Key>
      static Lookup<Key>& getLookup()
      {
        static thread_local Lookup<Key> lookup;
        return lookup;
      }

      /// @brief Owned reference coordinates in sample order.
      std::vector<Math::SpatialVector<Real>> m_points;

      /// @brief Owned reference weights in sample order.
      std::vector<Real> m_weights;
  };
}

#endif
