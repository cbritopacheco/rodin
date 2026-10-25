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
#include <span>
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
      /// @brief Read-only sample arrays borrowed from a concrete formula.
      struct Samples
      {
        /// @brief Contiguous reference coordinates in sample order.
        std::span<const Math::SpatialVector<Real>> points;

        /// @brief First weight, retained by the formula or immutable static storage.
        const Real* weights;

        /// @brief Scalar spacing between weights; zero repeats one uniform weight.
        size_t weightStride;
      };

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

      /**
       * @brief Borrows sample arrays when the concrete representation permits it.
       * @returns Array views valid until this formula is changed or destroyed,
       * or no value for formulas providing only pointwise access.
       */
      virtual Optional<Samples> getSamples() const
      {
        return {};
      }
  };

  /**
   * @brief Persistent immutable reference quadrature points and weights.
   * @ingroup RodinQuadrature
   *
   * # Architecture
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
   * concrete rules own their immutable samples; reading them requires no lock.
   *
   * Constructor arguments must determine all generated samples. Generators must
   * own their samples or borrow storage with process lifetime, and must not retain
   * caller-owned constructor arguments. They must not depend on unsynchronized
   * mutable state or cyclic factory requests.
   */
  class QuadratureFormula final : public QuadratureFormulaBase
  {
    /**
     * @brief One allocation owning the concrete rule and its published interface.
     * @tparam Formula Concrete formula constructed in place.
     */
    template <class Formula>
    struct Storage;

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
     * @tparam Formula Concrete formula retained by each entry.
     * @tparam Key Owned constructor-argument tuple.
     */
    template <class Formula, class Key>
    struct Pool
    {
      /// @brief Protects all accesses to the owning map.
      std::mutex mutex;

      /// @brief Published formulas, never erased or replaced.
      FlatMap<Key, std::unique_ptr<const Storage<Formula>>> formulas;
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

      /// @brief Destroys the published interface; its pool owns the concrete rule.
      ~QuadratureFormula() override = default;

      /// @returns Number of immutable reference samples.
      size_t getSize() const override
      {
        return m_size;
      }

      /**
       * @param index Reference sample index.
       * @returns Immutable reference weight.
       */
      Real getWeight(size_t index) const override
      {
        assert(index < m_size);
        if (m_samples)
          return m_samples->weights[index * m_samples->weightStride];
        return m_formula.get().getWeight(index);
      }

      /**
       * @param index Reference sample index.
       * @returns Reference coordinates retained until process termination.
       */
      const Math::SpatialVector<Real>& getPoint(size_t index) const override
      {
        assert(index < m_size);
        if (m_samples)
          return m_samples->points[index];
        return m_formula.get().getPoint(index);
      }

      /// @returns The retained rule's sample views, when available.
      Optional<Samples> getSamples() const override
      {
        return m_samples;
      }

      /**
       * @brief Selects a persistent quadrature formula for a reference geometry.
       * @param order Required polynomial degree of exactness.
       * @param geometry Reference polytope geometry.
       * @returns Immutable formula retained until process termination.
       *
       * Uses positive interior published tables where available and Gaussian
       * product rules beyond their supported degrees. Tensor factors are
       * persistent formulas from the same thread-safe factory.
       */
      static const QuadratureFormula& get(
        size_t order, Geometry::Polytope::Type geometry);

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
            pointer = &found->second->formula;
        }
        if (!pointer)
        {
          auto candidate =
            std::make_unique<const Storage<Formula>>(std::forward<Args>(args)...);
          {
            const std::lock_guard lock(pool.mutex);
            const auto found = pool.formulas.find(key);
            if (found != pool.formulas.end())
              pointer = &found->second->formula;
            else
            {
              const auto inserted = pool.formulas.emplace(key, std::move(candidate));
              pointer = &inserted.first->second->formula;
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
       * @brief Resolves a default request through the concrete-family factory.
       * @param order Required polynomial degree of exactness.
       * @param geometry Reference geometry.
       * @returns Permanent concrete-family formula, including persistent factors.
       */
#if __has_cpp_attribute(gnu::noinline)
      [[gnu::noinline]]
#endif
      static const QuadratureFormula& select(
        size_t order, Geometry::Polytope::Type geometry);

      /**
       * @brief Publishes the concrete rule owned by the same pool entry.
       * @param formula Immutable concrete rule retained with this interface.
       */
      explicit QuadratureFormula(const QuadratureFormulaBase& formula)
        : m_formula(formula),
          m_size(formula.getSize()),
          m_samples(formula.getSamples())
      {
        assert(!m_samples || m_samples->points.size() == m_size);
        assert(!m_samples || !m_size || m_samples->weights);
      }

      /**
       * @tparam Formula Concrete family selecting its own pool.
       * @tparam Key Owned constructor-argument tuple.
       * @returns Shared pool retained throughout static destruction.
       */
      template <class Formula, class Key>
      static Pool<Formula, Key>& getPool()
      {
        static auto* const pool = new Pool<Formula, Key>;
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

      /// @brief Concrete rule owned by the same permanent pool entry.
      std::reference_wrapper<const QuadratureFormulaBase> m_formula;

      /// @brief Sample count fixed when the interface is constructed.
      size_t m_size;

      /// @brief Borrowed arrays avoiding pointwise dispatch when available.
      Optional<Samples> m_samples;
  };

  template <class Formula>
  struct QuadratureFormula::Storage
  {
    /**
     * @brief Constructs the concrete rule before publishing its interface.
     * @tparam Args Constructor argument types.
     * @param args Arguments forwarded directly to the concrete rule.
     */
    template <class... Args>
    explicit Storage(Args&&... args)
      : generator(std::forward<Args>(args)...),
        formula(generator)
    {}

    /// @brief Immutable concrete rule owning its reference samples.
    const Formula generator;

    /// @brief Published interface borrowing this entry's concrete rule.
    const QuadratureFormula formula;
  };
}

#endif
