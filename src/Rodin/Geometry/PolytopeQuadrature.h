/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POLYTOPEQUADRATURE_H
#define RODIN_GEOMETRY_POLYTOPEQUADRATURE_H

/**
 * @file
 * @brief Cached quadrature data attached to mesh polytopes.
 *
 * This file defines:
 * - @ref Rodin::Geometry::PolytopeQuadrature, which stores the mapped
 *   quadrature points associated with a specific quadrature formula on a
 *   specific polytope.
 * - @ref Rodin::Geometry::PolytopeQuadratureIndex, a mesh-owned, thread-safe
 *   cache indexed by polytope dimension, polytope index, and quadrature
 *   formula identity.
 *
 * # Purpose
 *
 * Reference quadrature formulas belong to the @ref Rodin::QF module and live on
 * reference polytopes. In contrast, this file provides a geometry-side cache of
 * their realization on concrete mesh polytopes.
 *
 * For a polytope @f$ \tau @f$ and a reference quadrature formula
 * @f$ \{ (\hat x_q, w_q) \}_{q=1}^{n_q} @f$, the realized quadrature consists
 * of the mapped geometric points:
 *
 * @f[
 *   p_q = x_\tau(\hat x_q), \qquad q = 1, \dots, n_q,
 * @f]
 *
 * where @f$ x_\tau @f$ is the polytope transformation.
 *
 * # Cache organization
 *
 * The cache is organized as:
 *
 * @f[
 *   (d, i) \mapsto \text{Slot}_{d,i}
 * @f]
 *
 * where:
 * - @f$ d @f$ is the polytope dimension,
 * - @f$ i @f$ is the polytope index within that dimension.
 *
 * Each slot stores:
 * - one hot entry for the fastest repeated-hit access path,
 * - one small bounded ring buffer of owned cached quadratures.
 *
 * Lookup order is:
 * - hot entry,
 * - bounded per-slot cache,
 * - lazy construction on miss.
 *
 * Hits in the bounded cache are promoted to the hot slot. Misses are inserted
 * into the bounded cache in round-robin order and then promoted to the hot slot.
 *
 * # Key
 *
 * The cache key is the address of the quadrature formula object:
 *
 * @f[
 *   \texttt{const QF::QuadratureFormulaBase*}
 * @f]
 *
 * This design assumes that callers obtain quadrature formulas from a canonical
 * source such as @ref Rodin::QF::PolytopeQuadratureFormula, so that the
 * same logical formula has a stable object identity.
 */

#include <algorithm>
#include <atomic>
#include <array>
#include <cassert>
#include <deque>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

#include <boost/serialization/split_member.hpp>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Types.h"
#include "Rodin/QF/ForwardDecls.h"
#include "Point.h"

namespace Rodin::Geometry
{
  /**
   * @brief Cached quadrature attached to one specific polytope.
   *
   * A PolytopeQuadrature stores the mapped quadrature points corresponding to a
   * single quadrature formula on a single concrete polytope. It is:
   * - not a reference quadrature formula,
   * - not a variational quadrature rule,
   * - not a container of weights.
   *
   * Its purpose is to cache the geometric realization of a quadrature formula
   * on a polytope so that repeated integration passes can reuse mapped
   * @ref Rodin::Geometry::Point objects.
   */
  class PolytopeQuadrature
  {
    public:
      /**
       * @brief Constructs the realized quadrature on a polytope.
       * @param[in] polytope Concrete mesh polytope
       * @param[in] qf Reference quadrature formula
       */
      PolytopeQuadrature(
          const Polytope& polytope,
          const QF::QuadratureFormulaBase& qf);

      /**
       * @brief Gets the reference quadrature formula used to build this object.
       * @returns Reference to the quadrature formula
       */
      const QF::QuadratureFormulaBase& getQuadratureFormula() const
      {
        assert(m_qf);
        return *m_qf;
      }

      /**
       * @brief Gets the number of mapped quadrature points.
       * @returns Number of quadrature points
       */
      size_t getSize() const
      {
        return m_ps.size();
      }

      /**
       * @brief Gets a mapped quadrature point by index.
       * @param[in] i Quadrature point index
       * @returns Reference to the mapped point
       */
      const Point& getPoint(size_t i) const
      {
        return m_ps.at(i);
      }

    private:
      const QF::QuadratureFormulaBase* m_qf = nullptr; ///< Source quadrature formula
      std::vector<Point> m_ps;                         ///< Mapped quadrature points
  };

  /**
   * @brief Thread-safe cache of polytope quadratures.
   *
   * This class maintains a mesh-owned cache of @ref PolytopeQuadrature objects.
   * Within each dimension, it first resolves a quadrature formula and then
   * indexes a dense array by polytope index.
   *
   * # Cache design
   *
   * A formula cache owns one stable slot per polytope. Repeated assembly with a
   * fixed formula therefore resolves one atomically published formula pointer
   * and one atomically published quadrature pointer. Cached objects remain valid
   * until the mesh geometry is flushed.
   *
   * # Thread safety
   *
   * - Repeated cache hits are lock-free.
   * - Formula creation and storage growth are serialized per dimension.
   * - First quadrature construction is serialized through striped locks.
   *
   * # Lifecycle
   *
   * This cache is derived state. It should be cleared whenever the mesh
   * geometry changes.
   *
   * @note This cache is intentionally transient and is not serializable.
   */
  class PolytopeQuadratureIndex
  {
    public:
      friend class boost::serialization::access;

      /**
       * @brief Default constructor.
       */
      PolytopeQuadratureIndex() = default;

      /**
       * @brief Destructor.
       */
      ~PolytopeQuadratureIndex() = default;

      /**
       * @brief Copy constructor (deleted).
       */
      PolytopeQuadratureIndex(const PolytopeQuadratureIndex&) = delete;

      /**
       * @brief Copy assignment operator (deleted).
       */
      PolytopeQuadratureIndex& operator=(const PolytopeQuadratureIndex&) = delete;

      /**
       * @brief Move constructor.
       */
      PolytopeQuadratureIndex(PolytopeQuadratureIndex&& other) noexcept
        : m_dimensions(std::move(other.m_dimensions))
      {}

      /**
       * @brief Move assignment operator.
       */
      PolytopeQuadratureIndex& operator=(PolytopeQuadratureIndex&& other) noexcept
      {
        m_dimensions = std::move(other.m_dimensions);
        return *this;
      }

      /**
       * @brief Initializes the index for a mesh of given topological dimension.
       * @param[in] meshDim Topological mesh dimension
       *
       * Allocates dimension buckets for dimensions
       * @f$ 0, \dots, \texttt{meshDim} @f$.
       */
      void initialize(size_t meshDim)
      {
        m_dimensions.resize(meshDim + 1);
      }

      /**
       * @brief Resizes storage for polytopes of dimension @p d.
       * @param[in] d Polytope dimension
       * @param[in] count Number of polytope slots to reserve
       *
       * Ensures that the cache can store entries for at least @p count
       * polytopes of dimension @p d.
       *
       * @note The dimension bucket must already have been initialized.
       */
      void resize(size_t d, size_t count)
      {
        assert(d < m_dimensions.size());
        auto& dim = m_dimensions[d];
        std::lock_guard<std::mutex> lock(dim.mutex);
        for (auto& formula : dim.formulas)
          formula->resize(count);
      }

      /**
       * @brief Clears all cached quadratures.
       *
       * Releases all cached quadratures in all dimensions. The index remains
       * initialized, but contains no entries afterwards.
       */
      void clear()
      {
        for (auto& dim : m_dimensions)
        {
          std::lock_guard<std::mutex> lock(dim.mutex);
          dim.hot.store(nullptr, std::memory_order_relaxed);
          for (auto& entry : dim.lookup)
            entry.store(nullptr, std::memory_order_relaxed);
          dim.formulas.clear();
        }
      }

      /**
       * @brief Gets or creates a cached quadrature for a given polytope and formula.
       *
       * @tparam Factory Callable returning `std::unique_ptr<PolytopeQuadrature>`
       * @param[in] p Pair `(dimension, index)` identifying the polytope
       * @param[in] count Expected number of polytopes in that dimension
       * @param[in] qf Quadrature formula used as cache key
       * @param[in] factory Nullary factory called only on cache miss
       * @returns Reference to the cached polytope quadrature
       *
       * Formula and quadrature pointers are resolved atomically on repeated
       * accesses. A miss creates stable owned storage under a narrow lock.
       *
       * @note Cache mutation, including clear(), must not overlap use of a
       *       returned reference. Mesh geometry must remain unchanged during
       *       variational evaluation.
       *
       * @throws Rodin::Alert::Exception if:
       * - the dimension is outside the initialized range,
       * - the polytope index exceeds the supplied count.
       */
      template <class Factory>
      const PolytopeQuadrature& get(
          const std::pair<size_t, Index>& p,
          size_t count,
          const QF::QuadratureFormulaBase& qf,
          const Factory& factory) const
      {
        const size_t d = p.first;
        const Index idx = p.second;

        if (m_dimensions.size() <= d)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Requested dimension " << d
            << " is out of range for an index initialized with "
            << m_dimensions.size() << " dimension buckets."
            << Alert::Raise;
        }

        if (idx >= count)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Polytope index " << idx
            << " is out of range for dimension " << d
            << " with expected count " << count << "."
            << Alert::Raise;
        }

        auto& dim = m_dimensions[d];
        Formula* formula = dim.hot.load(std::memory_order_acquire);
        if (!formula || formula->qf != &qf)
        {
          formula = nullptr;
          for (const auto& entry : dim.lookup)
          {
            Formula* candidate = entry.load(std::memory_order_acquire);
            if (candidate && candidate->qf == &qf)
            {
              formula = candidate;
              break;
            }
          }

          if (!formula)
          {
            std::lock_guard<std::mutex> lock(dim.mutex);
            for (const auto& candidate : dim.formulas)
              if (candidate->qf == &qf)
              {
                formula = candidate.get();
                break;
              }

            if (!formula)
            {
              dim.formulas.emplace_back(std::make_unique<Formula>(&qf, count));
              formula = dim.formulas.back().get();
              for (auto& entry : dim.lookup)
                if (!entry.load(std::memory_order_relaxed))
                {
                  entry.store(formula, std::memory_order_release);
                  break;
                }
            }
          }
          dim.hot.store(formula, std::memory_order_release);
        }

        if (idx >= formula->publishedSize.load(std::memory_order_acquire))
          formula->resize(count);

        auto& slot = formula->slots[idx];
        if (auto* quadrature = slot.ptr.load(std::memory_order_acquire))
          return *quadrature;

        std::lock_guard<std::mutex> lock(
          formula->mutexes[static_cast<size_t>(idx) % formula->mutexes.size()]);
        PolytopeQuadrature* quadrature = slot.ptr.load(std::memory_order_relaxed);
        if (!quadrature)
        {
          slot.owner = factory();
          quadrature = slot.owner.get();
          slot.ptr.store(quadrature, std::memory_order_release);
        }
        return *quadrature;
      }

      /**
       * @brief Gets the number of initialized dimension buckets.
       * @returns Number of dimension buckets
       */
      size_t dimensions() const
      {
        return m_dimensions.size();
      }

    private:
      struct Slot
      {
          std::atomic<PolytopeQuadrature*> ptr{nullptr};
          std::unique_ptr<PolytopeQuadrature> owner;
      };

      struct Formula
      {
          static constexpr size_t MutexCount = 64;

          Formula(const QF::QuadratureFormulaBase* qf, size_t count)
            : qf(qf),
              slots(count),
              publishedSize(count)
          {}

          void resize(size_t count)
          {
            if (count <= publishedSize.load(std::memory_order_acquire))
              return;
            std::lock_guard<std::mutex> lock(resizeMutex);
            if (slots.size() < count)
              slots.resize(count);
            publishedSize.store(slots.size(), std::memory_order_release);
        }

        const QF::QuadratureFormulaBase* qf;
        std::deque<Slot> slots;
        std::atomic<size_t> publishedSize;
        std::mutex resizeMutex;
        std::array<std::mutex, MutexCount> mutexes;
      };

      struct Dimension
      {
          static constexpr size_t LookupCapacity = 16;

          Dimension() = default;
          Dimension(const Dimension&) = delete;
          Dimension& operator=(const Dimension&) = delete;

          Dimension(Dimension&& other) noexcept
            : formulas(std::move(other.formulas))
          {
            rebuildLookup();
          }

        Dimension& operator=(Dimension&& other) noexcept
        {
          formulas = std::move(other.formulas);
          rebuildLookup();
          return *this;
        }

        void rebuildLookup()
        {
          hot.store(nullptr, std::memory_order_relaxed);
          for (auto& entry : lookup)
            entry.store(nullptr, std::memory_order_relaxed);
          const size_t count = std::min(formulas.size(), lookup.size());
          for (size_t i = 0; i < count; ++i)
            lookup[i].store(formulas[i].get(), std::memory_order_relaxed);
        }

        std::vector<std::unique_ptr<Formula>> formulas;
        std::array<std::atomic<Formula*>, LookupCapacity> lookup{};
        std::atomic<Formula*> hot{nullptr};
        mutable std::mutex mutex;
      };

      template <class Archive>
      void save(Archive&, const unsigned int) const
      {}

      template <class Archive>
      void load(Archive&, const unsigned int)
      {
        clear();
      }

      BOOST_SERIALIZATION_SPLIT_MEMBER()

      mutable std::vector<Dimension> m_dimensions; ///< Dimension-partitioned storage
  };
}

#endif
