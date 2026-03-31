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
 * @brief Quadrature cached on a specific polytope.
 */

#include <shared_mutex>
#include <mutex>
#include <memory>
#include <array>
#include <utility>
#include <vector>

#include "Rodin/Alert/Exception.h"
#include "Rodin/Types.h"
#include "Rodin/QF/ForwardDecls.h"
#include "Point.h"

namespace Rodin::Geometry
{
  /**
   * @brief Quadrature rule cached on a polytope.
   *
   * Stores the geometric points associated with a specific
   * quadrature formula on a specific polytope.
   */
  class PolytopeQuadrature
  {
    public:
      PolytopeQuadrature(const Polytope& polytope, const QF::QuadratureFormulaBase& qf);

      const QF::QuadratureFormulaBase& getQuadratureFormula() const
      {
        return *m_qf;
      }

      size_t getSize() const
      {
        return m_ps.size();
      }

      const Point& getPoint(size_t i) const
      {
        return m_ps.at(i);
      }

      const std::vector<Point>& getPoints() const
      {
        return m_ps;
      }

    private:
      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Point> m_ps;
  };

  /**
   * @brief Thread-safe cache of PolytopeQuadrature objects.
   */
  class PolytopeQuadratureIndex
  {
    public:
      PolytopeQuadratureIndex() = default;

      void initialize(size_t meshDim)
      {
        m_dimensions.resize(meshDim + 1);
      }

      void resize(size_t d, size_t count)
      {
        if (m_dimensions.size() <= d)
          m_dimensions.resize(d + 1);

        auto& dim = m_dimensions[d];
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() < count)
          dim.slots.resize(count);
      }

      void clear()
      {
        for (auto& dim : m_dimensions)
        {
          std::unique_lock<std::shared_mutex> wr(dim.mutex);
          dim.slots.clear();
        }
      }

      template <class Factory>
      /**
       * @brief Gets or creates a cached quadrature for a given polytope and formula.
       * @throws std::out_of_range if @p d is outside initialized dimensions.
       */
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
          Alert::Exception()
            << "PolytopeQuadratureIndex: dimension out of range."
            << Alert::Raise;
        }

        auto& dim = m_dimensions[d];

        {
          std::shared_lock<std::shared_mutex> rd(dim.mutex);
          if (idx < dim.slots.size())
          {
            const auto& slot = dim.slots[idx];
            std::shared_lock<std::shared_mutex> rdslot(slot.mutex);
            for (size_t i = 0; i < slot.size; ++i)
            {
              if (slot.quadratures[i].first == &qf)
                return *slot.quadratures[i].second;
            }
          }
        }

        {
          std::unique_lock<std::shared_mutex> wr(dim.mutex);
          if (dim.slots.size() < count)
            dim.slots.resize(count);
        }

        assert(idx < dim.slots.size());
        auto& slot = dim.slots[idx];

        {
          std::unique_lock<std::shared_mutex> wrslot(slot.mutex);
          for (size_t i = 0; i < slot.size; ++i)
          {
            if (slot.quadratures[i].first == &qf)
              return *slot.quadratures[i].second;
          }

          if (slot.size >= Slot::s_maxQuadratures)
          {
            Alert::Exception()
              << "PolytopeQuadratureIndex: exceeded per-polytope cache capacity."
              << Alert::Raise;
          }

          const size_t entry = slot.size++;
          slot.quadratures[entry] = std::make_pair(&qf, factory());
          return *slot.quadratures[entry].second;
        }
      }

    private:
      struct Slot
      {
        /// Maximum number of distinct quadrature formulas cached per polytope slot.
        /// A linear scan is intentionally used because this cache is capped to 8 entries.
        static constexpr size_t s_maxQuadratures = 8;
        using CacheEntry = std::pair<const QF::QuadratureFormulaBase*, std::unique_ptr<PolytopeQuadrature>>;

        // The mutex is intentionally default-constructed in moved instances.
        Slot() = default;
        Slot(const Slot&) = delete;
        Slot& operator=(const Slot&) = delete;
        Slot(Slot&& other) noexcept
          : quadratures(std::move(other.quadratures)),
            size(std::exchange(other.size, 0))
        {}
        Slot& operator=(Slot&& other) noexcept
        {
          quadratures = std::move(other.quadratures);
          size = std::exchange(other.size, 0);
          return *this;
        }

        std::array<CacheEntry, s_maxQuadratures> quadratures;
        size_t size = 0;
        mutable std::shared_mutex mutex;
      };

      struct Dimension
      {
        // The mutex is intentionally default-constructed in moved instances.
        Dimension() = default;
        Dimension(const Dimension&) = delete;
        Dimension& operator=(const Dimension&) = delete;
        Dimension(Dimension&& other) noexcept
          : slots(std::move(other.slots))
        {}
        Dimension& operator=(Dimension&& other) noexcept
        {
          slots = std::move(other.slots);
          return *this;
        }

        std::vector<Slot> slots;
        mutable std::shared_mutex mutex;
      };

      mutable std::vector<Dimension> m_dimensions;
  };
}

#endif
