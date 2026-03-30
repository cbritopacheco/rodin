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
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/QF/QuadratureFormula.h"
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
          throw std::out_of_range("PolytopeQuadratureIndex: dimension out of range.");
        }

        auto& dim = m_dimensions[d];

        {
          std::shared_lock<std::shared_mutex> rd(dim.mutex);
          if (idx < dim.slots.size())
          {
            const auto& slot = dim.slots[idx];
            std::shared_lock<std::shared_mutex> rdslot(slot.mutex);
            if (auto it = slot.quadratures.find(&qf); it != slot.quadratures.end())
              return *it->second;
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
          auto it = slot.quadratures.find(&qf);
          if (it == slot.quadratures.end())
          {
            auto up = factory();
            it = slot.quadratures.emplace(&qf, std::move(up)).first;
          }
          return *it->second;
        }
      }

    private:
      struct Slot
      {
        FlatMap<const QF::QuadratureFormulaBase*, std::unique_ptr<PolytopeQuadrature>> quadratures;
        mutable std::shared_mutex mutex;
      };

      struct Dimension
      {
        std::vector<Slot> slots;
        mutable std::shared_mutex mutex;
      };

      std::vector<Dimension> m_dimensions;
  };
}

#endif
