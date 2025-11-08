#ifndef RODIN_GEOMETRY_POLYTOPETRANSFORMATIONINDEX_H
#define RODIN_GEOMETRY_POLYTOPETRANSFORMATIONINDEX_H

#include <atomic>
#include <deque>
#include <mutex>
#include <memory>
#include <vector>
#include <cassert>
#include <utility>
#include <shared_mutex>

#include <boost/serialization/access.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/vector.hpp>

#include "Rodin/Types.h"
#include "Rodin/Geometry/PolytopeTransformation.h"

namespace Rodin::Geometry
{
  class PolytopeTransformationIndex
  {
    friend class boost::serialization::access;

    struct Slot
    {
      std::atomic<PolytopeTransformation*> ptr{nullptr};
      std::unique_ptr<PolytopeTransformation> owner;

      template <class Archive>
      void save(Archive& ar, const unsigned int) const
      {
        ar & owner; // polymorphic unique_ptr
      }

      template <class Archive>
      void load(Archive& ar, const unsigned int)
      {
        ar & owner;
        ptr.store(owner.get(), std::memory_order_relaxed);
      }

      BOOST_SERIALIZATION_SPLIT_MEMBER()
    };

    struct Dimension
    {
      std::deque<Slot> slots;
      mutable std::shared_mutex mutex;

      Dimension() = default;

      Dimension(const Dimension&) = delete;

      Dimension& operator=(const Dimension&) = delete;

      Dimension(Dimension&& other) noexcept
        : slots(std::move(other.slots)) {}

      Dimension& operator=(Dimension&& other) noexcept
      {
        slots = std::move(other.slots);
        return *this;
      }

      template <class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & slots; // mutex is not serialized
      }
    };

  public:
    PolytopeTransformationIndex() = default;

    ~PolytopeTransformationIndex() = default;

    PolytopeTransformationIndex(const PolytopeTransformationIndex&) = delete;

    PolytopeTransformationIndex& operator=(const PolytopeTransformationIndex&) = delete;

    PolytopeTransformationIndex(PolytopeTransformationIndex&& other) noexcept
      : m_dimensions(std::move(other.m_dimensions))
    {}

    PolytopeTransformationIndex& operator=(PolytopeTransformationIndex&& other) noexcept
    {
      m_dimensions = std::move(other.m_dimensions);
      return *this;
    }

    void initialize(size_t meshDim)
    {
      m_dimensions.resize(meshDim + 1);
    }

    void resize(size_t d, size_t count)
    {
      assert(d < m_dimensions.size());
      auto& dim = m_dimensions[d];
      std::unique_lock<std::shared_mutex> wr(dim.mutex);
      if (dim.slots.size() < count)
        dim.slots.resize(count);
    }

    void set(const std::pair<size_t, Index>& p,
             size_t count,
             std::unique_ptr<PolytopeTransformation> obj)
    {
      const size_t d   = p.first;
      const Index  idx = p.second;

      assert(idx < count);

      assert(d < m_dimensions.size());
      auto& dim = m_dimensions[d];
      std::unique_lock<std::shared_mutex> wr(dim.mutex);
      if (dim.slots.size() < count)
        dim.slots.resize(count);

      assert(idx < dim.slots.size());
      Slot& s = dim.slots[idx];
      s.owner = std::move(obj);
      s.ptr.store(s.owner.get(), std::memory_order_release);
    }

    template <class Factory>
    const PolytopeTransformation&
    get(const std::pair<size_t, Index>& p, size_t count, const Factory& factory) const
    {
      const size_t d   = p.first;
      const Index  idx = p.second;

      assert(d < m_dimensions.size());
      auto& dim = m_dimensions[d];

      // Fast path: read under shared lock
      {
        std::shared_lock<std::shared_mutex> rd(dim.mutex);
        if (idx < dim.slots.size())
        {
          assert(idx < count);
          const Slot& s = dim.slots[idx];
          if (auto* q = s.ptr.load(std::memory_order_acquire)) return *q;
        }
      }

      // Slow path: exclusive lock, ensure size, init in-place
      {
        std::unique_lock<std::shared_mutex> wr(dim.mutex);

        if (dim.slots.size() < count)
          dim.slots.resize(count);

        assert(idx < dim.slots.size());
        Slot& s = dim.slots[idx];
        PolytopeTransformation* q = s.ptr.load(std::memory_order_relaxed);
        if (!q)
        {
          auto up = factory(d, idx);
          s.owner = std::move(up);
          q = s.owner.get();
          s.ptr.store(q, std::memory_order_release);
        }
        return *q; // still under wr; safe
      }
    }

    void clear()
    {
      for (auto& dim : m_dimensions)
      {
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        dim.slots.clear();
      }
    }

    size_t dimensions() const
    {
      return m_dimensions.size();
    }

    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
      ar & m_dimensions;
    }

    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
      clear();
      ar & m_dimensions;
      // ptr fields restored by Slot::load
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

  private:
    mutable std::vector<Dimension> m_dimensions;
  };
} // namespace Rodin::Geometry

#endif // RODIN_GEOMETRY_POLYTOPETRANSFORMATIONINDEX_H
