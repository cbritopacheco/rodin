#ifndef RODIN_GEOMETRY_ATTRIBUTEINDEX_H
#define RODIN_GEOMETRY_ATTRIBUTEINDEX_H

#include <mutex>
#include <vector>
#include <shared_mutex>

#include <boost/serialization/access.hpp>

#include "Rodin/Geometry/Types.h" // Attribute, Index, FlatSet, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE

namespace Rodin::Geometry
{
  class AttributeIndex
  {
    friend class boost::serialization::access;

    struct Dimension
    {
      friend class boost::serialization::access;

      std::vector<Attribute> slots;

      mutable std::shared_mutex mutex;

      Dimension() = default;

      Dimension(const Dimension& other)
        : slots(other.slots)
      {}

      Dimension& operator=(const Dimension& other)
      {
        if (this != &other)
        {
          slots = other.slots;
        }
        return *this;
      }

      Dimension(Dimension&& other) noexcept
        : slots(std::move(other.slots))
      {}

      Dimension& operator=(Dimension&& other) noexcept
      {
        slots = std::move(other.slots);
        return *this;
      }

      template <class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & slots;
      }
    };

    public:
      AttributeIndex() = default;

      ~AttributeIndex() = default;

      AttributeIndex(const AttributeIndex& other)
        : m_dimensions(other.m_dimensions)
      {}

      AttributeIndex& operator=(const AttributeIndex&) = delete;

      AttributeIndex(AttributeIndex&& other) noexcept
        : m_dimensions(std::move(other.m_dimensions))
      {}

      AttributeIndex& operator=(AttributeIndex&& other) noexcept
      {
        m_dimensions = std::move(other.m_dimensions);
        return *this;
      }

      // Call once before concurrency.
      void initialize(size_t meshDim)
      {
        m_dimensions.resize(meshDim + 1);
      }

      void resize(size_t d, size_t count)
      {
        auto& dim = m_dimensions.at(d);
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() < count)
          dim.slots.resize(count, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
      }

      void set(const std::pair<size_t, Index>& p, Attribute attr)
      {
        const auto& [d, idx] = p;

        assert(d < m_dimensions.size());

        auto& dim = m_dimensions.at(d);
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() <= idx)
          dim.slots.resize(idx + 1, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
        assert(idx < dim.slots.size());
        dim.slots[idx] = attr;
      }

      void set(const std::pair<size_t, Index>& p, size_t count, Attribute attr)
      {
        const auto& [d, idx] = p;

        assert(d < m_dimensions.size());
        assert(idx < count);

        auto& dim = m_dimensions.at(d);
        std::unique_lock<std::shared_mutex> wr(dim.mutex);
        if (dim.slots.size() < count)
          dim.slots.resize(count, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
        dim.slots[idx] = attr;
      }

      // Return current value. Grow to 'count' if needed, then ensure idx exists.
      Attribute get(const std::pair<size_t, Index>& p, size_t count) const
      {
        const auto& [d, idx] = p;

        assert(d < m_dimensions.size());
        assert(idx < count);

        // Read phase
        {
          const Dimension& dim = m_dimensions[d];                  // const ref
          std::shared_lock<std::shared_mutex> rd(dim.mutex);
          if (idx < dim.slots.size())
            return dim.slots[idx];
        }

        // Write phase
        {
          Dimension& dim = m_dimensions[d];                        // non-const ref
          std::unique_lock<std::shared_mutex> wr(dim.mutex);
          if (dim.slots.size() < count)
            dim.slots.resize(count, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
          if (idx >= dim.slots.size())
            dim.slots.resize(idx + 1, RODIN_DEFAULT_POLYTOPE_ATTRIBUTE);
          return dim.slots[idx];
        }
      }

      FlatSet<Attribute> getAttributes(size_t d) const
      {
        FlatSet<Attribute> out;
        const auto& dim = m_dimensions.at(d);
        std::shared_lock<std::shared_mutex> rd(dim.mutex);
        for (Attribute a : dim.slots)
          out.insert(a);
        return out;
      }

      template <class Archive>
      void serialize(Archive& ar, const unsigned int)
      {
        ar & m_dimensions;
      }

    private:
      mutable std::vector<Dimension> m_dimensions;
  };
}
#endif
