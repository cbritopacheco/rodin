/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_INDEXGENERATOR_H
#define RODIN_GEOMETRY_INDEXGENERATOR_H

/**
 * @file
 * @brief Index generation utilities for iterating over polytope indices.
 */

#include <set>
#include <memory>
#include <utility>

#include "Rodin/Copyable.h"
#include "Rodin/Moveable.h"

#include "ForwardDecls.h"
#include "Polytope.h"
#include "Types.h"

namespace Rodin
{
  /**
   * @brief Sentinel type for default-constructed iterators.
   */
  struct DefaultSentinelT {};
  
  /**
   * @brief Default sentinel value.
   */
  inline constexpr DefaultSentinelT DefaultSentinel;
}

namespace Rodin::Geometry
{
  /**
   * @brief Abstract base class for index generators.
   *
   * Index generators produce sequences of polytope indices for iteration.
   * They support standard iterator operations (increment, dereference, end
   * check) and can be copied or moved.
   *
   * @see EmptyIndexGenerator, BoundedIndexGenerator, IteratorIndexGenerator
   */
  class IndexGeneratorBase : public Copyable, public Moveable
  {
    public:
      /**
       * @brief Virtual destructor.
       */
      virtual ~IndexGeneratorBase() = default;
      
      /**
       * @brief Checks if the generator has reached the end.
       * @returns True if at the end, false otherwise
       */
      virtual bool end() const = 0;
      
      /**
       * @brief Advances to the next index.
       * @returns Reference to this generator
       */
      virtual IndexGeneratorBase& operator++() = 0;
      
      /**
       * @brief Dereferences to get the current index.
       * @returns Current index value
       */
      virtual Index operator*() const noexcept = 0;
      
      /**
       * @brief Creates a copy of this generator.
       * @returns Pointer to a new generator
       */
      virtual IndexGeneratorBase* copy() const noexcept override = 0;
      
      /**
       * @brief Creates a moved copy of this generator.
       * @returns Pointer to a new generator
       */
      virtual IndexGeneratorBase* move() noexcept override = 0;
  };

  /**
   * @brief Index generator that represents an empty sequence.
   *
   * This generator immediately reports end() as true and produces no indices.
   * It is useful as a placeholder or for representing empty sets.
   */
  class EmptyIndexGenerator final : public IndexGeneratorBase
  {
    public:
      /**
       * @brief Default constructor.
       */
      constexpr EmptyIndexGenerator() = default;

      /**
       * @brief Move constructor.
       */
      constexpr EmptyIndexGenerator(EmptyIndexGenerator&& other)
        :  IndexGeneratorBase(std::move(other))
      {}

      /**
       * @brief Copy constructor.
       */
      constexpr EmptyIndexGenerator(const EmptyIndexGenerator& other)
        :  IndexGeneratorBase(other)
      {}

      /**
       * @brief Checks if at end (always true).
       * @returns True
       */
      bool end() const override
      {
        return true;
      }

      /**
       * @brief Increment operator (should not be called).
       * @returns Reference to this generator
       */
      EmptyIndexGenerator& operator++() override
      {
        assert(false);
        return *this;
      }

      /**
       * @brief Dereference operator (should not be called).
       * @returns Dummy index value
       */
      Index operator*() const noexcept override
      {
        assert(false);
        return 0;
      }

      /**
       * @brief Creates a copy of this generator.
       * @returns Pointer to a new EmptyIndexGenerator
       */
      EmptyIndexGenerator* copy() const noexcept override
      {
        return new EmptyIndexGenerator(*this);
      }

      /**
       * @brief Creates a moved copy of this generator.
       * @returns Pointer to a new EmptyIndexGenerator
       */
      EmptyIndexGenerator* move() noexcept override
      {
        return new EmptyIndexGenerator(std::move(*this));
      }
  };

  /**
   * @brief Index generator for a bounded range of indices.
   *
   * Generates indices in the range [start, end), incrementing by 1 each step.
   * This is useful for iterating over consecutive polytope indices.
   */
  class BoundedIndexGenerator final : public IndexGeneratorBase
  {
    public:
      /**
       * @brief Constructs a generator for the range [start, end).
       * @param[in] start First index (inclusive)
       * @param[in] end Last index (exclusive)
       */
      constexpr
      BoundedIndexGenerator(Index start, Index end)
        : m_start(start), m_end(end), m_curr(start)
      {}

      /**
       * @brief Move constructor.
       */
      constexpr
      BoundedIndexGenerator(BoundedIndexGenerator&& other)
        :  IndexGeneratorBase(std::move(other)),
          m_start(other.m_start), m_end(other.m_end), m_curr(other.m_curr)
      {}

      /**
       * @brief Copy constructor.
       */
      constexpr
      BoundedIndexGenerator(const BoundedIndexGenerator& other)
        :  IndexGeneratorBase(other),
          m_start(other.m_start), m_end(other.m_end), m_curr(other.m_end)
      {}

      /**
       * @brief Checks if at end of range.
       * @returns True if current index equals end index
       */
      bool end() const override
      {
        return m_curr == m_end;
      }

      /**
       * @brief Advances to next index.
       * @returns Reference to this generator
       */
      BoundedIndexGenerator& operator++() override
      {
        ++m_curr;
        return *this;
      }

      /**
       * @brief Gets the current index.
       * @returns Current index value
       */
      Index operator*() const noexcept override
      {
        assert(!end());
        return m_curr;
      }

      /**
       * @brief Creates a copy of this generator.
       * @returns Pointer to a new BoundedIndexGenerator
       */
      BoundedIndexGenerator* copy() const noexcept override
      {
        return new BoundedIndexGenerator(*this);
      }

      /**
       * @brief Creates a moved copy of this generator.
       * @returns Pointer to a new BoundedIndexGenerator
       */
      BoundedIndexGenerator* move() noexcept override
      {
        return new BoundedIndexGenerator(std::move(*this));
      }

    private:
      const Index m_start;
      const Index m_end;
      Index m_curr;
  };

  /**
   * @brief Index generator that wraps an iterator pair.
   *
   * Adapts any iterator type to the IndexGeneratorBase interface.
   * The iterator must dereference to an Index value.
   *
   * @tparam Iterator Type of the wrapped iterator
   */
  template <class Iterator>
  class IteratorIndexGenerator : public IndexGeneratorBase
  {
    public:
      /**
       * @brief Constructs from an iterator range.
       * @param[in] it Begin iterator
       * @param[in] end End iterator
       */
      IteratorIndexGenerator(Iterator it, Iterator end)
        : m_it(it), m_end(end)
      {}

      /**
       * @brief Move constructor.
       */
      IteratorIndexGenerator(IteratorIndexGenerator&& other)
        : IndexGeneratorBase(std::move(other)),
          m_it(std::move(other.m_it)), m_end(std::move(other.m_end))
      {}

      /**
       * @brief Copy constructor.
       */
      IteratorIndexGenerator(const IteratorIndexGenerator& other)
        : IndexGeneratorBase(other),
          m_it(other.m_it), m_end(other.m_end)
      {}

      /**
       * @brief Checks if iterator has reached end.
       * @returns True if at end
       */
      bool end() const override
      {
        return m_it == m_end;
      }

      /**
       * @brief Advances the iterator.
       * @returns Reference to this generator
       */
      IteratorIndexGenerator& operator++() override
      {
        ++m_it;
        return *this;
      }

      /**
       * @brief Gets the current index.
       * @returns Index value from dereferencing the iterator
       */
      Index operator*() const noexcept override
      {
        assert(!end());
        return *m_it;
      }

      /**
       * @brief Creates a copy of this generator.
       * @returns Pointer to a new IteratorIndexGenerator
       */
      IteratorIndexGenerator* copy() const noexcept override
      {
        return new IteratorIndexGenerator(*this);
      }

      /**
       * @brief Creates a moved copy of this generator.
       * @returns Pointer to a new IteratorIndexGenerator
       */
      IteratorIndexGenerator* move() noexcept override
      {
        return new IteratorIndexGenerator(std::move(*this));
      }

    private:
      Iterator m_it;
      Iterator m_end;
  };

  /**
   * @brief Index generator backed by a vector of indices.
   *
   * Owns a vector of indices and iterates over them. The indices can be
   * in any order and need not be consecutive.
   */
  class VectorIndexGenerator : public IndexGeneratorBase
  {
    public:
      /**
       * @brief Constructs from a vector of indices (move).
       * @param[in] indices Vector of index values
       */
      VectorIndexGenerator(std::vector<Index>&& indices)
        : m_indices(std::move(indices)),
          m_it(m_indices.begin())
      {}

      /**
       * @brief Move constructor.
       */
      VectorIndexGenerator(VectorIndexGenerator&& other)
        : IndexGeneratorBase(std::move(other)),
          m_indices(std::move(other.m_indices)),
          m_it(std::move(other.m_it))
      {}

      /**
       * @brief Copy constructor.
       */
      VectorIndexGenerator(const VectorIndexGenerator& other)
        : IndexGeneratorBase(other),
          m_indices(other.m_indices),
          m_it(other.m_it)
      {}

      /**
       * @brief Checks if at end of vector.
       * @returns True if iterator equals end()
       */
      bool end() const override
      {
        return m_it == m_indices.end();
      }

      /**
       * @brief Advances to next index in vector.
       * @returns Reference to this generator
       */
      VectorIndexGenerator& operator++() override
      {
        ++m_it;
        return *this;
      }

      /**
       * @brief Gets the current index.
       * @returns Current index value from vector
       */
      Index operator*() const noexcept override
      {
        assert(!end());
        return *m_it;
      }

      /**
       * @brief Creates a copy of this generator.
       * @returns Pointer to a new VectorIndexGenerator
       */
      VectorIndexGenerator* copy() const noexcept override
      {
        return new VectorIndexGenerator(*this);
      }

      /**
       * @brief Creates a moved copy of this generator.
       * @returns Pointer to a new VectorIndexGenerator
       */
      VectorIndexGenerator* move() noexcept override
      {
        return new VectorIndexGenerator(std::move(*this));
      }

    private:
      std::vector<Index> m_indices;
      std::vector<Index>::iterator m_it;
  };

  /**
   * @brief Index generator backed by a set of indices.
   *
   * Owns a set of indices and iterates over them in sorted order.
   * Automatically eliminates duplicate indices.
   */
  class SetIndexGenerator : public IndexGeneratorBase
  {
    public:
      /**
       * @brief Constructs from a set of indices (copy).
       * @param[in] indices Set of index values
       */
      SetIndexGenerator(const IndexSet& indices)
        : m_indices(indices),
          m_it(m_indices.begin())
      {}

      /**
       * @brief Constructs from a set of indices (move).
       * @param[in] indices Set of index values
       */
      SetIndexGenerator(IndexSet&& indices)
        : m_indices(std::move(indices)),
          m_it(m_indices.begin())
      {}

      /**
       * @brief Move constructor.
       */
      SetIndexGenerator(SetIndexGenerator&& other)
        : IndexGeneratorBase(std::move(other)),
          m_indices(std::move(other.m_indices)),
          m_it(std::move(other.m_it))
      {}

      /**
       * @brief Copy constructor.
       */
      SetIndexGenerator(const SetIndexGenerator& other)
        : IndexGeneratorBase(other),
          m_indices(other.m_indices),
          m_it(other.m_it)
      {}

      /**
       * @brief Checks if at end of set.
       * @returns True if iterator equals end()
       */
      bool end() const override
      {
        return m_it == m_indices.end();
      }

      /**
       * @brief Advances to next index in set.
       * @returns Reference to this generator
       */
      SetIndexGenerator& operator++() override
      {
        ++m_it;
        return *this;
      }

      /**
       * @brief Gets the current index.
       * @returns Current index value from set
       */
      Index operator*() const noexcept override
      {
        assert(!end());
        return *m_it;
      }

      /**
       * @brief Creates a copy of this generator.
       * @returns Pointer to a new SetIndexGenerator
       */
      SetIndexGenerator* copy() const noexcept override
      {
        return new SetIndexGenerator(*this);
      }

      /**
       * @brief Creates a moved copy of this generator.
       * @returns Pointer to a new SetIndexGenerator
       */
      SetIndexGenerator* move() noexcept override
      {
        return new SetIndexGenerator(std::move(*this));
      }

    private:
      IndexSet m_indices;
      IndexSet::iterator m_it;
  };
}

#endif
