/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENTITERATOR_H
#define RODIN_VARIATIONAL_FINITEELEMENTITERATOR_H

#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
  class FiniteElementIterator
  {
    public:
      FiniteElementIterator(const Geometry::Simplex& mesh);

      FiniteElementIterator(const FiniteElementIterator&) = delete;

      FiniteElementIterator(FiniteElementIterator&&) = default;

      bool end() const;

      FiniteElementIterator& operator++();

      FiniteElement& operator*() const noexcept;

      FiniteElement* operator->() const noexcept;

      size_t getDimension() const
      {
        return m_dimension;
      }

    private:
      FiniteElement* generate() const;

      const size_t m_dimension;
      mutable bool m_dirty;
      mutable std::unique_ptr<FiniteElement> m_FiniteElement;
  };
}

#endif

