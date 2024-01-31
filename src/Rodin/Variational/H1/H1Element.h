/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H1ELEMENT_H
#define RODIN_VARIATIONAL_H1_H1ELEMENT_H

/**
 * @file
 * @brief Header contatining definitions for the class P1Element.
 */

/**
 * @ingroup RodinDirectives
 * @brief Indicates the maximum vector dimension a P1Element
 */
#define RODIN_H1_MAX_VECTOR_DIMENSION 16

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/FormLanguage/ForwardDecls.h"
#include "Rodin/Variational/FiniteElement.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   */
  template <class Range>
  struct Traits<Variational::H1Element<Range>>
  {
    using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  template <>
  class H1Element<Scalar> : public FiniteElementBase<H1Element<Scalar>>
  {
    public:
      /// Parent class
      using Parent = FiniteElementBase<H1Element<Scalar>>;

      /// Type of range
      using RangeType = Scalar;

      H1Element() = default;

      constexpr
      H1Element(Geometry::Polytope::Type geometry)
        : Parent(geometry)
      {}

      constexpr
      H1Element(const H1Element& other)
        : Parent(other)
      {}

      constexpr
      H1Element(H1Element&& other)
        : Parent(std::move(other))
      {}

      constexpr
      H1Element& operator=(const H1Element& other)
      {
        Parent::operator=(other);
        return *this;
      }

      constexpr
      H1Element& operator=(H1Element&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }
  };
}

#endif

