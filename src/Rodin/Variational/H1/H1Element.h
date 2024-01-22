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

#include "Rodin/FormLanguage/ForwardDecls.h"

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
}

#endif

