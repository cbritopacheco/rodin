/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Base class for integrators in finite element assembly.
 */

#ifndef RODIN_VARIATIONAL_INTEGRATOR_H
#define RODIN_VARIATIONAL_INTEGRATOR_H

#include "Rodin/FormLanguage/Base.h"

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Base class for integrators in finite element assembly.
   *
   * Integrators compute element-level contributions for forms (bilinear or linear).
   */
  class Integrator : public FormLanguage::Base
  {
    public:
      /// @brief Parent class type
      using Parent = FormLanguage::Base;

      /**
       * @brief Integrator type enumeration.
       */
      enum class Type
      {
        Linear,   ///< Linear form integrator
        Bilinear  ///< Bilinear form integrator
      };

      Integrator() = default;

      Integrator(const Integrator& other)
        : Parent(other)
      {}

      Integrator(Integrator&& other)
        : Parent(std::move(other))
      {}

      virtual ~Integrator() = default;

      virtual Type getType() const = 0;

      virtual Integrator* copy() const noexcept override = 0;
  };
}

#endif

