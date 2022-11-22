/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTEGRATOR_H
#define RODIN_VARIATIONAL_INTEGRATOR_H

#include <set>

#include "Rodin/Geometry/ForwardDecls.h"

#include "Rodin/FormLanguage/Base.h"


namespace Rodin::Variational
{
   class Integrator : public FormLanguage::Base
   {
      public:
         using Parent = FormLanguage::Base;

         enum class Type
         {
            Linear,
            Bilinear
         };

         Integrator() = default;

         Integrator(const Integrator& other)
            : Parent(other)
         {}

         Integrator(Integrator&& other)
            : Parent(std::move(other))
         {}

         virtual Type getType() const = 0;

         virtual Geometry::Region getRegion() const = 0;

         virtual Integrator* copy() const noexcept override = 0;
   };
}

#endif

