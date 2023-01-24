/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_BACKEND_BASES_BASETOPLEVELARTIST2D_H
#define RODIN_BACKEND_BASES_BASETOPLEVELARTIST2D_H

#include "BaseArtist2D.h"

namespace Rodin::Plot::Backend::Bases
{
  class BaseTopLevelArtist2D : public BaseArtist2D
  {
   public:
   explicit BaseTopLevelArtist2D() = default;
   virtual ~BaseTopLevelArtist2D() = default;

   bool isTopLevel() const override;

   virtual Renderer::DrawableGroup2D& getDrawableGroup() = 0;
   virtual const Renderer::DrawableGroup2D& getDrawableGroup() const = 0;
  };
}

#endif
