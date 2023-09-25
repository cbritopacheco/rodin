/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PLOT_BACKEND_BASES_BASEARTIST2D_H
#define RODIN_PLOT_BACKEND_BASES_BASEARTIST2D_H

#include <utility>
#include <optional>
#include <functional>

#include "Rodin/Plot/ForwardDecls.h"
#include "Rodin/Plot/Backend/Renderer/Common.h"

namespace Rodin::Plot::Backend::Bases
{
  /**
   * @brief Base class for artists in 2D.
   */
  class BaseArtist2D
  {
    public:
      BaseArtist2D() = default;

      BaseArtist2D(BaseArtist2D& parent);

      virtual ~BaseArtist2D() = default;

      template <class Drawable, class ... Args>
      Drawable& draw(Args&& ... args);

      template <class Artist, class ... Args>
      Artist& addArtist(Args&&... args);

      BaseArtist2D& setParent(BaseArtist2D& parent);
      Backend::Bases::BaseArtist2D& getParent();
      const Backend::Bases::BaseArtist2D& getParent() const;

      virtual bool isTopLevel() const;
      BaseTopLevelArtist2D& getTopLevelArtist();
      const BaseTopLevelArtist2D& getTopLevelArtist() const;

      virtual Renderer::Object2D& getObject2D() = 0;

    private:
      Renderer::DrawableGroup2D& getTopLevelDrawableGroup();
      const Renderer::DrawableGroup2D& getTopLevelDrawableGroup() const;

      std::optional<std::reference_wrapper<BaseArtist2D>> m_parent;
      std::vector<std::unique_ptr<BaseArtist2D>> m_artists;
      std::vector<BaseDrawable2D*> m_drawables;
  };
}

#include "BaseArtist2D.hpp"

#endif
