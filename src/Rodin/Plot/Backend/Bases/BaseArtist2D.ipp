#ifndef RODIN_PLOT_BACKEND_BASES_BASEARTIST2D_IPP
#define RODIN_PLOT_BACKEND_BASES_BASEARTIST2D_IPP

#include "Rodin/Plot/Backend/Bases/BaseTopLevelArtist2D.h"

#include "BaseArtist2D.h"

namespace Rodin::Plot::Backend::Bases
{
  template <class Drawable, class ... Args>
  Drawable& BaseArtist2D::draw(Args&& ... args)
  {
    // TODO: Who owns the Drawable?
    // Probably the top level artist. A drawable can still be available
    // after an artist dies
    // So probably push back the drawable into the top level artist
    // Furthermore, the drawable should probably belong to the top level
    // artist drawable group
    m_drawables.push_back(
        new Drawable(
          getObject2D().addChild<Renderer::Object2D>(),
          &getTopLevelDrawableGroup(),
          std::forward<Args>(args)...)
        );
    return static_cast<Drawable&>(*m_drawables.back());
  }

  template <class Artist, class ... Args>
  Artist& BaseArtist2D::addArtist(Args&&... args)
  {
    m_artists.push_back(
        std::make_unique<Artist>(*this, std::forward<Args>(args)...));
    return static_cast<Artist&>(*m_artists.back());
  }
}

#endif
