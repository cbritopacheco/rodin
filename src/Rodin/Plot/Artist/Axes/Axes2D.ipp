#ifndef RODIN_PLOT_ARTIST_AXES_AXES2D_IPP
#define RODIN_PLOT_ARTIST_AXES_AXES2D_IPP

#include "Rodin/Plot/Artist/Figure.h"
#include "Rodin/Plot/Artist/Lines/Line2D.h"
#include "Rodin/Plot/Backend/Renderer/Drawables/Line2D.h"

#include "Axes2D.h"

namespace Rodin::Plot::Artist::Axes
{
  template <class T, class ... Args>
  Lines::Line2D& Axes2D::plot(
      const Eigen::ArrayX<T>& x,
      const Eigen::ArrayX<T>& y,
      Args&&... args)
  {
    setXLimits({x.minCoeff(), x.maxCoeff()});
    setYLimits({y.minCoeff(), y.maxCoeff()});

    return addArtist<Lines::Line2D>(
        x.template cast<float>(), y.template cast<float>(),
        std::forward<Args>(args)...);
  }
}

#endif
