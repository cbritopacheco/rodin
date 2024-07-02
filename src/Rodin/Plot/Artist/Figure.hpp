#ifndef RODIN_PLOT_ARTIST_FIGURE_IPP
#define RODIN_PLOT_ARTIST_FIGURE_IPP

#include "Rodin/Plot/Artist/Axes/Axes2D.h"
#include "Rodin/Plot/Backend/Renderer/Drawables/Frame.h"

#include "Figure.h"

namespace Rodin::Plot::Artist
{
  template <class AxesT, class ... Args>
  AxesT& Figure::addAxes(Args&&... args)
  {
    static_assert(std::is_base_of<Backend::Bases::BaseAxes, AxesT>::value,
        "Axes must be derived from BaseAxes");
    auto ax = std::make_unique<AxesT>(*this, std::forward<Args>(args)...);
    auto& frame = draw<Backend::Renderer::Drawables::Frame>(
        Magnum::Math::Vector<Real>2<float>(ax->getBottomLeft()),
        Magnum::Math::Vector<Real>2<float>(ax->getSize()));
    auto wrapper = AxesWrapper{std::move(ax), std::ref(frame)};
    m_axes.push_back(std::move(wrapper));
    return static_cast<AxesT&>(*(m_axes.back().axes));
  }

  template <class AxesT>
  AxesT& Figure::addAxes()
  {
    auto size = getWindowSize();
    return addAxes<AxesT>(
        Magnum::Math::Vector<Real>2<int>(0.1 * size.x(), 0.1 * size.y()),
        Magnum::Math::Vector<Real>2<int>(0.8 * size.x(), 0.8 * size.y())
      );
  }
}

#endif
