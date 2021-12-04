/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cassert>
#include <functional>

#include "Rodin/Plot/Backend/Bases/BaseTopLevelArtist2D.h"

#include "BaseArtist2D.h"

namespace Rodin::Plot::Backend::Bases
{
  BaseArtist2D::BaseArtist2D(BaseArtist2D& parent)
    : m_parent(parent)
  {}

  BaseArtist2D& BaseArtist2D::setParent(BaseArtist2D& parent)
  {
    m_parent = parent;
    return *this;
  }

  BaseArtist2D& BaseArtist2D::getParent()
  {
    assert(m_parent.has_value());
    return *m_parent;
  }

  const BaseArtist2D& BaseArtist2D::getParent() const
  {
    assert(m_parent.has_value());
    return *m_parent;
  }

  bool BaseArtist2D::isTopLevel() const
  {
    return !m_parent.has_value();
  }

  Renderer::DrawableGroup2D& BaseArtist2D::getTopLevelDrawableGroup()
  {
    return getTopLevelArtist().getDrawableGroup();
  }

  const Renderer::DrawableGroup2D& BaseArtist2D::getTopLevelDrawableGroup() const
  {
    return getTopLevelArtist().getDrawableGroup();
  }

  BaseTopLevelArtist2D& BaseArtist2D::getTopLevelArtist()
  {
    // TODO: The artist tree is basically static, since when the parent goes
    // out of scope then it gets destroyed along with its children. So we
    // probably can do some kind of caching.
    if (isTopLevel())
      return static_cast<BaseTopLevelArtist2D&>(*this);
    else
    {
      std::reference_wrapper<BaseArtist2D> parent = getParent();
      while (!parent.get().isTopLevel())
        parent = parent.get().getParent();
      return static_cast<BaseTopLevelArtist2D&>(parent.get());
    }
  }

  const BaseTopLevelArtist2D& BaseArtist2D::getTopLevelArtist() const
  {
    // TODO: The artist tree is basically static, since when the parent goes
    // out of scope then it gets destroyed along with its children. So we
    // probably can do some kind of caching.
    if (isTopLevel())
      return static_cast<const BaseTopLevelArtist2D&>(*this);
    else
    {
      std::reference_wrapper<const BaseArtist2D> parent = getParent();
      while (!parent.get().isTopLevel())
        parent = parent.get().getParent();
      return static_cast<const BaseTopLevelArtist2D&>(parent.get());
    }
  }
}
