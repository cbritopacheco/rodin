#ifndef RODIN_MODELS_DITANCE_BASE_H
#define RODIN_MODELS_DITANCE_BASE_H

#include "Rodin/Geometry/Types.h"

namespace Rodin::Models::Distance
{
  template <class Derived>
  class Base
  {
    public:
      template <class A1, class ... As>
      Derived& setInterior(A1&& a1, As&& ... as)
      {
        m_interior =
          FlatSet<Geometry::Attribute>{std::forward<A1>(a1), std::forward<As>(as)...};
        return static_cast<Derived&>(*this);
      }

      Derived& setInterior(const FlatSet<Geometry::Attribute>& interior)
      {
        m_interior = interior;
        return static_cast<Derived&>(*this);
      }

      template <class A1, class ... As>
      Derived& setInterface(A1&& a1, As&& ... as)
      {
        m_interface =
          FlatSet<Geometry::Attribute>{std::forward<A1>(a1), std::forward<As>(as)...};
        return static_cast<Derived&>(*this);
      }

      Derived& setInterface(const FlatSet<Geometry::Attribute>& interface)
      {
        m_interface = interface;
        return static_cast<Derived&>(*this);
      }

      const auto& getInterior() const
      {
        return m_interior;
      }

      const auto& getInterface() const
      {
        return m_interface;
      }

    private:
      FlatSet<Geometry::Attribute> m_interior;
      FlatSet<Geometry::Attribute> m_interface;
  };
}

#endif

