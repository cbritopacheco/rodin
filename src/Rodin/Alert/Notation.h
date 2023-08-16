#ifndef RODIN_ALERT_NOTATION_H
#define RODIN_ALERT_NOTATION_H

#include "Text.h"

namespace Rodin::Alert
{
  using NotationForeground = MagentaT;

  class NotationT : public Text<NotationForeground>
  {
    public:
      using Parent = Text<NotationForeground>;

      NotationT(const std::string& text)
        : Parent(text)
      {}
  };

  namespace Notation
  {
    static const NotationT Arrow("---->");

    inline
    NotationT Text(const std::string& text)
    {
      return text;
    }

    inline
    NotationT Incidence(size_t d, size_t dp)
    {
      std::stringstream ss;
      ss << "Incidence(" << d << " " << Arrow << " " << dp << ")";
      return ss.str();
    }

    inline
    NotationT Polytope(size_t d, Index index)
    {
      std::stringstream ss;
      ss << "Polytope(" << d << ", " << index << ")";
      return ss.str();
    }

    inline
    NotationT True()
    {
      return NotationT("true");
    }

    inline
    NotationT False()
    {
      return NotationT("false");
    }

    inline
    NotationT Predicate(bool v, const std::string& pred)
    {
      std::stringstream ss;
      ss << "[" << pred <<  "] = " << (v ? "true" : "false");
      return ss.str();
    }
  }

}

#endif
