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
  }

}

#endif
