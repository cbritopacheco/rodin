#ifndef RODIN_ALERT_IDENTIFIER_H
#define RODIN_ALERT_IDENTIFIER_H

#include "Text.h"

namespace Rodin::Alert
{
  template <class Foreground>
  class IdentifierT : public Text<Foreground>
  {
    public:
      using Parent = Text<Foreground>;

      IdentifierT(const std::string& text)
        : Parent(text)
      {
        this->setBold();
      }
  };

  namespace Identifier
  {
    inline
    IdentifierT<CyanT> Class(const std::string& id)
    {
      return id;
    }

    inline
    IdentifierT<CyanT> Namespace(const std::string& id)
    {
      return id;
    }

    inline
    IdentifierT<ResetT> Function(const std::string& id)
    {
      return id;
    }
  }
}

#endif

