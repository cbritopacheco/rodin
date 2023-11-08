#ifndef RODIN_ALERT_NOTATION_H
#define RODIN_ALERT_NOTATION_H

#include "Text.h"

namespace Rodin::Alert
{
  using NotationForeground = MagentaT;

  class Notation : public Text<NotationForeground>
  {
    public:
      using Parent = Text<NotationForeground>;

      static const Notation Arrow;

      Notation(const std::string& text)
        : Parent(text)
      {}

      Notation(const char* out)
        : Parent(std::string(out))
      {}

      template <class T>
      Notation(const T& out)
        : Parent(std::to_string(out))
      {}

      Notation(const Notation& other)
        : Parent(other)
      {}

      Notation(Notation&& other)
        : Parent(std::move(other))
      {}

      inline
      static Notation Text(const std::string& text)
      {
        return text;
      }

      template <class T>
      inline
      static Notation Number(const T& out)
      {
        static_assert(std::is_arithmetic_v<T>);
        std::stringstream ss;
        ss << out;
        return ss.str();
      }

      template <class T>
      inline
      static Notation Print(const T& out)
      {
        std::stringstream ss;
        ss << out;
        return ss.str();
      }

      inline
      static Notation Incidence(size_t d, size_t dp)
      {
        std::stringstream ss;
        ss << "Incidence(" << d << " " << Arrow << " " << dp << ")";
        return ss.str();
      }

      inline
      static Notation Polytope(size_t d, Index index)
      {
        std::stringstream ss;
        ss << "Polytope(" << d << ", " << index << ")";
        return ss.str();
      }

      inline
      static Notation True()
      {
        return Notation("true");
      }

      inline
      static Notation False()
      {
        return Notation("false");
      }

      inline
      static Notation Predicate(bool v, const std::string& pred)
      {
        std::stringstream ss;
        ss << "[" << pred <<  "] = " << (v ? "true" : "false");
        return ss.str();
      }

      template <class Iterator>
      static Notation Set(Iterator first, Iterator last)
      {
        std::stringstream ss;
        ss << "{ ";
        for (Iterator it = first; it != last; ++it)
        {
          ss << *it;
          if (std::next(it) != last)
              ss << ", ";
        }
        ss << " }";
        return ss.str();
    }
  };
}

#endif
