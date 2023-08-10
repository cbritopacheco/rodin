/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_TEXT_H
#define RODIN_ALERT_TEXT_H

#include <string>
#include <termcolor/termcolor.hpp>

#include "Rodin/Types.h"

#include "ForwardDecls.h"
#include "Color.h"
#include "Reset.h"

namespace Rodin::Alert
{
  enum class Attribute
  {
    Bold,
    Dark,
    Italic,
    Underline,
    Blink,
    Reverse,
    Concealed,
    Crossed
  };

  inline
  std::ostream& operator<<(std::ostream& os, Attribute a)
  {
    switch (a)
    {
      case Attribute::Bold:
      {
        os << termcolor::bold;
        break;
      }
      case Attribute::Dark:
      {
        os << termcolor::dark;
        break;
      }
      case Attribute::Italic:
      {
        os << termcolor::italic;
        break;
      }
      case Attribute::Underline:
      {
        os << termcolor::underline;
        break;
      }
      case Attribute::Blink:
      {
        os << termcolor::blink;
        break;
      }
      case Attribute::Reverse:
      {
        os << termcolor::reverse;
        break;
      }
      case Attribute::Concealed:
      {
        os << termcolor::concealed;
        break;
      }
      case Attribute::Crossed:
      {
        os << termcolor::crossed;
        break;
      }
    }
    return os;
  }

  template <class F = NoColorT, class B = NoColorT>
  class Text
  {
    public:
      using Foreground = F;

      using Background = B;

      Text(const std::string& text)
        : m_string(text)
      {}

      Text(const std::string& text, const Foreground& fg)
        : m_string(text), m_fg(fg)
      {}

      Text(const std::string& text, const Foreground& fg, const Background& bg)
        : m_string(text), m_fg(fg), m_bg(bg)
      {}

      constexpr
      Text(const Text&) = default;

      constexpr
      Text(Text&&) = default;

      inline
      const std::string& getString() const
      {
        return m_string;
      }

      inline
      const Foreground& getForeground() const
      {
        return m_fg;
      }

      inline
      const Background& getBackground() const
      {
        return m_bg;
      }

      Text& setAttribute(Attribute a)
      {
        m_attributes.insert(a);
        return *this;
      }

      Text& setBold()
      {
        setAttribute(Attribute::Bold);
        return *this;
      }

      Text& setDark()
      {
        setAttribute(Attribute::Bold);
        return *this;
      }

      Text& setItalic()
      {
        setAttribute(Attribute::Italic);
        return *this;
      }

      Text& setUnderline()
      {
        setAttribute(Attribute::Underline);
        return *this;
      }

      Text& setBlink()
      {
        setAttribute(Attribute::Blink);
        return *this;
      }

      Text& setReverse()
      {
        setAttribute(Attribute::Reverse);
        return *this;
      }

      Text& setConcealed()
      {
        setAttribute(Attribute::Concealed);
        return *this;
      }

      Text& setCrossed()
      {
        setAttribute(Attribute::Crossed);
        return *this;
      }

      const FlatSet<Attribute>& getAttributes() const
      {
        return m_attributes;
      }

    private:
      std::string m_string;
      Foreground m_fg;
      Background m_bg;
      FlatSet<Attribute> m_attributes;
  };

  Text(const std::string&) -> Text<NoColorT, NoColorT>;

  template <class F>
  Text(const std::string&, const F&) -> Text<F, NoColorT>;

  template <class F, class B>
  Text(const std::string&, const F&, const B&) -> Text<F, B>;

  template <class Foreground, class Background>
  std::ostream& operator<<(std::ostream& os, const Text<Foreground, Background>& text)
  {
    os << Reset;
    os << text.getForeground() << text.getBackground();
    for (const auto& a : text.getAttributes())
      os << a;
    os <<  text.getString() << Reset;
    return os;
  }

  template <class Foreground>
  std::ostream& operator<<(std::ostream& os, const Text<Foreground, NoColorT>& text)
  {
    os << Reset;
    os << text.getForeground();
    for (const auto& a : text.getAttributes())
      os << a;
    os <<  text.getString() << Reset;
    return os;
  }

  template <class Background>
  std::ostream& operator<<(std::ostream& os, const Text<NoColorT, Background>& text)
  {
    os << Reset;
    os << text.getBackground();
    for (const auto& a : text.getAttributes())
      os << a;
    os <<  text.getString() << Reset;
    return os;
  }

  inline
  std::ostream& operator<<(std::ostream& os, const Text<NoColorT, NoColorT>& text)
  {
    os << Reset;
    for (const auto& a : text.getAttributes())
      os << a;
    os <<  text.getString() << Reset;
    return os;
  }
}

#endif

