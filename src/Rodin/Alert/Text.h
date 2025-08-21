/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_TEXT_H
#define RODIN_ALERT_TEXT_H

#include <string>
#include <cstdint>
#include <cstdlib>
#include <iostream>

#include <termcolor/termcolor.hpp>

#include "Rodin/Types.h"

#include "ForwardDecls.h"
#include "Color.h"
#include "Reset.h"

namespace Rodin::Alert
{
  /**
   * @brief Text formatting attributes for terminal output.
   * @ingroup AlertModule
   *
   * Enumeration of text formatting attributes that can be applied to terminal
   * output including bold, italic, underline, and other visual effects.
   */
  enum class Attribute
  {
    Bold,      ///< Bold text formatting
    Dark,      ///< Dark text formatting
    Italic,    ///< Italic text formatting  
    Underline, ///< Underlined text
    Blink,     ///< Blinking text effect
    Reverse,   ///< Reverse video (inverted colors)
    Concealed, ///< Concealed/hidden text
    Crossed    ///< Crossed-out text
  };

  /**
   * @brief Stream insertion operator for text attributes.
   * @ingroup AlertModule
   * @param os Output stream to write to.
   * @param a The attribute to apply.
   * @return Reference to the output stream.
   *
   * Applies the specified text formatting attribute to the output stream
   * using termcolor library functionality.
   */
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

  /**
   * @brief Templated text class with foreground and background color support.
   * @ingroup AlertModule
   * @tparam F Foreground color type (defaults to NoColorT).
   * @tparam B Background color type (defaults to NoColorT).
   *
   * A flexible text formatting class that supports foreground colors, background
   * colors, and various text attributes for rich terminal output. The class uses
   * template parameters to specify color types at compile time for type safety.
   *
   * Example usage:
   * @code
   * auto redText = Text<RedT>("Error message");
   * auto boldText = Text("Important").setBold();
   * auto styledText = Text<RedT, BlueT>(RedT{}, BlueT{}, "Styled text");
   * @endcode
   */
  template <class F = NoColorT, class B = NoColorT>
  class Text
  {
    public:
      /// @brief Foreground color type alias.
      using Foreground = F;

      /// @brief Background color type alias.
      using Background = B;

      /**
       * @brief Constructs text from a C-style string.
       * @param text The text content as a null-terminated string.
       */
      Text(const char* text)
        : Text(std::string(text))
      {}

      /**
       * @brief Constructs text with foreground color from a C-style string.
       * @param fg The foreground color.
       * @param text The text content as a null-terminated string.
       */
      Text(const Foreground& fg, const char* text)
        : Text(fg, std::string(text))
      {}

      /**
       * @brief Constructs text with foreground and background colors from a C-style string.
       * @param fg The foreground color.
       * @param bg The background color.
       * @param text The text content as a null-terminated string.
       */
      Text(const Foreground& fg, const Background& bg, const char* text)
        : Text(fg, bg, std::string(text))
      {}

      /**
       * @brief Constructs text from a string.
       * @param text The text content as a std::string.
       */
      Text(const std::string& text)
        : m_string(text)
      {}

      /**
       * @brief Constructs text with foreground color from a string.
       * @param fg The foreground color.
       * @param text The text content as a std::string.
       */
      Text(const Foreground& fg, const std::string& text)
        : m_fg(fg), m_string(text)
      {}

      /**
       * @brief Constructs text with foreground and background colors from a string.
       * @param fg The foreground color.
       * @param bg The background color.
       * @param text The text content as a std::string.
       */
      Text(const Foreground& fg, const Background& bg, const std::string& text)
        : m_fg(fg), m_bg(bg), m_string(text)
      {}

      /**
       * @brief Copy constructor.
       */
      constexpr
      Text(const Text&) = default;

      /**
       * @brief Move constructor.
       */
      constexpr
      Text(Text&&) = default;

      /**
       * @brief Gets the text string content.
       * @return Const reference to the text string.
       */
      inline
      const std::string& getString() const
      {
        return m_string;
      }

      /**
       * @brief Gets the foreground color.
       * @return Const reference to the foreground color.
       */
      inline
      const Foreground& getForeground() const
      {
        return m_fg;
      }

      /**
       * @brief Gets the background color.
       * @return Const reference to the background color.
       */
      inline
      const Background& getBackground() const
      {
        return m_bg;
      }

      /**
       * @brief Sets a text formatting attribute.
       * @param a The attribute to set.
       * @return Reference to this Text object for method chaining.
       */
      Text& setAttribute(Attribute a)
      {
        m_attributes.insert(a);
        return *this;
      }

      /**
       * @brief Sets bold text formatting.
       * @return Reference to this Text object for method chaining.
       */
      Text& setBold()
      {
        setAttribute(Attribute::Bold);
        return *this;
      }

      /**
       * @brief Sets dark text formatting.
       * @return Reference to this Text object for method chaining.
       */
      Text& setDark()
      {
        setAttribute(Attribute::Bold);
        return *this;
      }

      /**
       * @brief Sets italic text formatting.
       * @return Reference to this Text object for method chaining.
       */
      Text& setItalic()
      {
        setAttribute(Attribute::Italic);
        return *this;
      }

      /**
       * @brief Sets underlined text formatting.
       * @return Reference to this Text object for method chaining.
       */
      Text& setUnderline()
      {
        setAttribute(Attribute::Underline);
        return *this;
      }

      /**
       * @brief Sets blinking text effect.
       * @return Reference to this Text object for method chaining.
       */
      Text& setBlink()
      {
        setAttribute(Attribute::Blink);
        return *this;
      }

      /**
       * @brief Sets reverse video formatting (inverted colors).
       * @return Reference to this Text object for method chaining.
       */
      Text& setReverse()
      {
        setAttribute(Attribute::Reverse);
        return *this;
      }

      /**
       * @brief Sets concealed/hidden text formatting.
       * @return Reference to this Text object for method chaining.
       */
      Text& setConcealed()
      {
        setAttribute(Attribute::Concealed);
        return *this;
      }

      /**
       * @brief Sets crossed-out text formatting.
       * @return Reference to this Text object for method chaining.
       */
      Text& setCrossed()
      {
        setAttribute(Attribute::Crossed);
        return *this;
      }

      /**
       * @brief Gets the set of applied text attributes.
       * @return Const reference to the set of attributes.
       */
      const FlatSet<Attribute>& getAttributes() const
      {
        return m_attributes;
      }

    private:
      Foreground m_fg;
      Background m_bg;
      std::string m_string;
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

