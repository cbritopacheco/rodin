/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_COLOR_H
#define RODIN_ALERT_COLOR_H

#include <ostream>
#include <termcolor/termcolor.hpp>

namespace Rodin::Alert
{
  enum class Color16
  {
    Red,
    Green,
    Blue,
    Yellow,
    Magenta,
    Cyan,
    Gray,
    White,
    BrightGray,
    BrightWhite
  };

  inline
  std::ostream& operator<<(std::ostream& os, Color16 c)
  {
    switch (c)
    {
      case Color16::Red:
      {
        os << termcolor::red;
        break;
      }
      case Color16::Green:
      {
        os << termcolor::green;
        break;
      }
      case Color16::Blue:
      {
        os << termcolor::blue;
        break;
      }
      case Color16::Yellow:
      {
        os << termcolor::yellow;
        break;
      }
      case Color16::Magenta:
      {
        os << termcolor::magenta;
        break;
      }
      case Color16::Cyan:
      {
        os << termcolor::cyan;
        break;
      }
      case Color16::White:
      {
        os << termcolor::white;
        break;
      }
      case Color16::Gray:
      {
        os << termcolor::grey;
        break;
      }
      case Color16::BrightGray:
      {
        os << termcolor::bright_grey;
        break;
      }
      case Color16::BrightWhite:
      {
        os << termcolor::bright_white;
        break;
      }
    }
    return os;
  }

  template <size_t RED, size_t GREEN, size_t BLUE>
  struct RGB
  {
    static_assert(RED < 256);
    static_assert(GREEN < 256);
    static_assert(BLUE < 256);

    static constexpr size_t R = RED;
    static constexpr size_t G = GREEN;
    static constexpr size_t B = BLUE;

    static constexpr size_t HEX = (R << 16) + (G << 8) + B;

    constexpr
    RGB() = default;

    constexpr
    RGB(const RGB&) = default;

    constexpr
    RGB(RGB&&) = default;

    constexpr
    RGB& operator=(const RGB&) = default;

    constexpr
    RGB& operator=(RGB&&) = default;

    inline
    constexpr
    size_t r() const
    {
      return R;
    }

    inline
    constexpr
    size_t g() const
    {
      return G;
    }

    inline
    constexpr
    size_t b() const
    {
      return B;
    }

    inline
    constexpr
    size_t hex() const
    {
      return HEX;
    }
  };

  struct NoColorT
  {
    constexpr
    NoColorT() = default;

    constexpr
    NoColorT(const NoColorT&) = default;

    constexpr
    NoColorT(NoColorT&&) = default;
  };

  static constexpr NoColorT NoColor;

  template <class CodeT>
  class Color
  {
    public:
      using Code = CodeT;

      constexpr
      Color(const Code& code = Code())
        : m_code(code)
      {}

      constexpr
      Color(const Color&) = default;

      constexpr
      Color(Color&&) = default;

      constexpr
      Color& operator=(const Color&) = default;

      constexpr
      Color& operator=(Color&&) = default;

      inline
      constexpr
      const Code& getCode() const
      {
        return m_code;
      }

    private:
      Code m_code;
  };

  template <class Code>
  Color(const Code&) -> Color<Code>;

  template <class Code>
  inline
  std::ostream& operator<<(std::ostream& os, const Color<Code>&)
  {
    os << termcolor::color<Code::R, Code::G, Code::B>;
    return os;
  }

  struct RedT {};

  static constexpr RedT Red;

  inline
  std::ostream& operator<<(std::ostream& os, const RedT&)
  {
    os << Color16::Red;
    return os;
  }

  struct GreenT {};

  static constexpr GreenT Green;

  inline
  std::ostream& operator<<(std::ostream& os, const GreenT&)
  {
    os << Color16::Green;
    return os;
  }

  struct BlueT {};

  static constexpr BlueT Blue;

  inline
  std::ostream& operator<<(std::ostream& os, const BlueT&)
  {
    os << Color16::Blue;
    return os;
  }

  struct YellowT {};

  static constexpr YellowT Yellow;

  inline
  std::ostream& operator<<(std::ostream& os, const YellowT&)
  {
    os << Color16::Yellow;
    return os;
  }

  struct MagentaT {};

  static constexpr MagentaT Magenta;

  inline
  std::ostream& operator<<(std::ostream& os, const MagentaT&)
  {
    os << Color16::Magenta;
    return os;
  }

  struct CyanT {};

  static constexpr CyanT Cyan;

  inline
  std::ostream& operator<<(std::ostream& os, const CyanT&)
  {
    os << Color16::Cyan;
    return os;
  }

  struct WhiteT {};

  static constexpr WhiteT White;

  inline
  std::ostream& operator<<(std::ostream& os, const WhiteT&)
  {
    os << Color16::White;
    return os;
  }

  struct GrayT {};

  static constexpr GrayT Gray;

  inline
  std::ostream& operator<<(std::ostream& os, const GrayT&)
  {
    os << Color16::Gray;
    return os;
  }

  struct BrightGrayT {};

  static constexpr BrightGrayT BrightGray;

  inline
  std::ostream& operator<<(std::ostream& os, const BrightGrayT&)
  {
    os << Color16::BrightGray;
    return os;
  }

  struct BrightWhiteT {};

  static constexpr BrightWhiteT BrightWhite;

  inline
  std::ostream& operator<<(std::ostream& os, const BrightWhiteT&)
  {
    os << Color16::BrightWhite;
    return os;
  }
}

#endif

