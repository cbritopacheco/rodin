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
  /**
   * @brief Enumeration of basic 16-color terminal colors.
   * @ingroup AlertModule
   *
   * Provides an enumeration of the standard 16-color terminal palette
   * including primary colors, gray shades, and bright variants. These
   * colors are widely supported across terminal emulators.
   */
  enum class Color16
  {
    Red,          ///< Standard red color
    Green,        ///< Standard green color
    Blue,         ///< Standard blue color
    Yellow,       ///< Standard yellow color
    Magenta,      ///< Standard magenta color
    Cyan,         ///< Standard cyan color
    Gray,         ///< Standard gray color
    White,        ///< Standard white color
    BrightGray,   ///< Bright gray color
    BrightWhite   ///< Bright white color
  };

  /**
   * @brief Stream insertion operator for Color16.
   * @ingroup AlertModule
   * @param os Output stream to write to.
   * @param c The color to apply.
   * @return Reference to the output stream.
   *
   * Applies the specified 16-color terminal color to the output stream
   * using the termcolor library.
   */
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

  /**
   * @brief RGB color template with compile-time color values.
   * @ingroup AlertModule
   * @tparam RED Red component value (0-255).
   * @tparam GREEN Green component value (0-255).
   * @tparam BLUE Blue component value (0-255).
   *
   * Template struct for representing RGB colors with compile-time validation.
   * The color components are constrained to the range [0, 255] using static
   * assertions. This enables type-safe custom terminal colors.
   *
   * Example usage:
   * @code{.cpp}
   * using MyColor = RGB<128, 64, 255>;
   * @endcode
   */
  template <size_t RED, size_t GREEN, size_t BLUE>
  struct RGB
  {
    static_assert(RED < 256);
    static_assert(GREEN < 256);
    static_assert(BLUE < 256);

    static constexpr size_t R = RED;      ///< Red component value
    static constexpr size_t G = GREEN;    ///< Green component value
    static constexpr size_t B = BLUE;     ///< Blue component value

    /// @brief Combined hexadecimal representation of the RGB color.
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

    /**
     * @brief Gets the red component value.
     * @return The red component (0-255).
     */
    inline
    constexpr
    size_t r() const
    {
      return R;
    }

    /**
     * @brief Gets the green component value.
     * @return The green component (0-255).
     */
    inline
    constexpr
    size_t g() const
    {
      return G;
    }

    /**
     * @brief Gets the blue component value.
     * @return The blue component (0-255).
     */
    inline
    constexpr
    size_t b() const
    {
      return B;
    }

    /**
     * @brief Gets the hexadecimal representation of the color.
     * @return The color as a hexadecimal value.
     */
    inline
    constexpr
    size_t hex() const
    {
      return HEX;
    }
  };

  /**
   * @brief Tag type representing no color (transparent/default).
   * @ingroup AlertModule
   *
   * Empty tag type used to indicate the absence of a foreground or
   * background color. This allows the Text class to use default terminal
   * colors when no specific color is desired.
   */
  struct NoColorT
  {
    constexpr
    NoColorT() = default;

    constexpr
    NoColorT(const NoColorT&) = default;

    constexpr
    NoColorT(NoColorT&&) = default;
  };

  /**
   * @brief Instance of NoColorT tag type.
   * @ingroup AlertModule
   *
   * Constant instance of the NoColorT tag type for convenient usage.
   */
  static constexpr NoColorT NoColor;

  /**
   * @brief Type-safe wrapper for custom RGB colors.
   * @ingroup AlertModule
   * @tparam CodeT The RGB color code type.
   *
   * Template class providing a type-safe wrapper around RGB color codes.
   * Used to create custom terminal colors with compile-time color validation.
   */
  template <class CodeT>
  class Color
  {
    public:
      /// @brief The RGB color code type.
      using Code = CodeT;

      /**
       * @brief Constructs a Color with the given code.
       * @param code The RGB color code (defaults to default-constructed Code).
       */
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

      /**
       * @brief Gets the RGB color code.
       * @return Const reference to the color code.
       */
      inline
      constexpr
      const Code& getCode() const
      {
        return m_code;
      }

    private:
      Code m_code;
  };

  /**
   * @brief Deduction guide for Color.
   * @tparam Code The color code type.
   */
  template <class Code>
  Color(const Code&) -> Color<Code>;

  /**
   * @brief Stream insertion operator for custom RGB colors.
   * @ingroup AlertModule
   * @tparam Code The RGB color code type.
   * @param os Output stream to write to.
   * @return Reference to the output stream.
   *
   * Applies a custom RGB color to the output stream using the termcolor
   * library's color template.
   */
  template <class Code>
  inline
  std::ostream& operator<<(std::ostream& os, const Color<Code>&)
  {
    os << termcolor::color<Code::R, Code::G, Code::B>;
    return os;
  }

  /**
   * @brief Tag type for red terminal color.
   * @ingroup AlertModule
   */
  struct RedT {};

  /**
   * @brief Instance of RedT tag type.
   * @ingroup AlertModule
   */
  static constexpr RedT Red;

  /**
   * @brief Stream insertion operator for red color.
   * @ingroup AlertModule
   */
  inline
  std::ostream& operator<<(std::ostream& os, const RedT&)
  {
    os << Color16::Red;
    return os;
  }

  /**
   * @brief Tag type for green terminal color.
   * @ingroup AlertModule
   */
  struct GreenT {};

  /**
   * @brief Instance of GreenT tag type.
   * @ingroup AlertModule
   */
  static constexpr GreenT Green;

  /**
   * @brief Stream insertion operator for green color.
   * @ingroup AlertModule
   */
  inline
  std::ostream& operator<<(std::ostream& os, const GreenT&)
  {
    os << Color16::Green;
    return os;
  }

  /**
   * @brief Tag type for blue terminal color.
   * @ingroup AlertModule
   */
  struct BlueT {};

  /**
   * @brief Instance of BlueT tag type.
   * @ingroup AlertModule
   */
  static constexpr BlueT Blue;

  /**
   * @brief Stream insertion operator for blue color.
   * @ingroup AlertModule
   */
  inline
  std::ostream& operator<<(std::ostream& os, const BlueT&)
  {
    os << Color16::Blue;
    return os;
  }

  /**
   * @brief Tag type for yellow terminal color.
   * @ingroup AlertModule
   */
  struct YellowT {};

  /**
   * @brief Instance of YellowT tag type.
   * @ingroup AlertModule
   */
  static constexpr YellowT Yellow;

  /**
   * @brief Stream insertion operator for yellow color.
   * @ingroup AlertModule
   */
  inline
  std::ostream& operator<<(std::ostream& os, const YellowT&)
  {
    os << Color16::Yellow;
    return os;
  }

  /**
   * @brief Tag type for magenta terminal color.
   * @ingroup AlertModule
   */
  struct MagentaT {};

  /**
   * @brief Instance of MagentaT tag type.
   * @ingroup AlertModule
   */
  static constexpr MagentaT Magenta;

  /**
   * @brief Stream insertion operator for magenta color.
   * @ingroup AlertModule
   */
  inline
  std::ostream& operator<<(std::ostream& os, const MagentaT&)
  {
    os << Color16::Magenta;
    return os;
  }

  /**
   * @brief Tag type for cyan terminal color.
   * @ingroup AlertModule
   */
  struct CyanT {};

  /**
   * @brief Instance of CyanT tag type.
   * @ingroup AlertModule
   */
  static constexpr CyanT Cyan;

  /**
   * @brief Stream insertion operator for cyan color.
   * @ingroup AlertModule
   */
  inline
  std::ostream& operator<<(std::ostream& os, const CyanT&)
  {
    os << Color16::Cyan;
    return os;
  }

  /**
   * @brief Tag type for white terminal color.
   * @ingroup AlertModule
   */
  struct WhiteT {};

  /**
   * @brief Instance of WhiteT tag type.
   * @ingroup AlertModule
   */
  static constexpr WhiteT White;

  /**
   * @brief Stream insertion operator for white color.
   * @ingroup AlertModule
   */
  inline
  std::ostream& operator<<(std::ostream& os, const WhiteT&)
  {
    os << Color16::White;
    return os;
  }

  /**
   * @brief Tag type for gray terminal color.
   * @ingroup AlertModule
   */
  struct GrayT {};

  /**
   * @brief Instance of GrayT tag type.
   * @ingroup AlertModule
   */
  static constexpr GrayT Gray;

  /**
   * @brief Stream insertion operator for gray color.
   * @ingroup AlertModule
   */
  inline
  std::ostream& operator<<(std::ostream& os, const GrayT&)
  {
    os << Color16::Gray;
    return os;
  }

  /**
   * @brief Tag type for bright gray terminal color.
   * @ingroup AlertModule
   */
  struct BrightGrayT {};

  /**
   * @brief Instance of BrightGrayT tag type.
   * @ingroup AlertModule
   */
  static constexpr BrightGrayT BrightGray;

  /**
   * @brief Stream insertion operator for bright gray color.
   * @ingroup AlertModule
   */
  inline
  std::ostream& operator<<(std::ostream& os, const BrightGrayT&)
  {
    os << Color16::BrightGray;
    return os;
  }

  /**
   * @brief Tag type for bright white terminal color.
   * @ingroup AlertModule
   */
  struct BrightWhiteT {};

  /**
   * @brief Instance of BrightWhiteT tag type.
   * @ingroup AlertModule
   */
  static constexpr BrightWhiteT BrightWhite;

  /**
   * @brief Stream insertion operator for bright white color.
   * @ingroup AlertModule
   */
  inline
  std::ostream& operator<<(std::ostream& os, const BrightWhiteT&)
  {
    os << Color16::BrightWhite;
    return os;
  }
}

#endif

