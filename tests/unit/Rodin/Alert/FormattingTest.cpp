/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <array>
#include <complex>
#include <sstream>
#include <utility>

#include <gtest/gtest.h>

#include <Rodin/Alert/Color.h>
#include <Rodin/Alert/MemberFunctionWarning.h>
#include <Rodin/Alert/Notation.h>
#include <Rodin/Alert/Reset.h>
#include <Rodin/Alert/Stylize.h>
#include <Rodin/Alert/Text.h>

namespace Rodin::Tests::Unit
{
  using namespace Alert;

  /// @brief Generic notation owns temporary strings and handles empty and populated ranges.
  TEST(AlertFormatting, GenericNotationValues)
  {
    const Notation text(std::string("owned text"));
    EXPECT_EQ(text.getString(), "owned text");
    EXPECT_EQ(Notation::Print(std::complex<int>(2, -3)).getString(), "(2,-3)");
    const std::array<int, 3> values{2, -3, 7};
    EXPECT_EQ(Notation::Set(values.begin(), values.end()).getString(), "{ 2, -3, 7 }");
    EXPECT_EQ(Notation::Set(values.end(), values.end()).getString(), "{  }");
  }

  /// @brief Member warnings retain function/class context and emit without throwing.
  TEST(AlertFormatting, MemberWarningContext)
  {
    struct WarningOrigin
    {};
    const WarningOrigin origin;
    MemberFunctionWarning warning(origin, "evaluate");
    warning << "legacy input";
    const std::string expected = "In member function evaluate of class " +
      boost::typeindex::type_id_with_cvr<WarningOrigin>().pretty_name() +
      ": legacy input";
    EXPECT_EQ(std::string(warning.what()), expected);
    std::ostringstream output, expectedOutput;
    warning.setOutputStream(output);
    expectedOutput << Stylize << WarningPrefix() << ": In member function "
                   << Identifier::Function("evaluate") << " of class "
                   << Identifier::Class(
                        boost::typeindex::type_id_with_cvr<WarningOrigin>().pretty_name())
                   << ": legacy input\n";
    EXPECT_NO_THROW(warning.raise());
    EXPECT_EQ(output.str(), expectedOutput.str());
  }

  /// @brief Every palette enum and named tag emits its prescribed ANSI code.
  TEST(AlertFormatting, PaletteAndNamedTags)
  {
    const std::array palette{Color16::Red, Color16::Green, Color16::Blue, Color16::Yellow,
      Color16::Magenta, Color16::Cyan, Color16::Gray, Color16::White, Color16::BrightGray,
      Color16::BrightWhite};
    const std::array codes{"\x1b[31m", "\x1b[32m", "\x1b[34m", "\x1b[33m", "\x1b[35m",
      "\x1b[36m", "\x1b[30m", "\x1b[37m", "\x1b[90m", "\x1b[97m"};
    for (size_t index = 0; index < palette.size(); ++index)
    {
      SCOPED_TRACE(index);
      std::ostringstream stream;
      EXPECT_EQ(&(stream << Stylize), &stream);
      EXPECT_TRUE(stream.str().empty());
      EXPECT_EQ(&(stream << palette[index]), &stream);
      EXPECT_EQ(stream.str(), codes[index]);
    }
    std::ostringstream named;
    named << Stylize << Red << Green << Blue << Yellow << Magenta << Cyan << Gray << White
          << BrightGray << BrightWhite;
    EXPECT_EQ(named.str(),
      "\x1b[31m\x1b[32m\x1b[34m\x1b[33m\x1b[35m\x1b[36m\x1b[30m\x1b[37m\x1b[90m\x1b[97m");
  }

  /// @brief RGB value semantics preserve components, packed code and output.
  TEST(AlertFormatting, RGBValuesAndLifecycle)
  {
    using Code = RGB<0, 128, 255>;
    const Code original;
    Code copied(original), moved(std::move(copied));
    copied = original;
    moved = std::move(copied);
    EXPECT_EQ(moved.r(), 0u);
    EXPECT_EQ(moved.g(), 128u);
    EXPECT_EQ(moved.b(), 255u);
    EXPECT_EQ(moved.hex(), 0x0080ffu);
    const Color originalColor(original);
    Color copiedColor(originalColor), movedColor(std::move(copiedColor));
    copiedColor = originalColor;
    movedColor = std::move(copiedColor);
    EXPECT_EQ(movedColor.getCode().hex(), original.hex());
    std::ostringstream stream;
    stream << Stylize;
    EXPECT_EQ(&(stream << movedColor), &stream);
    EXPECT_EQ(stream.str(), "\x1b[38;2;0;128;255m");
    stream << Color(RGB<0, 0, 0>{}) << Color(RGB<255, 255, 255>{});
    EXPECT_EQ(stream.str(), "\x1b[38;2;0;128;255m\x1b[38;2;0;0;0m\x1b[38;2;255;255;255m");
  }

  /// @brief Attributes and reset obey independent, explicitly enabled streams.
  TEST(AlertFormatting, AttributesResetAndStreamIsolation)
  {
    const std::array attributes{Attribute::Bold, Attribute::Dark, Attribute::Italic,
      Attribute::Underline, Attribute::Blink, Attribute::Reverse, Attribute::Concealed,
      Attribute::Crossed};
    const std::array codes{"\x1b[1m", "\x1b[2m", "\x1b[3m", "\x1b[4m", "\x1b[5m",
      "\x1b[7m", "\x1b[8m", "\x1b[9m"};
    std::ostringstream plain, styled;
    styled << Stylize;
    for (size_t index = 0; index < attributes.size(); ++index)
    {
      SCOPED_TRACE(index);
      styled.str("");
      EXPECT_EQ(&(styled << attributes[index]), &styled);
      EXPECT_EQ(styled.str(), codes[index]);
      EXPECT_EQ(&(plain << attributes[index]), &plain);
    }
    plain << Red << Color(RGB<0, 128, 255>{}) << Reset << "plain";
    EXPECT_EQ(plain.str(), "plain");
    styled.str("");
    EXPECT_EQ(&(styled << Reset), &styled);
    styled << Red << "styled" << Reset;
    EXPECT_EQ(styled.str(), "\x1b[00m\x1b[31mstyled\x1b[00m");
    plain << Stylize << Green;
    EXPECT_EQ(plain.str(), "plain\x1b[32m");
    styled << termcolor::nocolorize << Blue << Reset << " tail";
    EXPECT_EQ(styled.str(), "\x1b[00m\x1b[31mstyled\x1b[00m tail");
  }

  /// @brief Each named text setter selects exactly its documented attribute.
  TEST(AlertFormatting, TextStyleSetters)
  {
    using TextType = Text<>;
    using Setter = TextType& (TextType::*)();
    const std::array<Setter, 8> setters{&TextType::setBold, &TextType::setDark,
      &TextType::setItalic, &TextType::setUnderline, &TextType::setBlink,
      &TextType::setReverse, &TextType::setConcealed, &TextType::setCrossed};
    const std::array attributes{Attribute::Bold, Attribute::Dark, Attribute::Italic,
      Attribute::Underline, Attribute::Blink, Attribute::Reverse, Attribute::Concealed,
      Attribute::Crossed};
    const std::array codes{"\x1b[1m", "\x1b[2m", "\x1b[3m", "\x1b[4m", "\x1b[5m",
      "\x1b[7m", "\x1b[8m", "\x1b[9m"};
    for (size_t index = 0; index < setters.size(); ++index)
    {
      SCOPED_TRACE(index);
      TextType text("message");
      EXPECT_EQ(&(text.*setters[index])(), &text);
      EXPECT_EQ(&(text.*setters[index])(), &text);
      ASSERT_EQ(text.getAttributes().size(), 1u);
      EXPECT_EQ(*text.getAttributes().begin(), attributes[index]);
      TextType copied(text);
      TextType moved(std::move(copied));
      EXPECT_EQ(moved.getString(), "message");
      std::ostringstream stream;
      stream << Stylize;
      EXPECT_EQ(&(stream << moved), &stream);
      EXPECT_EQ(stream.str(), std::string("\x1b[00m") + codes[index] + "message\x1b[00m");
    }
  }
}
