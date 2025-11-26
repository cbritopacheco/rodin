/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_MESSAGE_H
#define RODIN_ALERT_MESSAGE_H

#include <string>
#include <sstream>
#include <cstdint>
#include <cstdlib>
#include <iostream>

#include <iomanip>
#include <type_traits>

#include "ForwardDecls.h"
#include "Text.h"
#include "Color.h"
#include "Raise.h"
#include "NewLine.h"
#include "Notation.h"
#include "Stylize.h"

namespace Rodin::Alert
{
  /**
   * @internal
   * @brief Internal helper traits for detecting streamable types.
   *
   * Namespace containing implementation details for the Alert module.
   */
  namespace Internal
  {
    /**
     * @brief Type trait to check if a type can be output to a stream.
     * @tparam T The type to check.
     *
     * Uses SFINAE to detect if a type T can be used with operator<<
     * on an output stream. The Value member is true if T is streamable,
     * false otherwise.
     */
    template <typename T>
    class CanBeOutput
    {
      template <
        class U,
        class = decltype(std::declval<std::ostream&>() << std::declval<const U&>())>
      static std::true_type test(U*);

      template <typename>
      static std::false_type test(...);

    public:
        /// @brief True if T can be streamed to an ostream, false otherwise.
        static constexpr bool Value = decltype(test<T>(nullptr))::value;
    };
  }

  /**
   * @brief Base class for message prefixes with colored text.
   * @ingroup AlertModule
   * @tparam Foreground The foreground color type.
   *
   * Template class that provides colored, bolded prefix text for messages.
   * This class is used as a base for specific message type prefixes like
   * ExceptionPrefix, WarningPrefix, etc.
   */
  template <class Foreground>
  class MessagePrefix : public Text<Foreground>
  {
    public:
      /// @brief Parent class type alias.
      using Parent = Text<Foreground>;

      /**
       * @brief Constructs a message prefix with foreground color and text.
       * @param fg The foreground color.
       * @param prefix The prefix text string.
       */
      MessagePrefix(const Foreground& fg, const std::string& prefix)
        : Parent(fg, prefix)
      {
        this->setBold();
      }

      /**
       * @brief Constructs a message prefix with text using default color.
       * @param prefix The prefix text string.
       */
      MessagePrefix(const std::string& prefix)
        : Parent(prefix)
      {
        this->setBold();
      }

      /**
       * @brief Copy constructor.
       */
      MessagePrefix(const MessagePrefix& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       */
      MessagePrefix(MessagePrefix&& other)
        : Parent(std::move(other))
      {}
  };

  /**
   * @brief Base template class for formatted message output.
   * @ingroup AlertModule
   * @tparam Prefix The prefix type for the message.
   *
   * Base class for creating formatted messages with customizable prefixes,
   * colors, and styles. Messages support stream-like insertion operators
   * for building content incrementally. This class is used as the base
   * for Exception, Warning, Info, and Success message types.
   *
   * The Message class handles:
   * - Formatted output with colored prefixes
   * - Automatic indentation for multi-line messages
   * - Stream-like insertion operators for content
   * - Integration with the Alert tag types (NewLine, Raise, etc.)
   *
   * Example usage (through derived classes):
   * @code{.cpp}
   * Exception() << "Error in function: " << funcName << NewLine
   *             << "Invalid parameter value: " << value << Raise;
   * @endcode
   */
  template <class Prefix>
  class Message
  {
    public:
      /**
       * @brief Constructs a message with output stream and prefix.
       * @param os The output stream to write the message to.
       * @param prefix The message prefix object.
       */
      Message(std::ostream& os, const Prefix& prefix) noexcept
        : m_os(os),
          m_prefix(prefix),
          m_newline(false)
      {
        m_styled << Stylize << prefix << ": ";
      }

      /**
       * @brief Copy constructor.
       * @param other Object to copy.
       *
       * Performs a copy of the message's state including the accumulated
       * message content.
       */
      Message(const Message& other)
        : m_os(other.m_os),
          m_prefix(other.m_prefix),
          m_newline(other.m_newline)
      {
        m_ss << other.m_ss.rdbuf();
        m_styled << other.m_ss.rdbuf();
      }

      /**
       * @brief Move constructor.
       * @param other Object to move.
       */
      Message(Message&& other) = default;

      /**
       * @brief Virtual destructor.
       */
      virtual ~Message() = default;

      /**
       * @brief Gets the message content as a C-string.
       * @return Null-terminated string containing the message.
       *
       * Returns the accumulated message content without formatting or
       * color codes. Useful for logging or storing the message text.
       */
      const char* what() const noexcept
      {
        s_what = m_ss.str();
        return s_what.c_str();
      }

      /**
       * @brief Stream insertion operator for arbitrary streamable types.
       * @tparam T The type to insert (must be streamable).
       * @param v The value to insert into the message.
       * @return Reference to this Message object for method chaining.
       *
       * Appends content to the message. Handles automatic indentation
       * when content follows a newline. Only enabled for types that
       * can be output to an ostream.
       */
      template <class T>
      std::enable_if_t<Internal::CanBeOutput<T>::Value, Message&>
      operator<<(const T& v) noexcept
      {
        if (m_newline)
        {
          const auto indent = std::string(m_prefix.getString().size() + 2, ' ');
          m_ss << indent;
          m_styled << indent;
        }
        m_ss << v;
        m_styled << v;
        m_newline = false;
        return *this;
      }

      /**
       * @brief Stream insertion operator for NewLine tag.
       * @param The NewLine tag (unused).
       * @return Reference to this Message object for method chaining.
       *
       * Inserts a newline character and marks the next insertion for
       * automatic indentation to align with the message prefix.
       */
      Message& operator<<(const NewLineT&)
      {
        operator<<('\n');
        m_newline = true;
        return *this;
      }

      /**
       * @brief Stream insertion operator for Raise tag.
       * @param The Raise tag (unused).
       *
       * Triggers the raise() method to output the message and perform
       * any associated actions (such as program termination for exceptions).
       */
      void operator<<(const RaiseT&)
      {
        this->raise();
      }

      /**
       * @brief Raises (outputs) the message to the user.
       *
       * Default behavior outputs the formatted message to the configured
       * output stream. Derived classes may override this to add additional
       * behavior (e.g., Exception terminates the program).
       */
      virtual void raise() const
      {
        m_os.get() << m_styled.rdbuf() << NewLine;
      }

      /**
       * @brief Sets the output stream for this message.
       * @param os The new output stream.
       *
       * Changes where the message will be output when raised.
       */
      void setOutputStream(std::ostream& os)
      {
        m_os = os;
      }

    private:
      /// @brief Thread-local storage for the what() string.
      static thread_local std::string s_what;

      std::reference_wrapper<std::ostream> m_os;
      Prefix m_prefix;
      std::stringstream m_ss;
      std::stringstream m_styled;
      bool m_newline;
  };

  /// @brief Thread-local storage initialization for Message::s_what.
  template <class Prefix>
  thread_local std::string Message<Prefix>::s_what;
}

#endif
