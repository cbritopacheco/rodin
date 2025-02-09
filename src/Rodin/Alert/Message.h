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
  namespace Internal
  {
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
        static constexpr bool Value = decltype(test<T>(nullptr))::value;
    };
  }

  template <class Foreground>
  class MessagePrefix : public Text<Foreground>
  {
    public:
      using Parent = Text<Foreground>;

      MessagePrefix(const Foreground& fg, const std::string& prefix)
        : Parent(fg, prefix)
      {
        this->setBold();
      }

      MessagePrefix(const std::string& prefix)
        : Parent(prefix)
      {
        this->setBold();
      }

      MessagePrefix(const MessagePrefix& other)
        : Parent(other)
      {}

      MessagePrefix(MessagePrefix&& other)
        : Parent(std::move(other))
      {}
  };

  /**
   * @brief Base class for objects which represents output messages.
   *
   * Represents a message to output to the user with potential visible effects.
   */
  template <class Prefix>
  class Message
  {
    public:
      Message(std::ostream& os, const Prefix& prefix) noexcept
        : m_os(os),
          m_prefix(prefix),
          m_newline(false)
      {
        m_styled << Stylize << prefix << ": ";
      }

      /**
       * @brief Performs a copy of the Alert's message.
       * @param[in] other Object to copy.
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
       * @brief Performs a move of the Alert's message.
       * @param[in] other Object to move.
       */
      Message(Message&& other) = default;

      /**
       * @brief Virtual destructor.
       */
      virtual ~Message() = default;

      /**
       * @brief Gets the description (or reason) for the alert.
       * @returns String containing the message.
       */
      const char* what() const noexcept
      {
        s_what = m_ss.str();
        return s_what.c_str();
      }

      /**
       * @brief Operator overload to aid in the construction of Alert
       * messages.
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

      Message& operator<<(const NewLineT&)
      {
        operator<<('\n');
        m_newline = true;
        return *this;
      }

      /**
       * @brief Operator overload to raise the Alert from a stream.
       *
       * This method will call @ref raise().
       */
      void operator<<(const RaiseT&)
      {
        this->raise();
      }

      /**
       * @brief Raises the Alert to the user.
       *
       * The actual behaviour for raising the Alert is specified in its
       * subclasses by overriding this function.
       */
      virtual void raise() const
      {
        m_os.get() << m_styled.rdbuf() << NewLine;
      }

      void setOutputStream(std::ostream& os)
      {
        m_os = os;
      }

    private:
      static thread_local std::string s_what;

      std::reference_wrapper<std::ostream> m_os;
      Prefix m_prefix;
      std::stringstream m_ss;
      std::stringstream m_styled;
      bool m_newline;
  };

  template <class Prefix>
  thread_local std::string Message<Prefix>::s_what;
}

#endif
