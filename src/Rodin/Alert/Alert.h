/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ALERT_ALERT_H
#define RODIN_ALERT_ALERT_H

#include <string>
#include <sstream>
#include <iomanip>
#include <type_traits>

namespace Rodin::Alert
{
  /**
   * @brief Empty class tag type to raise an Alert.
   *
   * RaiseT is an empty class tag type to specify that a derived object of
   * Alert should be raised.
   */
  struct RaiseT
  {
    explicit RaiseT() = default;
  };

  /**
   * @brief Instance of RaiseT
   *
   * Instance of the empty struct tag type RaiseT.
   */
  static constexpr RaiseT Raise;

  static constexpr char NewLine = '\n';

  /**
   * @brief Base class for objects which represents alerts.
   *
   * Represents an Alert which can be raised and will have visible effects.
   * These effects may include but are not limited to:
   *   - Outputting to stdout or stderr
   *   - Asking user for input through stdin
   *   - Writing to a log file
   *   - Aborting the program
   *   - Any combination of the above.
   *
   * Every Alert has a message containing its description or reason for being
   * raised, which may be obtained by the @ref what() method.
   * Furthermore, an Alert may be raised via its @ref raise() method.
   */
  class Alert
  {
    public:
      /**
       * @brief Initializes an Alert with an empty message.
       */
      Alert() noexcept;

      /**
       * @brief Initializes an Alert with the given message.
       * @param[in] what Description or reason for the Alert being raised.
       */
      Alert(const std::string& what);

      /**
       * @brief Performs a copy of the Alert's message.
       * @param[in] other Object to copy.
       */
      Alert(const Alert& other);

      /**
       * @brief Performs a move of the Alert's message.
       * @param[in] other Object to move.
       */
      Alert(Alert&& other);

      /**
       * @brief Virtual destructor.
       */
      virtual ~Alert() = default;

      /**
       * @brief Gets the description (or reason) for the alert.
       * @returns String containing the message.
       */
      const char* what() const noexcept;

      /**
       * @brief Operator overload to aid in the construction of Alert
       * messages.
       */
      template <class T>
      std::enable_if_t<!std::is_same_v<RaiseT, T>, Alert&>
      operator<<(T&& v) noexcept
      {
        std::stringstream ss;
        if (m_what.back() == NewLine)
          ss << std::string(m_indent, ' ');
        ss << std::forward<T>(v);
        m_what += ss.str();
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
      virtual void raise() const = 0;

    protected:
      Alert(int setw) noexcept;

    private:
      std::string m_what;
      size_t m_indent;
  };
}

#endif
