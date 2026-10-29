/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ATRIUM_LOG_H
#define RODIN_ATRIUM_LOG_H

#include <chrono>
#include <exception>
#include <iomanip>
#include <Rodin/Alert.h>
#include <Rodin/Types.h>

namespace Rodin::Examples::Heart::Atrium
{
  class Stage;

  /**
   * @brief Component-labelled progress message with message-local numeric formatting.
   *
   * Callers emit messages on the communicator root after collective diagnostics
   * have been evaluated. Slash-separated labels identify solver levels, such as
   * Solid/Newton and Fluid/KSP. Native PETSc monitors follow their solve heading.
   */
  class Log : public Alert::Message<Alert::MessagePrefix<Alert::BlueT>>
  {
    public:
      using Prefix = Alert::MessagePrefix<Alert::BlueT>;
      using Parent = Alert::Message<Prefix>;

      /**
       * @brief Starts a component-labelled message with six decimal places.
       * @param component Component or solver path displayed in the prefix.
       * @param stage Whether to draw a stage separator.
       */
      explicit Log(const std::string& component, bool stage = false)
        : Parent(std::cout,
            Prefix(std::string(2 * s_depth, ' ') + (stage ? "|--------- " : "") + "[" +
              component + "]"))
      {
        *this << std::scientific << std::setprecision(6);
      }

    private:
      friend class Stage;
      inline static thread_local size_t s_depth = 0;
  };

  /**
   * @brief Nested stage with inclusive root-rank elapsed wall time.
   *
   * ## Architecture
   *
   * A thread-local nesting level controls message indentation.
   * The stage prints its beginning and end at the same indentation level.
   * Enclosed messages and child stages are indented. Timing is local and adds
   * no MPI synchronization; parents include child work and communication waits.
   * Stack unwinding closes failed stages without changing the numerical state.
   */
  class Stage
  {
    public:
      /**
       * @brief Starts a timed stage and increases the message indentation.
       * @param name Stage name shown on both separator lines.
       * @param root Whether this rank emits the stage messages.
       */
      Stage(const std::string& name, bool root)
        : m_name(name),
          m_root(root),
          m_start(std::chrono::steady_clock::now()),
          m_exceptions(std::uncaught_exceptions())
      {
        if (m_root)
          Log(m_name, true) << "begin" << Alert::Raise;
        ++Log::s_depth;
      }

      /// @brief Closes the stage and restores its enclosing indentation.
      ~Stage()
      {
        --Log::s_depth;
        if (m_root)
        {
          Log(m_name, true) << (std::uncaught_exceptions() > m_exceptions ? "failed"
                                                                          : "end")
                            << " wall_time="
                            << std::chrono::duration<Real>(
                                 std::chrono::steady_clock::now() - m_start)
                                 .count()
                            << " s (root)" << Alert::Raise;
        }
      }

      Stage(const Stage&) = delete;
      Stage& operator=(const Stage&) = delete;

    private:
      std::string m_name;
      bool m_root;
      std::chrono::steady_clock::time_point m_start;
      int m_exceptions;
  };
}

#endif
