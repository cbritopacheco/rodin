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
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <Rodin/Alert.h>
#include <Rodin/Types.h>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

namespace Rodin::Examples::Heart::Atrium
{
  class Stage;

  /**
   * @brief Component-labelled progress message with message-local numeric formatting.
   *
   * Callers emit messages on the communicator root after collective diagnostics
   * have been evaluated. Slash-separated labels identify solver levels, such as
   * Solid/Newton and Fluid/KSP. Native PETSc monitors follow their solve heading.
   * Component and status accents affect only decorative markers; labels and values
   * retain the default foreground and remain meaningful without color.
   */
  class Log : public Alert::Message<Alert::MessagePrefix<Alert::NoColorT>>
  {
    public:
      using Prefix = Alert::MessagePrefix<Alert::NoColorT>;
      using Parent = Alert::Message<Prefix>;

      /// @brief Text-labelled solver outcome, independent of terminal colors.
      enum class Status
      {
        Normal,
        Success,
        Retry,
        Failure
      };

      /**
       * @brief Starts a component-labelled message with six decimal places.
       * @param component Component or solver path displayed in the prefix.
       * @param stage Whether to draw a stage separator.
       * @param status Solver outcome used to color only the decorative marker.
       */
      explicit Log(
        const std::string& component, bool stage = false, Status status = Status::Normal)
        : Parent(std::cout,
            Prefix(std::string(2 * s_depth, ' ') + (stage ? "|--------- " : "| ") + "[" +
              component + "]")),
          m_indent(2 * s_depth, ' '),
          m_marker(stage ? "|---------" : "|"),
          m_label("[" + component + "]"),
          m_color(color(component, status))
      {
        *this << std::scientific << std::setprecision(6);
      }

      /**
       * @brief Prints readable labels and values in the terminal's default foreground.
       *
       * Only the decorative marker is colored. Non-terminal output, NO_COLOR,
       * CLICOLOR=0, and TERM=dumb select plain text. CLICOLOR_FORCE overrides
       * terminal detection (useful for MPI launchers), but never NO_COLOR.
       */
      void raise() const override
      {
        const char* noColor = std::getenv("NO_COLOR");
        const char* force = std::getenv("CLICOLOR_FORCE");
        const char* clicolor = std::getenv("CLICOLOR");
        const char* term = std::getenv("TERM");
#ifdef _WIN32
        const bool terminal = _isatty(_fileno(stdout));
#else
        const bool terminal = isatty(fileno(stdout));
#endif
        const bool forced = force && *force && std::string(force) != "0";
        const bool plain = (noColor && *noColor) ||
          (clicolor && std::string(clicolor) == "0") ||
          (term && std::string(term) == "dumb");
        const bool styled = !plain && (terminal || forced);
        std::cout << m_indent;
        if (styled && m_color)
        {
          std::cout << Alert::Stylize
                    << Alert::Text<Alert::Color16>(*m_color, m_marker).setBold();
        }
        else
          std::cout << m_marker;
        std::cout << ' ';
        if (styled)
          std::cout << Alert::Text<>(m_label).setBold();
        else
          std::cout << m_label;
        if (m_label.starts_with("[Timestep "))
        {
          const std::string text(what());
          const auto split = text.find(' ');
          constexpr size_t BannerWidth = 100;
          const std::string event = text.substr(0, split);
          const size_t used =
            m_indent.size() + m_marker.size() + 1 + m_label.size() + 2 + event.size() + 3;
          const size_t dashes = BannerWidth > used ? BannerWidth - used : 1;
          std::cout << ": " << event << ' ' << std::string(dashes, '-') << " |";
          if (split != std::string::npos)
            std::cout << text.substr(split);
          std::cout << std::endl;
        }
        else
          std::cout << ": " << what() << std::endl;
      }

    private:
      /**
       * @brief Selects a standard terminal-palette accent without prescribing a background.
       * @param component Component label or stage name.
       * @param status Solver outcome, taking precedence over component accents.
       * @returns Accent color, or no accent for run headings and diagnostics.
       */
      static Optional<Alert::Color16> color(StringView component, Status status)
      {
        if (status == Status::Failure)
          return Alert::Color16::Red;
        if (status == Status::Retry || component.find("retry") != StringView::npos)
          return Alert::Color16::Yellow;
        if (status == Status::Success)
          return Alert::Color16::Green;
        if (component.find("diagnostics") != StringView::npos ||
          component.find("export") != StringView::npos ||
          component.ends_with("/Timing") || component.ends_with("/State") ||
          component.ends_with("/Cache"))
          return {};
        if (component.starts_with("ALE"))
          return Alert::Color16::Yellow;
        if (component.starts_with("Solid") || component.starts_with("Passive solid") ||
          component == "Construct solid")
          return Alert::Color16::Magenta;
        if (component.starts_with("Fluid") || component == "Construct fluid" ||
          component.starts_with("Wall shear") || component.starts_with("Species") ||
          component.starts_with("VMS") || component.starts_with("Ventricular"))
          return Alert::Color16::Blue;
        if (component.starts_with("FSI") || component.starts_with("Coupling") ||
          component.starts_with("Interface") ||
          component.starts_with("Accepted interface"))
          return Alert::Color16::Cyan;
        return {};
      }

      std::string m_indent, m_marker, m_label;
      Optional<Alert::Color16> m_color;
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
          Log(m_name, true,
            std::uncaught_exceptions() > m_exceptions ? Log::Status::Failure
                                                      : Log::Status::Normal)
            << (std::uncaught_exceptions() > m_exceptions ? "failed" : "end")
            << " wall_time="
            << std::chrono::duration<Real>(std::chrono::steady_clock::now() - m_start)
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
