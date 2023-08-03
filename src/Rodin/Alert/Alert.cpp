/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>

#include "Alert.h"

namespace Rodin::Alert
{
  Alert::Alert() noexcept
    : m_indent(0),
      m_newline(false)
  {}

  Alert::Alert(int indent) noexcept
    : m_indent(indent),
      m_newline(false)
  {}

  Alert::Alert(const std::string& what, int indent)
    : m_ss(what),
      m_indent(indent),
      m_newline(false)
  {}

  Alert::Alert(const Alert& other)
    : m_indent(other.m_indent),
      m_newline(other.m_newline)
  {
    m_ss << other.m_ss.rdbuf();
  }

  Alert::Alert(Alert&& other)
    : m_ss(std::move(other.m_ss)),
      m_indent(std::move(other.m_indent)),
      m_newline(std::move(other.m_newline))
  {}

  const char* Alert::what() const noexcept
  {
    m_what = m_ss.str();
    return m_what.c_str();
  }
}
