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
    : m_indent(0)
  {}

  Alert::Alert(int setw) noexcept
    : m_indent(setw)
  {}

  Alert::Alert(const std::string& what)
    : m_what(what)
  {}

  Alert::Alert(const Alert& other)
    : m_what(other.m_what),
      m_indent(other.m_indent)
  {}

  Alert::Alert(Alert&& other)
    : m_what(std::move(other.m_what)),
      m_indent(std::move(other.m_indent))
  {}

  const char* Alert::what() const noexcept
  {
    return m_what.c_str();
  }
}
