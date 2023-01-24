/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include "Exceptions.h"

#include "Composition.h"

namespace Rodin::Variational
{
  Composition<std::function<double(double)>, FunctionBase>
  ::Composition(std::function<double(double)> f, const FunctionBase& g)
    : m_f(f), m_g(g.copy())
  {
    if (g.getRangeType() != RangeType::Scalar)
    {
      UnexpectedRangeTypeException(RangeType::Scalar, g.getRangeType()).raise();
    }
  }

  Composition<std::function<double(double)>, FunctionBase>
  ::Composition(const Composition& other)
    :  ScalarFunctionBase(other),
      m_f(other.m_f), m_g(other.m_g->copy())
  {}

  Composition<std::function<double(double)>, FunctionBase>
  ::Composition(Composition&& other)
    :  ScalarFunctionBase(std::move(other)),
      m_f(std::move(other.m_f)), m_g(std::move(other.m_g))
  {}
}

