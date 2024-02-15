/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Mult.h"

namespace Rodin::Variational
{
  Mult<Scalar, LocalBilinearFormIntegratorBase>::Mult(Scalar lhs, const LocalBilinearFormIntegratorBase& rhs)
    : LocalBilinearFormIntegratorBase(rhs),
      m_lhs(lhs), m_rhs(rhs.copy())
  {}

  Mult<Scalar, LocalBilinearFormIntegratorBase>::Mult(const Mult& other)
    : LocalBilinearFormIntegratorBase(other),
      m_lhs(other.m_lhs), m_rhs(other.m_rhs->copy())
  {}

  Mult<Scalar, LocalBilinearFormIntegratorBase>::Mult(Mult&& other)
    : LocalBilinearFormIntegratorBase(std::move(other)),
      m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
  {}

  Integrator::Region
  Mult<Scalar, LocalBilinearFormIntegratorBase>::getRegion() const
  {
    return m_rhs->getRegion();
  }

  void Mult<Scalar, LocalBilinearFormIntegratorBase>::assemble(const Geometry::Polytope& element)
  {
    m_rhs->assemble(element);
    auto& res = getMatrix();
    res = std::move(m_rhs->getMatrix());
    res *= m_lhs;
  }

  Mult<Scalar, LocalBilinearFormIntegratorBase> operator*(
      Scalar lhs, const LocalBilinearFormIntegratorBase& rhs)
  {
    return Mult(lhs, rhs);
  }

  Mult<Scalar, LinearFormIntegratorBase>::Mult(Scalar lhs, const LinearFormIntegratorBase& rhs)
    : LinearFormIntegratorBase(rhs),
      m_lhs(lhs), m_rhs(rhs.copy())
  {}

  Mult<Scalar, LinearFormIntegratorBase>::Mult(const Mult& other)
    : LinearFormIntegratorBase(other),
      m_lhs(other.m_lhs), m_rhs(other.m_rhs->copy())
  {}

  Mult<Scalar, LinearFormIntegratorBase>::Mult(Mult&& other)
    : LinearFormIntegratorBase(std::move(other)),
      m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
  {}

  Integrator::Region
  Mult<Scalar, LinearFormIntegratorBase>::getRegion() const
  {
    return m_rhs->getRegion();
  }

  void Mult<Scalar, LinearFormIntegratorBase>::assemble(const Geometry::Polytope& element)
  {
    m_rhs->assemble(element);
    auto& res = getVector();
    res = std::move(m_rhs->getVector());
    res *= m_lhs;
  }

  Mult<Scalar, LinearFormIntegratorBase> operator*(
      Scalar lhs, const LinearFormIntegratorBase& rhs)
  {
    return Mult(lhs, rhs);
  }
}
