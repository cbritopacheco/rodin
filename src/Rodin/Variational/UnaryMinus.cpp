/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "UnaryMinus.h"
#include "RangeShape.h"

namespace Rodin::Variational
{
  // ---- LinearFormIntegratorBase ------------------------------------------
  UnaryMinus<LinearFormIntegratorBase>::UnaryMinus(const LinearFormIntegratorBase& op)
    : LinearFormIntegratorBase(op),
      m_op(op.copy())
  {}

  UnaryMinus<LinearFormIntegratorBase>::UnaryMinus(const UnaryMinus& other)
    :  LinearFormIntegratorBase(other),
      m_op(other.m_op->copy())
  {}

  UnaryMinus<LinearFormIntegratorBase>::UnaryMinus(UnaryMinus&& other)
    : LinearFormIntegratorBase(std::move(other)),
      m_op(std::move(other.m_op))
  {}

  Integrator::Region UnaryMinus<LinearFormIntegratorBase>::getRegion() const
  {
    return m_op->getRegion();
  }

  void
  UnaryMinus<LinearFormIntegratorBase>::assemble(const Geometry::Polytope& polytope)
  {
    m_op->assemble(polytope);
    auto& res = getVector();
    res = std::move(m_op->getVector());
    res *= -1.0;
  }

  UnaryMinus<LinearFormIntegratorBase>
  operator-(const LinearFormIntegratorBase& op)
  {
    return UnaryMinus(op);
  }

  // ---- BilinearFormIntegratorBase ----------------------------------------
  UnaryMinus<LocalBilinearFormIntegratorBase>::UnaryMinus(const LocalBilinearFormIntegratorBase& op)
    : LocalBilinearFormIntegratorBase(op),
      m_op(op.copy())
  {}

  UnaryMinus<LocalBilinearFormIntegratorBase>::UnaryMinus(const UnaryMinus& other)
    : LocalBilinearFormIntegratorBase(other),
      m_op(other.m_op->copy())
  {}

  UnaryMinus<LocalBilinearFormIntegratorBase>::UnaryMinus(UnaryMinus&& other)
    : LocalBilinearFormIntegratorBase(std::move(other)),
      m_op(std::move(other.m_op))
  {}

  Integrator::Region
  UnaryMinus<LocalBilinearFormIntegratorBase>::getRegion() const
  {
    return m_op->getRegion();
  }

  void UnaryMinus<LocalBilinearFormIntegratorBase>::assemble(const Geometry::Polytope& element)
  {
    m_op->assemble(element);
    auto& res = getMatrix();
    res = std::move(m_op->getMatrix());
    res *= -1.0;
  }

  UnaryMinus<LocalBilinearFormIntegratorBase> operator-(const LocalBilinearFormIntegratorBase& op)
  {
    return UnaryMinus(op);
  }

  UnaryMinus<FormLanguage::List<LocalBilinearFormIntegratorBase>>
  operator-(const FormLanguage::List<LocalBilinearFormIntegratorBase>& op)
  {
    return UnaryMinus<FormLanguage::List<LocalBilinearFormIntegratorBase>>(op);
  }

  UnaryMinus<FormLanguage::List<LinearFormIntegratorBase>>
  operator-(const FormLanguage::List<LinearFormIntegratorBase>& op)
  {
    return UnaryMinus<FormLanguage::List<LinearFormIntegratorBase>>(op);
  }
}
