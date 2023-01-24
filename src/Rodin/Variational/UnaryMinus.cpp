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
  // ---- FunctionBase ------------------------------------------------------

  UnaryMinus<FunctionBase>::UnaryMinus(const FunctionBase& op)
    :  FunctionBase(op),
      m_op(op.copy())
  {}

  UnaryMinus<FunctionBase>::UnaryMinus(const UnaryMinus& other)
    :  FunctionBase(other),
      m_op(other.m_op->copy())
  {}

  UnaryMinus<FunctionBase>::UnaryMinus(UnaryMinus&& other)
    : FunctionBase(std::move(other)),
      m_op(std::move(other.m_op))
  {}

  RangeShape UnaryMinus<FunctionBase>::getRangeShape() const
  {
    return m_op->getRangeShape();
  }

  UnaryMinus<FunctionBase> operator-(const FunctionBase& op)
  {
    return UnaryMinus(op);
  }

  // ---- LinearFormIntegratorBase ------------------------------------------
  UnaryMinus<LinearFormIntegratorBase>::UnaryMinus(const LinearFormIntegratorBase& op)
    :  LinearFormIntegratorBase(op),
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

  mfem::Vector
  UnaryMinus<LinearFormIntegratorBase>::getVector(const Geometry::Simplex& element)
  const
  {
    mfem::Vector vec = m_op->getVector(element);
    vec *= -1.0;
    return vec;
  }

  UnaryMinus<LinearFormIntegratorBase>
  operator-(const LinearFormIntegratorBase& op)
  {
    return UnaryMinus(op);
  }

  // ---- BilinearFormIntegratorBase ----------------------------------------
  UnaryMinus<BilinearFormIntegratorBase>::UnaryMinus(const BilinearFormIntegratorBase& op)
    :  BilinearFormIntegratorBase(op),
      m_op(op.copy())
  {}

  UnaryMinus<BilinearFormIntegratorBase>::UnaryMinus(const UnaryMinus& other)
    :  BilinearFormIntegratorBase(other),
      m_op(other.m_op->copy())
  {}

  UnaryMinus<BilinearFormIntegratorBase>::UnaryMinus(UnaryMinus&& other)
    : BilinearFormIntegratorBase(std::move(other)),
      m_op(std::move(other.m_op))
  {}

  Integrator::Region
  UnaryMinus<BilinearFormIntegratorBase>::getRegion() const
  {
    return m_op->getRegion();
  }

  mfem::DenseMatrix
  UnaryMinus<BilinearFormIntegratorBase>
  ::getMatrix(const Geometry::Simplex& element) const
  {
    mfem::DenseMatrix op = m_op->getMatrix(element);
    op.Neg();
    return op;
  }

  UnaryMinus<BilinearFormIntegratorBase> operator-(const BilinearFormIntegratorBase& op)
  {
    return UnaryMinus(op);
  }

  UnaryMinus<FormLanguage::List<BilinearFormIntegratorBase>>
  operator-(const FormLanguage::List<BilinearFormIntegratorBase>& op)
  {
    return UnaryMinus<FormLanguage::List<BilinearFormIntegratorBase>>(op);
  }

  UnaryMinus<FormLanguage::List<LinearFormIntegratorBase>>
  operator-(const FormLanguage::List<LinearFormIntegratorBase>& op)
  {
    return UnaryMinus<FormLanguage::List<LinearFormIntegratorBase>>(op);
  }
}
