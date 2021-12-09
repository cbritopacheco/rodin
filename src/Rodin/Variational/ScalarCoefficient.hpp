/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SCALARCOEFFICIENT_HPP
#define RODIN_VARIATIONAL_SCALARCOEFFICIENT_HPP

#include "ScalarCoefficient.h"

#include "FormLanguage/ScalarCoefficientSum.h"

namespace Rodin::Variational
{
   // ---- T (Arithmetic type) -----------------------------------------------
   // ------------------------------------------------------------------------
   template <class T>
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>
   ::ScalarCoefficient(const T& x)
      : m_x(x)
   {}

   template <class T>
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>
   ::ScalarCoefficient(const ScalarCoefficient& other)
      : m_x(other.m_x)
   {}

   template <class T>
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>&
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>::toggleSign()
   {
      m_x = -m_x;
      return *this;
   }

   template <class T>
   bool ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>::isEvaluated() const
   {
      return static_cast<bool>(m_coeff);
   }

   template <class T>
   void
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>::buildMFEMCoefficient()
   {
      track(new mfem::ConstantCoefficient(m_x));
   }

   template <class T>
   mfem::Coefficient&
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>::getMFEMCoefficient()
   {
      assert(m_coeff);
      return *m_coeff;
   }

   template <class T>
   void
   ScalarCoefficient<T,
      std::enable_if_t<std::is_arithmetic_v<T>>>::track(mfem::Coefficient* ptr)
   {
      m_coeff = std::unique_ptr<mfem::Coefficient>(ptr);
   }

   template <class T>
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>*
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>::copy()
   const noexcept
   {
      return new ScalarCoefficient(*this);
   }

   template <class T>
   template <class ... Args>
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>*
   ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>::create(Args&&... args)
   noexcept
   {
      return new ScalarCoefficient(std::forward<Args>(args)...);
   }

   // ---- ScalarCoefficient<T> ----------------------------------------------------------
   // ------------------------------------------------------------------------
   template <class T>
   ScalarCoefficient<ScalarCoefficient<T>>
   ::ScalarCoefficient(const ScalarCoefficient<T>& nested)
      : m_nested(nested.copy())
   {}

   template <class T>
   ScalarCoefficient<ScalarCoefficient<T>>
   ::ScalarCoefficient(const ScalarCoefficient& other)
      : m_nested(other.m_nested->copy())
   {}

   template <class T>
   void ScalarCoefficient<ScalarCoefficient<T>>::buildMFEMCoefficient()
   {
      m_nested->eval();
   }

   template <class T>
   mfem::Coefficient&
   ScalarCoefficient<ScalarCoefficient<T>>::getMFEMCoefficient()
   {
      return m_nested->coeff();
   }

   template <class T>
   ScalarCoefficient<ScalarCoefficient<T>>&
   ScalarCoefficient<ScalarCoefficient<T>>::toggleSign()
   {
      m_nested->toggleSign();
      return *this;
   }

   template <class T>
   bool ScalarCoefficient<ScalarCoefficient<T>>::isEvaluated() const
   {
      return m_nested->isEvaluated();
   }

   template <class T>
   template <class ... Args>
   ScalarCoefficient<ScalarCoefficient<T>>*
   ScalarCoefficient<ScalarCoefficient<T>>::create(Args&&... args) noexcept
   {
      return new ScalarCoefficient(std::forward<Args>(args)...);
   }

   template <class T>
   ScalarCoefficient<ScalarCoefficient<T>>* ScalarCoefficient<ScalarCoefficient<T>>::copy() const noexcept
   {
      return new ScalarCoefficient(*this);
   }

   // ---- FormLanguage::CoeffSum<Lhs, Rhs> ----------------------------------
   // ------------------------------------------------------------------------
   template <class Lhs, class Rhs>
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>
   ::ScalarCoefficient(const FormLanguage::ScalarCoefficientSum<Lhs, Rhs>& expr)
      : m_expr(expr.copy())
   {}

   template <class Lhs, class Rhs>
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>
   ::ScalarCoefficient(const ScalarCoefficient& other)
      :  m_expr(other.m_expr->copy())
   {}

   template <class Lhs, class Rhs>
   bool ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>::isEvaluated() const
   {
      return static_cast<bool>(m_coeff);
   }

   template <class Lhs, class Rhs>
   void
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>::buildMFEMCoefficient()
   {
      m_lhs = std::make_unique<ScalarCoefficient<Lhs>>(m_expr->lhs());
      m_rhs = std::make_unique<ScalarCoefficient<Rhs>>(m_expr->rhs());

      m_lhs->eval();
      m_rhs->eval();

      track(new mfem::SumCoefficient(m_lhs->coeff(), m_rhs->coeff()));
   }

   template <class Lhs, class Rhs>
   mfem::Coefficient&
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>::getMFEMCoefficient()
   {
      assert(m_coeff);
      return *m_coeff;
   }

   template <class Lhs, class Rhs>
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>&
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>::toggleSign()
   {
      m_lhs->toggleSign();
      m_rhs->toggleSign();
      return *this;
   }

   template <class Lhs, class Rhs>
   void
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>::track(mfem::Coefficient* ptr)
   {
      m_coeff = std::unique_ptr<mfem::Coefficient>(ptr);
   }

   template <class Lhs, class Rhs>
   template <class ... Args>
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>*
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>::create(Args&&... args)
   noexcept
   {
      return new ScalarCoefficient(std::forward<Args>(args)...);
   }

   template <class Lhs, class Rhs>
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>*
   ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>::copy()
   const noexcept
   {
      return new ScalarCoefficient(*this);
   }

   // ---- CoeffUnaryMinus<T> ------------------------------------------------
   // ------------------------------------------------------------------------
   template <class T>
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>
   ::ScalarCoefficient(const FormLanguage::ScalarCoefficientUnaryMinus<T>& expr)
      : m_expr(expr.copy())
   {}

   template <class T>
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>
   ::ScalarCoefficient(const ScalarCoefficient& other)
      : m_expr(other.m_expr->copy())
   {}

   template <class T>
   void
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>::buildMFEMCoefficient()
   {
      m_expr->toggleSign().eval();
   }

   template <class T>
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>&
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>::toggleSign()
   {
      m_expr->toggleSign();
      return *this;
   }

   template <class T>
   mfem::Coefficient&
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>::getMFEMCoefficient()
   {
      return m_expr->coeff();
   }

   template <class T>
   bool
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>::isEvaluated()
   const
   {
      return m_expr->isEvaluated();
   }

   template <class T>
   template <class ... Args>
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>*
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>::create(Args&&... args) noexcept
   {
      return new ScalarCoefficient(std::forward<Args>(args)...);
   }

   template <class T>
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>*
   ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>::copy()
   const noexcept
   {
      return new ScalarCoefficient(*this);
   }
}

#endif
