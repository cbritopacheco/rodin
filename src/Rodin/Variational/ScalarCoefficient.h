/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_COEFF_H
#define RODIN_VARIATIONAL_COEFF_H

#include <memory>
#include <type_traits>

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "FormLanguage/RodinBase.h"

namespace Rodin::Variational
{
   class ScalarCoefficientBase : public FormLanguage::RodinBase
   {
      public:
         virtual ScalarCoefficientBase& toggleSign() = 0;
         virtual mfem::Coefficient& coeff() = 0;
   };

   template <class T, class Enable>
   class ScalarCoefficient
   {
      public:
         ScalarCoefficient(T&)
         {}

         ScalarCoefficient(const T&)
         {}
   };

   /**
    * Represents a coefficient of arithmetic type `T`.
    *
    * @see [std::is_arithmetic](https://en.cppreference.com/w/cpp/types/is_arithmetic)
    */
   template <class T>
   class ScalarCoefficient<T, std::enable_if_t<std::is_arithmetic_v<T>>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const T& x);
         ScalarCoefficient(const ScalarCoefficient& other);

         void eval() override;

         bool isEvaluated() const;
         mfem::Coefficient& coeff() override;
         ScalarCoefficient& toggleSign() override;

         template <class ... Args>
         static ScalarCoefficient* create(Args&&... args) noexcept;
         virtual ScalarCoefficient* copy() const noexcept override;

      protected:
         void track(mfem::Coefficient* ptr);

      private:
         T m_x;
         std::unique_ptr<mfem::Coefficient> m_coeff;
   };

   /**
    * Represents a coefficient which can be constructed from derived types of
    * `mfem::Coefficient`.
    *
    * @see [mfem::Coefficient](https://mfem.github.io/doxygen/html/classmfem_1_1Coefficient.html)
    */
   template <class T>
   class ScalarCoefficient<T,
         std::enable_if_t<std::is_base_of_v<mfem::Coefficient, T>>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const T& coeff);
         ScalarCoefficient(const ScalarCoefficient& other);

         void eval() override;
         bool isEvaluated() const;
         mfem::Coefficient& coeff() override;

         ScalarCoefficient& toggleSign() override;

         template <class ... Args>
         static ScalarCoefficient* create(Args&&... args) noexcept;

         virtual ScalarCoefficient* copy() const noexcept override;

      protected:
         void track(mfem::Coefficient* ptr);

      private:
         bool m_sign;
         std::unique_ptr<mfem::Coefficient> m_copy;
         std::unique_ptr<mfem::Coefficient> m_coeff;
   };

   /**
    * Represents the construction from another `Coeff`. The evaluation
    * behaviour is to the copy the nested coefficient and return its
    * evaluation.
    *
    */
   template <class T>
   class ScalarCoefficient<ScalarCoefficient<T>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const ScalarCoefficient<T>& nested);
         ScalarCoefficient(const ScalarCoefficient& other);
         ScalarCoefficient(ScalarCoefficient&&) = default;

         /**
          * Evaluates the coefficient into an `mfem::Coefficient`.
          *
          * @returns The evaluated mfem::Coefficient reference.
          */
         void eval() override;
         bool isEvaluated() const;
         mfem::Coefficient& coeff() override;
         ScalarCoefficient& toggleSign() override;

         template <class ... Args>
         static ScalarCoefficient* create(Args&&... args) noexcept;

         virtual ScalarCoefficient* copy() const noexcept override;

      private:
         std::unique_ptr<ScalarCoefficient<T>> m_nested;
   };

   template <class Lhs, class Rhs>
   class ScalarCoefficient<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const FormLanguage::ScalarCoefficientSum<Lhs, Rhs>& expr);
         ScalarCoefficient(const ScalarCoefficient& other);
         ScalarCoefficient(ScalarCoefficient&&) = default;

         void eval() override;
         mfem::Coefficient& coeff() override;

         ScalarCoefficient& toggleSign() override;

         bool isEvaluated() const;

         template <class ... Args>
         static ScalarCoefficient* create(Args&&... args) noexcept;

         virtual ScalarCoefficient* copy() const noexcept override;

      protected:
         void track(mfem::Coefficient* ptr);

      private:
         std::unique_ptr<FormLanguage::ScalarCoefficientSum<Lhs, Rhs>> m_expr;

         std::unique_ptr<ScalarCoefficient<Lhs>> m_lhs;
         std::unique_ptr<ScalarCoefficient<Rhs>> m_rhs;

         std::unique_ptr<mfem::Coefficient> m_coeff;
   };

   template <class T>
   class ScalarCoefficient<FormLanguage::ScalarCoefficientUnaryMinus<T>>
      : public ScalarCoefficientBase
   {
      public:
         ScalarCoefficient(const FormLanguage::ScalarCoefficientUnaryMinus<T>& expr);
         ScalarCoefficient(const ScalarCoefficient& other);
         ScalarCoefficient(ScalarCoefficient&&) = default;

         void eval() override;
         mfem::Coefficient& coeff() override;

         ScalarCoefficient& toggleSign() override;

         bool isEvaluated() const;

         template <class ... Args>
         static ScalarCoefficient* create(Args&&... args) noexcept;

         virtual ScalarCoefficient* copy() const noexcept override;

      private:
         std::unique_ptr<FormLanguage::ScalarCoefficientUnaryMinus<T>> m_expr;

         std::unique_ptr<ScalarCoefficient<T>> m_v;

         std::unique_ptr<mfem::Coefficient> m_coeff;
   };
}

#include "ScalarCoefficient.hpp"

#endif
